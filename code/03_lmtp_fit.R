#!/usr/bin/env Rscript
# Fit the longitudinal modified treatment policy and report positivity diagnostics.
#
# Runs in GATED STAGES, and the gating is the point:
#
#   0  setup      load the frame, build the shift functions, assert the never-clamped
#                 property before any model touches the data
#   1  smoke      one cheap fit (SL.glm, folds = 2). Does lmtp accept this frame at all,
#                 including the zero-variance post-window exposure columns?
#   2  gate       natural course + the primary policy, full learner library. Reports EVERY
#                 diagnostic and NO effect estimate. Stops.
#   3  expand     the full delta ladder x {S1, S2} x {SDR, TMLE}. Separate invocation,
#                 run only once stage 2 has been read by a human.
#
# DIAGNOSTICS BEFORE ESTIMATES. Once an effect estimate has been seen, every subsequent
# decision about trimming or covariates is contaminated. config/lmtp_design.json already
# states that the bicarbonate comparison is "a diagnostic, not a covariate-selection rule";
# this script enforces the same discipline structurally rather than by good intentions.
#
# Every protocol value is READ from config/lmtp_design.json. Nothing here decides one.
#
# Usage:
#   Rscript code/03_lmtp_fit.R smoke
#   Rscript code/03_lmtp_fit.R gate
#   Rscript code/03_lmtp_fit.R expand
#
# DATA SAFETY: reads patient-level data and produces fitted objects whose influence
# functions are per-observation. Fits stay in output/intermediate_phi/. Only aggregates
# reach output/final_no_phi/.

suppressPackageStartupMessages({
  library(lmtp)
  library(nanoparquet)
  library(jsonlite)
  library(SuperLearner)
})

STAGE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(STAGE)) STAGE <- "smoke"
stopifnot(STAGE %in% c("smoke", "gate", "expand"))

REPO <- normalizePath(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = FALSE)
if (!dir.exists(file.path(REPO, "config"))) REPO <- normalizePath(".")

# ---------------------------------------------------------------------------
# CONFIG BLOCK  —  read, never decide
# ---------------------------------------------------------------------------
config <- fromJSON(file.path(REPO, "config", "config.json"))
design <- fromJSON(file.path(REPO, "config", "lmtp_design.json"), simplifyVector = TRUE)

EST      <- design$estimation
POL      <- design$policy
TIME     <- design$time
COV      <- design$covariates

DELTA_PRIMARY <- POL$delta_primary
DELTA_LADDER  <- POL$delta_ladder
FLOOR         <- POL$floor
GRID          <- TIME$outcome_grid_days
TAU           <- length(GRID)
EXPO_PERIODS  <- TIME$exposure_periods
K             <- if (identical(EST$k, "Inf")) Inf else as.numeric(EST$k)
SEED          <- EST$seed

cat(sprintf("stage: %s\n", STAGE))
cat(sprintf("protocol: tau=%d over days %s; policy intervenes in periods %s\n",
            TAU, paste(GRID, collapse = "/"), paste(EXPO_PERIODS, collapse = ",")))
cat(sprintf("  delta primary %s, ladder {%s}, floor %s, k=%s, folds=%d, seed=%d\n",
            DELTA_PRIMARY, paste(DELTA_LADDER, collapse = ", "), FLOOR,
            format(K), EST$folds, SEED))

# ---------------------------------------------------------------------------
# Column layout. Derived from the config, not written down.
# ---------------------------------------------------------------------------
A_S1 <- paste0("a", seq_len(TAU), "_s1")
A_S2 <- paste0("a", seq_len(TAU), "_s2")
Y    <- paste0("y_d", GRID)
CENS <- paste0("c_d", GRID)
COMP <- paste0("d_d", GRID)
EXPO_COLS_S1 <- A_S1[EXPO_PERIODS]
EXPO_COLS_S2 <- A_S2[EXPO_PERIODS]

# ---------------------------------------------------------------------------
# Stage 0: load and assert
# ---------------------------------------------------------------------------
stage_0_load <- function() {
  f <- file.path(REPO, "output", "intermediate_phi", "lmtp_df.parquet")
  d <- as.data.frame(read_parquet(f))
  cat(sprintf("\n  frame: %s rows x %d cols\n", format(nrow(d), big.mark = ","), ncol(d)))

  # lmtp refuses a data.table and requires no NA in trt/outcome/cens/compete.
  stopifnot(!inherits(d, "data.table"))
  need <- c(A_S1, A_S2, Y, CENS, COMP, EST$id)
  missing <- setdiff(need, names(d))
  if (length(missing)) stop("frame is missing columns: ", paste(missing, collapse = ", "))
  stopifnot(!anyNA(d[, c(A_S1, A_S2, Y, CENS, COMP)]))

  # trt must be length 1 or exactly tau (lmtp R/checks.R:29-40), and tau for a survival
  # outcome is the number of outcome columns (R/Task.R:125-130). Assert it here so the
  # failure is legible rather than an assertion from inside the package.
  stopifnot(length(A_S1) == length(Y), length(A_S2) == length(Y))
  cat(sprintf("  trt length %d == outcome length %d, tau = %d\n",
              length(A_S1), length(Y), TAU))
  d
}

# ---------------------------------------------------------------------------
# The shift function
#
# Self-limiting and NEVER CLAMPED. `max(a - delta, floor)` would map every subject in
# [floor, floor + delta) onto the single value floor, collapsing an interval to a point:
# no inverse, so no change-of-variables, so no post-intervention density (Diaz eq. 6,
# p.222), and a point mass where the observed data has only continuous density, which
# makes the density ratio infinite rather than merely large.
#
# Node-conditional: `shift` is applied one column at a time and receives the column NAME
# (lmtp R/shift.R:39-45), so the policy can act on the exposure periods and return the
# observed value everywhere else. Unintervened periods then carry a density ratio of
# exactly 1 and cost nothing.
# ---------------------------------------------------------------------------
make_shift <- function(delta, floor_, exposure_cols) {
  force(delta); force(floor_); force(exposure_cols)
  # The shift FORM is protocol. Read it and refuse anything this function does not
  # implement, rather than silently applying the self-limiting form under a config that
  # now says something else.
  if (!identical(POL$shift_form, "self_limiting_never_clamped")) {
    stop("policy.shift_form is '", POL$shift_form, "'; this script implements only ",
         "'self_limiting_never_clamped'. Implement it or revert the config.")
  }
  function(data, trt) {
    a <- data[[trt]]
    if (!trt %in% exposure_cols) return(a)
    ifelse(a - delta >= floor_, a - delta, a)
  }
}

assert_shift_is_sane <- function(d, arm_cols, expo_cols) {
  for (delta in DELTA_LADDER) {
    sh <- make_shift(delta, FLOOR, expo_cols)
    for (cl in arm_cols) {
      a <- d[[cl]]; s <- sh(d, cl)
      if (!cl %in% expo_cols) {
        stopifnot(identical(a, s))                       # policy silent outside the window
        next
      }
      # every value is either shifted by exactly delta, or untouched. Never in between,
      # and never parked on the floor.
      moved <- abs(s - (a - delta)) < 1e-9
      kept  <- abs(s - a) < 1e-9
      stopifnot(all(moved | kept))
      stopifnot(all(!(moved & (a - delta < FLOOR))))     # never shifts below the floor
      stopifnot(all(!(kept & (a - delta >= FLOOR) & (a > 0))))  # never declines a legal shift
    }
  }
  cat(sprintf("  shift asserted for %d deltas x %d columns: self-limiting, never clamped\n",
              length(DELTA_LADDER), length(arm_cols)))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Covariate layout: baseline W, and time_vary as one entry per period.
#
# Periods beyond the exposure window introduce NO new covariates, which lmtp permits
# (R/Vars.R:18 allows NULL entries). With k = Inf those periods still condition on the
# node-3 covariates through the history. THAT COUPLING IS NOT OPTIONAL: under a truncated
# k they would condition on nothing.
# ---------------------------------------------------------------------------
build_covariates <- function(d) {
  suffix_for <- function(node) if (node == 1) "0" else as.character(node)
  tv_names <- COV$time_varying_Lt
  present <- function(v) v[v %in% names(d)]

  L <- vector("list", TAU)
  for (p in seq_len(TAU)) {
    if (p <= length(EXPO_PERIODS)) {
      L[[p]] <- present(paste0(tv_names, "_", suffix_for(p)))
    } else {
      # `L[[p]] <- NULL` would DELETE the element and shrink the list, which is the
      # classic R trap. `L[p] <- list(NULL)` stores NULL as the value, which is what
      # lmtp wants: "this period introduces no new covariates" (R/Vars.R:18 permits it).
      L[p] <- list(NULL)
    }
  }
  # PARTITION, not union. `covariates.baseline_L0` mixes two kinds of term and lmtp wants
  # them in different arguments:
  #
  #   truly time-invariant (age, sex, weight, comorbidity)  -> `baseline`, adjusted for at
  #                                                            every time point
  #   measured over [-24h, 0), i.e. the "_0" columns        -> `time_vary[[1]]`, the first
  #                                                            time-varying set, measured
  #                                                            before A_1
  #
  # Putting a "_0" column in BOTH is not merely redundant: the same name then appears twice
  # in the constructed history and lmtp's variable renaming produces a predictor set that
  # its own predict step cannot satisfy ("undefined columns selected"). Each declared name
  # must land in exactly one place.
  cci <- COV$definitions$cci_components$entered
  bl  <- c(COV$baseline_L0, cci)
  static_terms <- bl[bl %in% names(d)]                       # resolve bare -> baseline
  node0_terms  <- bl[!(bl %in% names(d)) & paste0(bl, "_0") %in% names(d)]
  unresolved   <- setdiff(bl, c(static_terms, node0_terms))
  W <- static_terms
  L[[1]] <- unique(c(L[[1]], paste0(node0_terms, "_0")))

  stopifnot(length(L) == TAU)

  # EVERY COVARIATE THE CONFIG DECLARES MUST RESOLVE TO A COLUMN. `present()` above
  # silently drops names the frame does not carry, which is exactly how a config key
  # becomes decorative: 02 stops building a covariate, 03 quietly fits without it, and
  # nothing fails. Fail loudly instead.
  want_tv <- unlist(lapply(seq_along(EXPO_PERIODS),
                           function(p) paste0(tv_names, "_", suffix_for(p))))
  missing <- c(setdiff(want_tv, names(d)), unresolved)
  if (length(missing)) {
    stop("config declares covariates the frame does not carry, so they would be ",
         "silently dropped: ", paste(missing, collapse = ", "))
  }
  cat(sprintf("  baseline W: %d terms | time_vary: %s\n", length(W),
              paste(sapply(L, function(x) if (is.null(x)) "NULL" else length(x)),
                    collapse = ", ")))
  cat(sprintf("  config coverage: all %d declared covariate columns present\n",
              length(want_tv) + length(bl)))
  list(W = W, L = L)
}

# ---------------------------------------------------------------------------
# Imputation and encoding
#
# 02 deliberately leaves BASELINE covariates missing, on the grounds that imputation
# belongs with the fit rather than with dataset construction. 37.2% of rows carry at
# least one missing covariate, concentrated in P/F, S/F, lactate and IMV status.
#
# For the pilot: median for continuous, mode for binary, plus an EXPLICIT missingness
# indicator for every imputed column. The indicator is not bookkeeping -- it is
# informative. Arterial sampling tracks severity, so "was a blood gas drawn" carries
# signal of its own, and SuperLearner can use it.
#
# The honest cost, stated rather than hidden: imputation is treated as fixed
# preprocessing, so its uncertainty does not propagate into the influence-function
# variance and standard errors are mildly optimistic. fluid_ARDS S4.4 specifies MICE;
# that choice is deferred to the definitive analysis and must be pre-specified there.
# ---------------------------------------------------------------------------
impute_and_encode <- function(d, vars) {
  cols <- unique(c(vars$W, unlist(vars$L)))
  added <- character(0)

  # sex is a character column; SuperLearner needs numerics.
  if ("sex" %in% cols && is.character(d$sex)) {
    d$sex <- as.integer(tolower(d$sex) == "female")
    cat("  encoded sex as binary (1 = female)\n")
  }

  # An indicator is only worth carrying if the missingness is common enough to survive
  # cross-fitting. Below the threshold it is near-constant overall and EXACTLY constant
  # inside some training fold, at which point SuperLearner drops the variable and its
  # predict step asks for a column that no longer exists ("undefined columns selected").
  MV <- COV$missing_values
  min_n <- max(MV$indicator_min_count, ceiling(MV$indicator_min_pct / 100 * nrow(d)))

  n_imputed <- 0L; skipped <- character(0)
  for (cl in cols) {
    x <- d[[cl]]
    if (!is.numeric(x)) x <- as.numeric(as.factor(x))
    na <- is.na(x)
    if (!any(na)) { d[[cl]] <- x; next }
    binary <- all(x[!na] %in% c(0, 1))
    fill <- if (binary) as.numeric(names(which.max(table(x[!na])))) else median(x, na.rm = TRUE)
    if (sum(na) >= min_n) {
      ind <- paste0(cl, "_miss")
      d[[ind]] <- as.integer(na)
      added <- c(added, ind)
    } else {
      skipped <- c(skipped, sprintf("%s (%d)", cl, sum(na)))
    }
    x[na] <- fill
    d[[cl]] <- x
    n_imputed <- n_imputed + sum(na)
  }
  if (length(skipped)) {
    cat(sprintf("  no indicator for %d columns below the %d-row threshold: %s\n",
                length(skipped), min_n, paste(skipped, collapse = ", ")))
  }
  # The indicators join the same layer as the covariate they describe.
  for (ind in added) {
    base <- sub("_miss$", "", ind)
    if (base %in% vars$W) vars$W <- c(vars$W, ind)
    for (p in seq_along(vars$L)) {
      if (!is.null(vars$L[[p]]) && base %in% vars$L[[p]]) vars$L[[p]] <- c(vars$L[[p]], ind)
    }
  }
  cat(sprintf("  imputed %s cells; added %d missingness indicators\n",
              format(n_imputed, big.mark = ","), length(added)))

  # Any covariate with no variance left carries no information and will be dropped by
  # some learners and not others, which is exactly the silent inconsistency that breaks
  # cross-fitting. Drop them here, loudly, rather than letting each learner decide.
  drop_constant <- function(v) {
    if (is.null(v)) return(v)
    keep <- v[sapply(v, function(cl) var(d[[cl]], na.rm = TRUE) > 0)]
    dropped <- setdiff(v, keep)
    if (length(dropped)) cat("  dropped zero-variance covariates:",
                            paste(dropped, collapse = ", "), "\n")
    keep
  }
  vars$W <- drop_constant(vars$W)
  vars$L <- lapply(vars$L, drop_constant)

  stopifnot(!anyNA(d[unique(c(vars$W, unlist(vars$L)))]))
  list(d = d, vars = vars)
}

# ---------------------------------------------------------------------------
# Diagnostics. These are the deliverable of the pilot, not a footnote to it.
# ---------------------------------------------------------------------------
kish_ess <- function(w) sum(w)^2 / sum(w^2)

diagnose <- function(fit, label, d, arm_cols) {
  dr <- fit$density_ratios
  stopifnot(!is.null(dr))
  cum <- t(apply(dr, 1, cumprod))
  rows <- list()
  add <- function(metric, period, value)
    rows[[length(rows) + 1]] <<- data.frame(policy = label, metric = metric,
                                            period = period, value = as.numeric(value))

  for (p in seq_len(ncol(dr))) {
    q <- quantile(dr[, p], c(.5, .95, .99, 1), na.rm = TRUE)
    add("dr_median", p, q[1]); add("dr_p95", p, q[2])
    add("dr_p99", p, q[3]);    add("dr_max", p, q[4])
    cq <- quantile(cum[, p], c(.5, .95, .999, 1), na.rm = TRUE)
    add("cum_median", p, cq[1]); add("cum_p95", p, cq[2])
    add("cum_p999", p, cq[3]);   add("cum_max", p, cq[4])
    w <- cum[, p]
    add("ess_kish", p, kish_ess(w))
    add("ess_pct_of_n", p, 100 * kish_ess(w) / length(w))
    top <- sort(w, decreasing = TRUE)[seq_len(max(1, floor(0.01 * length(w))))]
    add("weight_share_top1pct", p, 100 * sum(top) / sum(w))
  }

  # THE conditional positivity check: node-2 density ratio stratified by A1, not
  # marginally. This is the risk Tier 1 structurally could not see.
  if (ncol(dr) >= 2) {
    a1 <- d[[arm_cols[1]]]
    on <- d$node_status_1 == "on_crrt"
    qs <- quantile(a1[on], seq(0, 1, .2), na.rm = TRUE)
    bin <- cut(a1, unique(qs), include.lowest = TRUE, labels = FALSE)
    for (b in sort(unique(bin[!is.na(bin)]))) {
      sel <- which(bin == b)
      add(sprintf("dr_node2_given_A1_q%d_p95", b), 2, quantile(dr[sel, 2], .95, na.rm = TRUE))
      add(sprintf("dr_node2_given_A1_q%d_max", b), 2, max(dr[sel, 2], na.rm = TRUE))
    }
  }
  do.call(rbind, rows)
}

sl_coefficients <- function(fit, label) {
  grab <- function(fits, nuisance) {
    out <- list()
    if (is.null(fits)) return(NULL)
    for (i in seq_along(fits)) {
      f <- fits[[i]]
      co <- tryCatch(f$coef, error = function(e) NULL)
      if (is.null(co)) next
      out[[length(out) + 1]] <- data.frame(policy = label, nuisance = nuisance,
                                           fold = i, learner = names(co),
                                           coefficient = as.numeric(co))
    }
    if (length(out)) do.call(rbind, out) else NULL
  }
  rbind(grab(fit$fits_m, "outcome"), grab(fit$fits_r, "treatment"))
}

# ---------------------------------------------------------------------------
# Provenance and writing, mirroring 01 and 02 stage 8.
# ---------------------------------------------------------------------------
code_version <- function() {
  v <- tryCatch(system2("git", c("-C", REPO, "describe", "--always", "--dirty"),
                        stdout = TRUE, stderr = NULL), error = function(e) "unknown")
  if (length(v) == 0 || is.na(v[1])) "unknown" else v[1]
}

write_shareable <- function(df, name) {
  prov <- list(site_id = config$site_name, code_version = code_version(),
               clif_version = config$clif_version,
               definition_version = design$definition_version,
               generated = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"))
  for (k in names(prov)) df[[k]] <- prov[[k]]
  ID_LIKE <- c("encounter_block", "hospitalization_id", "patient_id", "_dttm", "eif")
  leaked <- names(df)[sapply(names(df), function(n) any(sapply(ID_LIKE, grepl, x = n)))]
  if (length(leaked)) stop("identifier-like columns in a shareable file: ",
                           paste(leaked, collapse = ", "))
  share <- file.path(REPO, "output", "final_no_phi")
  dir.create(share, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(share, sprintf("%s_%s.csv", config$site_name, name))
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("  wrote final_no_phi/%s  %d rows, PHI-checked\n", basename(path), nrow(df)))
}

save_fits <- function(obj, name) {
  phi <- file.path(REPO, "output", "intermediate_phi")
  dir.create(phi, showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, file.path(phi, sprintf("%s.rds", name)))
  cat(sprintf("  wrote intermediate_phi/%s.rds  (per-observation EIFs: never share)\n", name))
}

# ---------------------------------------------------------------------------
# The fit
# ---------------------------------------------------------------------------
fit_policy <- function(d, arm_cols, expo_cols, vars, delta, learners, folds,
                       estimator = "sdr", trim = EST$trim) {
  sh <- if (is.null(delta)) NULL else make_shift(delta, FLOOR, expo_cols)
  # learners_trt is read separately from learners_outcome. They are identical in the
  # shipped config, but the config documents that learners_trt "must be BINARY
  # CLASSIFIERS", so the two are allowed to diverge and passing one for both would
  # silently ignore that. `learners` overrides both only for the cheap smoke test.
  lo <- if (!is.null(learners)) learners else EST$learners_outcome
  lt <- if (!is.null(learners)) learners else EST$learners_trt
  args <- list(data = d, trt = arm_cols, outcome = Y, baseline = vars$W,
               time_vary = vars$L, cens = CENS, compete = COMP,
               shift = sh, mtp = EST$mtp, k = K, outcome_type = "survival",
               id = EST$id, folds = folds,
               learners_outcome = lo, learners_trt = lt,
               control = lmtp_control(.trim = trim,
                                      .bound = EST$bound,
                                      .discrete = EST$discrete_superlearner,
                                      .return_full_fits = TRUE))
  set.seed(SEED)
  do.call(if (estimator == "tmle") lmtp_tmle else lmtp_sdr, args)
}

# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------
d <- stage_0_load()
assert_shift_is_sane(d, A_S1, EXPO_COLS_S1)
assert_shift_is_sane(d, A_S2, EXPO_COLS_S2)
vars <- build_covariates(d)
imp <- impute_and_encode(d, vars); d <- imp$d; vars <- imp$vars

if (STAGE == "smoke") {
  cat("\n=== STAGE 1: smoke test (SL.glm, folds = 2) ===\n")
  cat("Only question: does lmtp accept this frame, including the constant a4-a6 columns?\n")
  t0 <- Sys.time()
  fit <- fit_policy(d, A_S1, EXPO_COLS_S1, vars, DELTA_PRIMARY,
                    learners = "SL.glm", folds = 2)
  cat(sprintf("  completed in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  dr <- fit$density_ratios
  cat(sprintf("  density_ratios: %d x %d\n", nrow(dr), ncol(dr)))
  cat("  per-period max density ratio: ",
      paste(sprintf("%.3f", apply(dr, 2, max, na.rm = TRUE)), collapse = "  "), "\n")
  # In an unintervened period the shifted and natural treatments are identical, so the
  # density ratio is 1 -- but ONLY for subjects still at risk. The weight also carries
  # the at-risk indicator R_t = 1{D_t = 0, Y_t = 0} (Diaz p.223), so anyone who has
  # already died or been discharged gets 0. Those zeros are the competing-risk structure
  # zeroing out the recursion, not a positivity failure, and checking for "all exactly 1"
  # would wrongly flag them.
  unint <- setdiff(seq_len(TAU), EXPO_PERIODS)
  sub <- as.matrix(dr[, unint, drop = FALSE])
  nz <- sub[sub > 1e-12]
  ok <- length(nz) > 0 && all(abs(nz - 1) < 1e-6)
  cat(sprintf("  unintervened periods (%s): ratio is 1 for every at-risk subject: %s\n",
              paste(unint, collapse = ","), ok))
  cat(sprintf("    (%.0f%% of those cells are 0, i.e. already dead or discharged)\n",
              100 * mean(sub <= 1e-12)))
  cat("\nSmoke test passed. Frame is accepted. Run `gate` next.\n")
  quit(status = 0)
}

if (STAGE == "gate") {
  cat("\n=== STAGE 2: DIAGNOSTIC GATE ===\n")
  cat("Natural course and the primary policy, full learner library.\n")
  cat("Reports diagnostics ONLY. No effect estimate is printed at this stage.\n\n")
  learners <- EST$learners_outcome
  t0 <- Sys.time()

  cat("  fitting natural course (shift = NULL) ...\n")
  nat <- fit_policy(d, A_S1, EXPO_COLS_S1, vars, NULL, learners, EST$folds)
  cat("  fitting delta =", DELTA_PRIMARY, "under S1 ...\n")
  shf <- fit_policy(d, A_S1, EXPO_COLS_S1, vars, DELTA_PRIMARY, learners, EST$folds)
  cat(sprintf("  both fits completed in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  diag <- rbind(diagnose(nat, "natural_course", d, A_S1),
                diagnose(shf, sprintf("delta_%s_s1", DELTA_PRIMARY), d, A_S1))

  # Shifted fraction, per period, as a headline
  for (p in EXPO_PERIODS) {
    on <- d[[paste0("node_status_", p)]] == "on_crrt"
    a <- d[[A_S1[p]]][on]
    diag <- rbind(diag, data.frame(policy = sprintf("delta_%s_s1", DELTA_PRIMARY),
                                   metric = "shifted_fraction_pct", period = p,
                                   value = 100 * mean(a - DELTA_PRIMARY >= FLOOR)))
  }

  cat("\n---- cumulative density ratio by period (the compounding is the risk) ----\n")
  sub <- diag[diag$policy != "natural_course" &
              diag$metric %in% c("cum_median", "cum_p95", "cum_p999", "cum_max"), ]
  print(reshape(sub[, c("metric", "period", "value")], idvar = "metric",
                timevar = "period", direction = "wide"), row.names = FALSE)

  cat("\n---- effective sample size (Kish) as % of n ----\n")
  ess <- diag[diag$policy != "natural_course" & diag$metric == "ess_pct_of_n", ]
  cat("  ", paste(sprintf("p%d: %.1f%%", ess$period, ess$value), collapse = "   "), "\n")
  cat("   S4.6 scale: ~60% comfortable, ~30% the limit, ~10% means tau=3 territory\n")

  cat("\n---- weight concentration: share of total weight in the top 1% ----\n")
  wc <- diag[diag$policy != "natural_course" & diag$metric == "weight_share_top1pct", ]
  cat("  ", paste(sprintf("p%d: %.1f%%", wc$period, wc$value), collapse = "   "), "\n")
  cat("   above ~20% means a handful of patients carry the estimate\n")

  cat("\n---- CONDITIONAL positivity at node 2, by A1 quintile ----\n")
  cq <- diag[diag$policy != "natural_course" & grepl("given_A1", diag$metric), ]
  print(cq[, c("metric", "value")], row.names = FALSE)

  write_shareable(diag, "lmtp_diagnostics_gate")
  sl <- rbind(sl_coefficients(nat, "natural_course"),
              sl_coefficients(shf, sprintf("delta_%s_s1", DELTA_PRIMARY)))
  if (!is.null(sl)) write_shareable(sl, "sl_coefficients_gate")
  save_fits(list(natural = nat, shifted = shf), "lmtp_fits_gate")

  cat("\n=== GATE REACHED ===\n")
  cat("Read the diagnostics above before running `expand`. No effect estimate has been\n")
  cat("printed, deliberately: seeing one now would contaminate any later decision about\n")
  cat("trimming or covariates.\n")
  quit(status = 0)
}

if (STAGE == "expand") {
  cat("\n=== STAGE 3: full ladder ===\n")
  learners <- EST$learners_outcome
  results <- list(); diags <- list(); sls <- list()
  for (arm in c("s1", "s2")) {
    cols <- if (arm == "s1") A_S1 else A_S2
    expo <- if (arm == "s1") EXPO_COLS_S1 else EXPO_COLS_S2
    nat <- fit_policy(d, cols, expo, vars, NULL, learners, EST$folds)
    for (delta in DELTA_LADDER) {
      for (est in c("sdr", "tmle")) {
        lab <- sprintf("delta_%s_%s_%s", delta, arm, est)
        cat("  fitting", lab, "...\n")
        f <- fit_policy(d, cols, expo, vars, delta, learners, EST$folds, estimator = est)
        ct <- lmtp_contrast(f, ref = nat, type = "additive")
        results[[lab]] <- data.frame(policy = lab, arm = arm, delta = delta, estimator = est,
                                     estimate = ct$vals$theta, std_error = ct$vals$std.error,
                                     conf_low = ct$vals$conf.low, conf_high = ct$vals$conf.high)
        if (est == "sdr") diags[[lab]] <- diagnose(f, lab, d, cols)
        sls[[lab]] <- sl_coefficients(f, lab)
        results[[paste0(lab, "_fit")]] <- NULL
      }
    }
    results[[paste0("natural_", arm)]] <- data.frame(
      policy = paste0("natural_", arm), arm = arm, delta = NA, estimator = "sdr",
      estimate = tidy(nat)$estimate, std_error = tidy(nat)$std.error,
      conf_low = tidy(nat)$conf.low, conf_high = tidy(nat)$conf.high)
  }
  est_df <- do.call(rbind, results)
  write_shareable(est_df, "lmtp_estimates")
  write_shareable(do.call(rbind, diags), "lmtp_diagnostics_full")
  if (length(sls)) write_shareable(do.call(rbind, sls), "sl_coefficients_full")
  cat("\nDone.\n")
}
