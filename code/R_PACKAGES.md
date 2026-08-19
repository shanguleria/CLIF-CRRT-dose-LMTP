# R environment for the LMTP fit

**Status: installed, locked, and verified on R 4.3.1.** `renv.lock` pins 74 packages
including `lmtp` 1.5.4. The question this file used to flag as "the first real risk to
retire on the R side" is retired: `lmtp` 1.5.4 runs on R 4.3.1, and a survival fit with
competing risks completes end to end.

## What was verified, 2026-08-18

| Check | Result |
|---|---|
| `lmtp` current CRAN version | **1.5.4**, the version the design targets |
| `lmtp` declared requirement | `R (>= 2.10)` |
| Recursive dependency closure | 64 packages, **none** requiring more than R 4.3.1 |
| `lmtp` needs compilation | **No.** Pure R |
| Install on R 4.3.1 | **Succeeds** |
| Survival fit with `compete` | **Runs**, returns an estimate with CI |

## The install strategy, and why it is not the obvious one

CRAN still serves macOS arm64 binaries for R 4.3, but that repo is **frozen**: R 4.3 has
fallen out of the current-plus-previous window, so its binaries stopped advancing. The
newest `lmtp` binary for R 4.3 is **1.5.2**, not the 1.5.4 the design pins. Same story
across the stack (`glmnet` binary 4.1-9 against source 5.0, `xgboost` 1.7.11.1 against
3.2.1.1, `SuperLearner` 2.0-29 against 2.0-40).

Building everything from source is not an option here: 32 of the 63 closure packages need
compilation and **`gfortran` is not installed**, which sinks `glmnet`, `earth` and
`survival` immediately.

So the strategy is split, and it works because `lmtp` itself is pure R:

```r
# 1. binaries for the compiled stack (frozen 4.3 versions, but they build)
renv::install(c("mlr3superlearner","earth","nnls","SuperLearner","glmnet",
                "ranger","xgboost","data.table","future","progressr","isotone"),
              type = "binary", prompt = FALSE)

# 2. lmtp at the pinned version from source; no compiler needed
renv::install("lmtp@1.5.4", prompt = FALSE)

renv::snapshot(type = "all")
```

Installed in 45 seconds.

### Why 1.5.2 would not have been an acceptable fallback

Two changes between 1.5.2 and 1.5.4 bear directly on this analysis:

- **1.5.4**: risk ratios and odds ratios in `lmtp_contrast()` are now computed on the log
  scale and exponentiated. On 1.5.2 a ratio contrast and its confidence interval are
  built on the natural scale, which is the wrong scale for a ratio.
- **1.5.4**: added a check that Super Learner cross-validation folds do not exceed the
  number of clusters. This design passes `id`, so that check is live for us.

`1.5.3`'s fixes concern multivariate exposure and do not apply.

## Locked versions

R 4.3.1, 74 packages. The ones that matter:

| Package | Version | | Package | Version |
|---|---|---|---|---|
| `lmtp` | **1.5.4** | | `SuperLearner` | 2.0-29 |
| `ife` | 0.2.5 | | `glmnet` | 4.1-9 |
| `mlr3superlearner` | 0.1.2 | | `ranger` | 0.17.0 |
| `earth` | 5.3.4 | | `xgboost` | 1.7.11.1 |
| `nnls` | 1.6 | | `data.table` | 1.17.8 |
| `isotone` | 1.1-2 | | `future` / `progressr` | 1.40.0 / 0.15.1 |

`ife` is `lmtp`'s influence-function class and arrived as a dependency; 1.5.4 requires
`> 0.1.2`.

> **Federated consequence, flagged not solved.** These are the newest versions that exist
> as R 4.3 binaries, so at a site running R 4.5 they are *old* versions that `renv::restore()`
> will have to fetch from the CRAN archive, in some cases building from source. The lock is
> reproducible, which is the point; it is not necessarily convenient at a site with a newer
> R. Decide before shipping whether sites are asked to match R 4.3.1 or whether the lock is
> regenerated per major R version.

## The R version landmine, still true

| | Version |
|---|---|
| System R on this machine | **4.3.1** (2023-06-16) |
| Sibling `CLIF-epidemiology-of-CRRT` `renv.lock` targets | **4.5.2** |

`renv::restore()` against the sibling's lock **fails here**: no matching binaries for
4.3.1, so it compiles from source and dies on a batch (`arrow`, `Matrix`,
`RcppArmadillo`, `V8`, `systemfonts`), emptying `renv/library` on the way out. **Build
our own lock, which is now done. Do not copy, symlink, or restore the sibling's.**

## Two package facts that change the design

### 1. IPW and g-computation are defunct, not deprecated

`lmtp_ipw()` and `lmtp_sub()` live in `R/defunct.R` and call
`lifecycle::deprecate_stop()`. They **raise an error**, they do not warn:

> IPW requires the use of correctly pre-specified parametric models for valid statistical
> inference. Use `lmtp_tmle()` or `lmtp_sdr()` instead.

Removed in 1.5.0, so no reachable version of `lmtp` offers them.
`config/lmtp_design.json` lists both under `estimation.diagnostics`. That entry is not
implementable and needs amending.

### 2. `trt` must be length 1 or exactly `tau`, and `tau` comes from the outcome vector

Verified by running it (`lmtp` 1.5.4, synthetic data):

| Test | Result |
|---|---|
| `trt` length 3, `outcome` length 4 | **`Assertion on 'trt' failed: 'trt' should be of length 1 or 4`** |
| `trt` length 4, `outcome` length 4 | runs |
| `trt` length 1, `outcome` length 4 | runs (recycled) |

So exposure nodes at 0/24/48h against an outcome grid of day 3/7/14/30 **cannot be passed
as written**. `lmtp_survival()` confirms the intended shape: it loops over time and slices
`trt[seq_len(time)]`, `time_vary[seq_len(time)]`, `cens[seq_len(time)]`,
`compete[seq_len(time)]`, so the package expects **one common time grid** in which every
period carries an exposure, covariates, censoring, competing and outcome column.

**The fix, verified to work.** `shift` is called as `.f(data, column_name)`
(`R/shift.R:shift_trt_character`), so it can behave differently per node:

```r
EXPOSURE_NODES <- c("a1", "a2", "a3")
shift <- function(data, trt) {
  a <- data[[trt]]
  if (!trt %in% EXPOSURE_NODES) return(a)        # post-window: policy does nothing
  ifelse(a - DELTA >= FLOOR, a - DELTA, a)        # self-limiting, never clamped
}
```

Run on a 6-period synthetic grid intervening only in periods 1-3, the estimated density
ratios are **exactly 1.000 at every quantile in periods 4, 5 and 6**. Carrying
unintervened periods therefore costs nothing in positivity, which is what makes this a
real solution rather than a workaround.

## What the fitted object gives you

`class "lmtp"`, fields `estimator`, `estimate`, `shift`, `outcome_reg`,
**`density_ratios`**, `fits_m`, `fits_r`, `outcome_type`. `lmtp::tidy()` returns
`estimate`, `std.error`, `conf.low`, `conf.high`.

`density_ratios` is an **n x tau matrix**, which is the positivity diagnostic the design
asks for: the trimmed fraction, the per-node distribution, and the pre-specified node-2
comparison with and without bicarbonate all read straight off it.

## Still open

`lmtp_control()` defaults worth setting explicitly: `.trim = 0.999` and
`.discrete = TRUE`. The latter is *discrete* Super Learner, a single best learner rather
than the weighted ensemble; set `.discrete = FALSE` for true stacking.
