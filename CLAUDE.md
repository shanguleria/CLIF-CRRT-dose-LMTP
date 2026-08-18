# CLAUDE.md

Guidance for Claude Code working in this repository.

## Project overview

Apply **longitudinal modified treatment policies (LMTP)** to CRRT dose as a
**time-varying** exposure at t = 0 / 24 / 48h, with time-varying covariates,
estimating the effect of a fixed dose *reduction* at each interval on 30-day
in-hospital mortality, with discharge alive as a competing event. Multi-site,
CLIF 2.1.0, federated (no person-level data leaves a site).

This is the successor to the time-varying MSM that was cut from the CJASN
manuscript for positivity violations. It is **not** a re-run of that analysis:
the MTP framing changes the estimand so that positivity is structurally
attainable rather than something to be truncated into submission.

## Current status

`code/01_build_cohort.py` is complete and runs end to end, from raw CLIF 2.1.0
tables to written artifacts:

```
166,814 hospitalizations -> 166,677 encounter blocks -> 3,152 with CRRT
-> 2,246 after ESRD -> 2,147 with a t=0 -> 2,144 with a weight
-> 2,144 with a dose series -> 2,144 with an outcome
```

It writes `cohort.parquet` and `dose_series.parquet` to `output/intermediate_phi/`,
and a provenance-stamped STROBE table to `output/final_no_phi/`.

`02_build_lmtp_df.py` and `03_lmtp_fit.R` do not exist, and no R packages are
installed. The decision that used to block 02 was made on 2026-08-18, so **02 is
now writable**.

## The estimand decision that used to block dataset code

**How to handle CRRT discontinuation.** The exposure is undefined once CRRT
stops, and 39% of encounter blocks at the coordinating site have no third node.
**RESOLVED 2026-08-18.** `discontinuation_handling: "dose_zero_unshifted"`.

CRRT liberation is neither a competing event nor censoring. It is an exposure value
of **zero, left unshifted by the policy**.

- Not a competing event: liberation does not preclude in-hospital death, and
  `compete` is correctly occupied by discharge alive.
- Not censoring: that would estimate the effect in a world where nobody is liberated
  within 72h, and positivity fails for recovering patients.
- Safe because the self-limiting shift leaves zero unchanged (`0 - delta < floor` for
  every delta in the ladder), so the policy can never move a patient onto or off
  therapy. That is what stops this being a covert "dose or no dose" contrast.

The paired decision is the **node statistic**: a time-weighted mean over charted CRRT
records in each node, `sum(dose_h) / n_charted_hours`. Charted zero-dose intervals
enter; uncharted gaps are not imputed. One rule for every cause of a zero-dose
interval: liberation, SCUF-only, machine downtime, filter clotting. Full reasoning
and the measured comparison live in `exposure._node_statistic_why` and
`_discontinuation_handling_RESOLVED` in `config/lmtp_design.json`.

Charted zero versus uncharted gap is pre-specified as a **bracket**, not a single
choice (`exposure._node_statistic_sensitivity_LADDER`): S1 excludes gaps from the
denominator (primary, upper bound), S2 counts gaps inside the charted span as
downtime, S3 counts all node time (blocked on the discharge-timestamp defect). S1 is
primary because gaps are not random: mortality falls from 73% at no gap to 53% at
>6h of gap, so S2 loads the exposure with recovery status. **Consequence for 02: the
node covariates must carry the liberation predictors** (urine output, vasopressors,
fluid balance, SOFA), or S2 is confounded in a direction already known.

## Methods are inherited BY REFERENCE. Do not copy them.

The methodological source of truth is the sibling **fluid_ARDS** project. Those
are live documents under active edit; a copy here would drift and then two
projects would disagree about the same estimator.

| Path | What to read it for |
|---|---|
| `~/Desktop/Research/CLIF/fluid_ARDS/plans/lmtp_discussion.md` (965 lines) | S3.1 estimand vs estimator; S4 Assumption 5 and the self-limiting shift; S5 competing risks; S6 `lmtp` v1.5.4 mechanics and `lmtp_sdr` vs `lmtp_tmle` |
| `~/Desktop/Research/CLIF/fluid_ARDS/plans/analysis_plan.md` (1,518 lines) | S4.5 estimation; S4.9 the difference curve; S4.10 federated export contract and site-admission layers |
| `~/Desktop/Research/CLIF/fluid_ARDS/references/` | Diaz/Schenk MTP paper, `lmtp` package manual |

**Record only CRRT-specific deviations here.** Everything else is a pointer.

## Design decisions

Machine-readable in `config/lmtp_design.json`; that file is the
`definition_version` source for output provenance. Summary of the reasoning:

| Decision | Value | Why |
|---|---|---|
| delta primary | **5 mL/kg/hr** | Clean at 8 of 10 sites, 58-98% coverage. delta=10 is diluted at two low-dose sites and unstable at two tight-dosing sites |
| delta ladder | **{2.5, 5, 10}** | The ladder *is* the dose-response result. Report all three |
| FLOOR | **15**, sensitivity 20 | FLOOR=20 collapses coverage to 4-11% at two sites |
| Shift form | **self-limiting, never clamped** | `lmtp` supplies no guard and will not warn you. Clamping breaks invertibility |
| tau | **3 nodes at 0/24/48h** | Retention >= 38% at 48h everywhere; 72h thin at one site. Do not exceed tau=4 |
| Node statistic | **time-weighted mean** | Delivered dose. Charted zeros in, uncharted gaps not imputed. Matches Quickfall/Koyner 2026 at the same institution. Supersedes windowed median, whose charting-lag rationale does not survive a 24h node |
| Estimator | `lmtp_sdr` primary, `lmtp_tmle` secondary | Diaz's own applied choice; g-comp and IPW as diagnostics |
| Competing event | **discharge alive**, via `compete` | Mortality is in-hospital by construction |
| Site admission | **band occupancy under the shift** | Outcome-blind and pre-fit; replaces the `<100 high-dose-arm` rule |
| `EXCLUDE_SHORT_CRRT` | **dropped** | Only existed because a point-treatment design needs 24h to define exposure |

### The result worth remembering

At delta = 2.5 the binned density ratio never exceeds 6.0 at any of the ten
sites, and no site has more than 3.1% of its cohort in a thin bin. The archived
MSM had 29.5% of patients above propensity 0.95 at t = 12h and abandoned
truncation because it "did not help". That is a structural win, not a tuning win.

### The residual risk Tier 1 could not see

Marginal support is not conditional support. Assumption 1 is conditional on
history, and with A2 | A1 near-deterministic (r = 0.75) the conditional density
at the shifted value is thinner than the marginal audit suggests. **Settle this
in the single-site pilot** by fitting the treatment mechanism at node 2 and
inspecting the density-ratio distribution conditional on A1, not marginally.

Related: 91% of exposure variance at the coordinating site is between-patient
(within-patient SD 0.87 against between-patient 14.93). The time-varying design
is carried by the roughly 20-30% of patients whose dose actually moves, so
effective sample size at nodes 2 and 3 is far below n.

## The cohort is written here, not vendored

An earlier plan was to copy `00_cohort.py` and four dependencies verbatim from
`CLIF-epidemiology-of-CRRT`. **That was abandoned on 2026-08-16.** That file is a
converted Jupyter notebook (55 cells, 43% blanks and comments, ~900 lines of
diagnostics) with no callable API, so vendoring it meant adopting 3,281 lines to
reach roughly 600 of actual logic. Separately, `clifpy` 0.4.9 already provides
`stitch_encounters`, `compute_sofa_polars`, `create_wide_dataset` and
`apply_outlier_handling`, which covered most of what the manifest listed.

`code/01_build_cohort.py` is therefore written against clifpy primitives.
`code/vendor/` retains the pin and the byte-equality test for the two **config**
files taken from that repo at SHA `ee4774b`
(`clif_data_requirements.yaml`, `outlier_config.json`).

**The consequence to keep in mind:** because the cohort logic is reimplemented
rather than copied, matching the sibling's N is a **claim to be tested**, not a
fact. Current standing is 2,144 blocks against its 2,145, with one block of the
difference explained by a documented ESRD-grain divergence. See
`config/lmtp_design.json` under `cohort._must_reconcile`.

**Do not inherit `04_build_causal_df.py`'s schema.** It reads
`tableone_analysis_df.parquet`, so its baseline block is coupled to the Table 1
script, and it **violates the `L_t -> A_t` ordering** by pairing a 0-24h mean dose
with covariates measured at 24h.

## Protocol settings versus site settings

Every configurable value belongs to exactly one of two files, and putting one in
the wrong place is a scientific error rather than a tidiness problem.

| | File | Contains | Differs by site? |
|---|---|---|---|
| **Site** | `config/config.json` | `data_directory`, `filetype`, `timezone`, `site_name`, `has_crrt_settings` | **Yes** |
| **Protocol** | `config/lmtp_design.json` | delta ladder, floor, node schedule, study window, estimator, competing event, cohort rules | **No.** Identical everywhere |

The test: **if two sites set this differently, is the pooled result still
meaningful?** If no, it is protocol.

Never hardcode a protocol value in a script. Scripts **read** from
`lmtp_design.json`; they do not decide. A hardcoded constant is invisible to
everyone who did not write it and was agreed by nobody.

`lmtp_design.json` is the `definition_version` source stamped onto every
shareable output, so any change to it is a protocol amendment that bumps the
version and leaves earlier results identifiable.

Worked example: the study window (2018-2024) is protocol, not a site preference.
Calendar time confounds here, since CRRT practice and mortality moved through
2020-2021; per-site windows would make between-site heterogeneity partly reflect
*when* a site contributed. A site with shorter coverage still participates, and
the cohort script reports its actual range so partial coverage is stated rather
than silent.

## Repo conventions

Read `.gitignore`'s header before adding a rule to it. This repo is **private**,
This repo is **public**, because it ships to consortium sites. `.claude/` and
`docs/` are therefore **gitignored symlinks** into the private
`crrt-manuscript-tools` repo (`lmtp-claude/`, `lmtp-docs/`), which holds the
feasibility audit, the internal planning, and the cohort tutorial. All three carry
per-site or coordinating-site numbers that cannot be public. Both were purged from
this repo's history on 2026-08-18; do not restore them by un-ignoring.

`CLAUDE.md` and `README.md` ship, and are written with no site named.

- Python: **uv**, `requires-python >=3.11,<3.12`, `clifpy==0.4.9` exact. The
  compute stack is pinned exactly for cross-site numeric reproducibility;
  presentation libraries float.
- R: this repo builds **its own** `renv.lock` against system R 4.3.1. The
  sibling's lock targets 4.5.2 and `renv::restore()` fails against it here. See
  `code/R_PACKAGES.md`.
- `has_crrt_settings: true` is **required**; vendored `00_cohort.py` needs flow
  rates and mode columns to compute dose at all.
- Outputs split into `output/final_no_phi/` (aggregate, shareable) and
  `output/intermediate_phi/` (patient-level, never leaves the site).

## Data safety

**Never read raw data files** (parquet, csv) under any CLIF data directory. They
hold protected patient data that must not enter the conversation context. Write
code that reads the data and print only aggregated summaries. This applies to
Read, Glob, Grep, Bash, and to any subagent.

## Inherited facts that bear on this analysis

Memory is path-keyed, so these did not follow the code from the sibling project.
They are also written as memory files under this repo's memory directory.

- **Dose is the median of the first 3 hours, not the first charted value.** At
  one participating institution the first value runs about 12 mL/kg/hr low from
  charting lag; 62% of patients there ramp up. Other sites move by under 2.
- **Dose is modality-agnostic as of 2026-08-13** (`e9763a8`): it sums dialysate
  plus pre- and post-filter replacement for every dose-eligible mode, counting
  whichever are charted. SCUF excluded. **The coordinating site cannot validate
  this change** (CVVHD-only, never charts replacement fluid, so old and new
  formulas are bit-identical there). Confirming it by running that site confirms
  nothing.
- **Mortality is in-hospital by construction.** Death comes from discharge
  disposition, hospice counts as a death (2.7%), `death_dttm` is missing for 37%
  of counted deaths with last vital sign as fallback. **The user has committed to
  stating these explicitly in any Methods text and asked to be reminded.**
- **Weight denominator is unresolved.** `weight_kg` was the top treatment
  predictor at all five sites examined and negative at every one, most likely a
  denominator artifact (prescribed on ideal or adjusted weight, computed on
  actual). If so, part of the exposure contrast is adiposity rather than
  delivered therapy.
- **Pooling method depends on the estimand.** Inverse-variance DerSimonian-Laird
  random effects treats each site as one unit and is the only form carrying CIs
  and I-squared; crude / N-weighted treats each patient as one unit and is for
  baseline descriptors. Mixing them in one sentence is the easy error.
- **SOFA here is SOFA-1** (Vincent 1996). The trailing `-1` in "(SOFA)-1" is
  intentional disambiguation from SOFA-2, not a typo. Do not "fix" it.

## Audience-facing writing

- **Anonymize consortium sites** in anything audience-facing: manuscripts,
  abstracts, posters, slides, public dashboards. Refer to the *number* of sites,
  which is expected to change between drafts. Per-site figures use "Site 1,
  Site 2, ...". Author affiliations are the only place real names appear.
- **No em-dashes** in audience-facing content. Use commas, semicolons, or
  periods. En-dashes in numeric ranges are fine. Internal working files under
  `.claude/` are exempt.
- **Cite a source for every number** in any `.md` drafting file, as an inline
  `<!-- src: ... -->` comment or a parenthetical pointer, so the user can verify
  a value in seconds without re-deriving it.
