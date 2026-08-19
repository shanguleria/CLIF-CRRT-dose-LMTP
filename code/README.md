# code

| Path | Language | Status | Description |
|---|---|---|---|
| `vendor/00_cohort.py` | Python | not vendored | Cohort identification, ESRD exclusion, CRRT initiation per encounter block, outcomes |
| `vendor/01_create_wide_df.py` | Python | not vendored | Wide time-indexed labs, vitals, meds, respiratory support, NEE |
| `vendor/sofa_calculator.py` | Python | not vendored | Polars SOFA-1; accepts an arbitrary start/end per encounter, so a node-specific SOFA is a parameter, not new code |
| `vendor/pipeline_helpers.py` | Python | not vendored | Config loading and validation, intermediate IO |
| `vendor/utils.py` | Python | not vendored | CRRT outlier handling |
| `01_build_cohort.py` | Python | **runs** | Cohort, ESRD exclusion, CRRT initiation, weight, dose series, outcomes |
| `02_build_lmtp_df.py` | Python | **runs** | Exposure and covariate nodes at 0/24/48h, wide frame for `lmtp` |
| `03_lmtp_fit.R` | R | not built | `lmtp_sdr` over the delta ladder, influence-function exports |

## Vendored code

Files under `vendor/` are copied verbatim from `CLIF-epidemiology-of-CRRT` at the
commit pinned in `vendor/VENDOR_SHA`, and are **never edited in place**.
`tests/test_vendor_integrity.py` asserts byte-equality against that pin.

To take an upstream change: move the SHA, run `bash vendor/sync_vendor.sh`, run
the tests, and commit all of it together.

## Two constraints step 02 is built around

1. **`L_t -> A_t` ordering.** Covariates at a node must precede the exposure at
   that node. The sibling repo's `04_build_causal_df.py` violates this by pairing
   a 0-24h mean dose with covariates measured at 24h. 02 does not merely avoid
   this: it **asserts** the window bounds against the node schedule at import,
   before any data is read, and the assertion has been seen to fire.
2. **Node statistic is the time-weighted mean** over charted CRRT records,
   `sum(dose_h) / n_charted_hours`, computed on an hourly-bin reconstruction.
   This superseded the windowed median on 2026-08-18: the exposure is *delivered*
   dose, and the charting-lag artefact that motivated a median is bounded at
   -0.63 mL/kg/hr over a 24h node. Charted zeros enter; uncharted gaps are not
   imputed. The charted-zero versus uncharted-gap ambiguity is carried as a
   pre-specified bracket (S1 primary, S2 alongside), not a single choice.

Step 02's build walkthrough is `docs/lmtp_df_build_notes.md` (private repo).
`.claude/feasibility_code/tier2_within_patient.py` was the node-assembly prototype;
it targets the sibling repo's intermediate files and is superseded.

R environment: see `R_PACKAGES.md`.
