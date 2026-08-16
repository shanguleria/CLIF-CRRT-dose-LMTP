# code

| Path | Language | Status | Description |
|---|---|---|---|
| `vendor/00_cohort.py` | Python | not vendored | Cohort identification, ESRD exclusion, CRRT initiation per encounter block, outcomes |
| `vendor/01_create_wide_df.py` | Python | not vendored | Wide time-indexed labs, vitals, meds, respiratory support, NEE |
| `vendor/sofa_calculator.py` | Python | not vendored | Polars SOFA-1; accepts an arbitrary start/end per encounter, so a node-specific SOFA is a parameter, not new code |
| `vendor/pipeline_helpers.py` | Python | not vendored | Config loading and validation, intermediate IO |
| `vendor/utils.py` | Python | not vendored | CRRT outlier handling |
| `02_build_lmtp_df.py` | Python | not built | Exposure and covariate nodes at 0/24/48h |
| `03_lmtp_fit.R` | R | not built | `lmtp_sdr` over the delta ladder, influence-function exports |

## Vendored code

Files under `vendor/` are copied verbatim from `CLIF-epidemiology-of-CRRT` at the
commit pinned in `vendor/VENDOR_SHA`, and are **never edited in place**.
`tests/test_vendor_integrity.py` asserts byte-equality against that pin.

To take an upstream change: move the SHA, run `bash vendor/sync_vendor.sh`, run
the tests, and commit all of it together.

## Two constraints for step 02

1. **Enforce `L_t -> A_t` ordering.** Covariates at a node must precede the
   exposure at that node. The sibling repo's `04_build_causal_df.py` violates
   this by pairing a 0-24h mean dose with covariates measured at 24h. Do not
   reproduce that.
2. **Node statistic is the windowed median, not the mean.** The archived MSM used
   means, which reintroduce the dose charting-lag artefact that the median exists
   to defuse.

`.claude/feasibility_code/tier2_within_patient.py` is a working prototype of the
node assembly and is the starting point for step 02.

R environment: see `R_PACKAGES.md`.
