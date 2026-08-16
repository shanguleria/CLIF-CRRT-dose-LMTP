"""Build the CRRT cohort for the LMTP analysis, from raw CLIF 2.1.0 tables.

Produces, per encounter block:
  - the encounter_block <-> hospitalization_id mapping
  - CRRT initiation datetime (the index event, t = 0)
  - weight at initiation (the dose denominator)
  - the effluent dose time series over the exposure window
  - outcomes: death flag, event datetime, 30-day mortality anchored to initiation

Deliberately NOT produced here: baseline labs, SOFA, Table 1 covariates. Those
are computed per exposure node in 02_build_lmtp_df.py, because the L_t -> A_t
ordering requires them measured at each node rather than once at baseline.

Built stage by stage. See docs/clif_cohort_tutorial.md for the walkthrough.

DATA SAFETY: this script reads protected patient data. Print aggregates only,
never rows. Outputs split into output/intermediate_phi/ (patient-level, stays at
the site) and output/final_no_phi/ (aggregate, shareable).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Stage 0: What do we actually have?
# ---------------------------------------------------------------------------

def stage_0_inspect():
    """Load clif_crrt_therapy and report its shape, modes, and completeness.

    Answers, before any cohort logic is written:
      1. rows, distinct hospitalizations, date range
      2. crrt_mode_category distribution (counts and percentages)
      3. missingness of the flow-rate columns, ideally broken down BY MODE
      4. whether recorded_dttm carries a timezone

    Returns whatever you find useful to carry forward; at this stage printing is
    the point, so a return value is optional.
    """
    raise NotImplementedError("Stage 0: see docs/clif_cohort_tutorial.md")


# ---------------------------------------------------------------------------
# Stage 1: Adults and study years
# ---------------------------------------------------------------------------

def stage_1_base_population():
    """Restrict to adults (age >= 18) admitted within the study years."""
    raise NotImplementedError("Stage 1: not yet designed")


# ---------------------------------------------------------------------------
# Stage 2: CRRT encounters and the encounter block
# ---------------------------------------------------------------------------

def stage_2_encounter_blocks():
    """Stitch hospitalizations into encounter blocks; keep blocks with CRRT."""
    raise NotImplementedError("Stage 2: not yet designed")


# ---------------------------------------------------------------------------
# Stage 3: The ESRD exclusion
# ---------------------------------------------------------------------------

def stage_3_exclude_esrd():
    """Drop encounter blocks with pre-existing end-stage renal disease."""
    raise NotImplementedError("Stage 3: not yet designed")


# ---------------------------------------------------------------------------
# Stage 4: CRRT initiation, the index event
# ---------------------------------------------------------------------------

def stage_4_crrt_initiation():
    """Define t = 0 as the first crrt_therapy record per encounter block."""
    raise NotImplementedError("Stage 4: not yet designed")


# ---------------------------------------------------------------------------
# Stage 5: Weight at initiation
# ---------------------------------------------------------------------------

def stage_5_weight():
    """Find the weight closest to initiation. This is the dose denominator."""
    raise NotImplementedError("Stage 5: not yet designed")


# ---------------------------------------------------------------------------
# Stage 6: The effluent dose
# ---------------------------------------------------------------------------

def stage_6_dose():
    """Compute delivered effluent dose in mL/kg/hr, modality-agnostic.

    Sums dialysate + pre-filter + post-filter replacement for every dose-eligible
    mode, counting whichever are charted. SCUF excluded.
    """
    raise NotImplementedError("Stage 6: not yet designed")


# ---------------------------------------------------------------------------
# Stage 7: Outcomes
# ---------------------------------------------------------------------------

def stage_7_outcomes():
    """Death, event datetime, and 30-day mortality anchored to CRRT initiation.

    In-hospital by construction: death comes from discharge disposition, hospice
    counts as a death, and death_dttm is missing for a large share of deaths.
    """
    raise NotImplementedError("Stage 7: not yet designed")


# ---------------------------------------------------------------------------
# Stage 8: Write output and STROBE counts
# ---------------------------------------------------------------------------

def stage_8_write():
    """Write patient-level artifacts and the shareable STROBE count table."""
    raise NotImplementedError("Stage 8: not yet designed")


if __name__ == "__main__":
    stage_0_inspect()
