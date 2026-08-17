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
# %%
# ---------------------------------------------------------------------------
# IMPORT BLOCK
# ---------------------------------------------------------------------------
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd
from clifpy import CrrtTherapy, Hospitalization

try:
    REPO_ROOT = Path(__file__).resolve().parent.parent
except NameError:               # no __file__ in an interactive session
    REPO_ROOT = Path.cwd()

# %%
# ---------------------------------------------------------------------------
# CONFIG BLOCK
# ---------------------------------------------------------------------------

config = json.loads((REPO_ROOT / "config" / "config.json").read_text())
design = json.loads((REPO_ROOT / "config" / "lmtp_design.json").read_text())

STUDY_YEAR_START = design["cohort"]["study_year_start"]
STUDY_YEAR_END = design["cohort"]["study_year_end"]
print(f"study window: {STUDY_YEAR_START}-{STUDY_YEAR_END}")

_kw = dict(
    data_directory=config["data_directory"],
    filetype=config["filetype"],
    timezone=config["timezone"],
    output_directory=str(REPO_ROOT / "output"),
)

crrt = CrrtTherapy.from_file(**_kw).df
hosp = Hospitalization.from_file(**_kw).df
print(f"crrt: {len(crrt):,} rows   hosp: {len(hosp):,} rows")

# ---------------------------------------------------------------------------
# Stage 0: What do we actually have?
# ---------------------------------------------------------------------------

def stage_0_inspect(df):
    """Report the shape, modality mix, and completeness of clif_crrt_therapy.
    Prints aggregates only; nothing here may show patient-level rows.
    """
    # Exploring rows, hospitalizations, date ranges, time zone, modalities
    print(f"rows: {len(df):,}")
    print(f"distinct hospitalizations: {df['hospitalization_id'].nunique():,}")
    print(f"date range: {df['recorded_dttm'].min()} to {df['recorded_dttm'].max()}")
    print(f"recorded_dttm timezone: {df['recorded_dttm'].dtype}")

    counts = df['crrt_mode_category'].value_counts(dropna=False)
    pct = df['crrt_mode_category'].value_counts(dropna=False, normalize=True)*100
    modes = pd.DataFrame({"n": counts, "pct": pct.round(1)})
    print(f"CRRT modalities:\n{modes}")

    # Missingness
    flow_cols = [
        "blood_flow_rate",
        "dialysate_flow_rate",
        "pre_filter_replacement_fluid_rate",
        "post_filter_replacement_fluid_rate",
        "ultrafiltration_out",
    ]
    # dropna=False so a null modality becomes a visible row rather than rows
    # that silently vanish from the table.
    by_mode = (
        df.groupby("crrt_mode_category", dropna=False)[flow_cols]
        .apply(lambda g: g.isna().mean())
    )
    print(f"Missingness by Modality:\n{by_mode.T.round(3)}")
    print("Note that DFR will be missing for convective-only modalities \n"
    "and fluid replacement will be missing from diffusive-only modalities")

stage_0_inspect(crrt)

# ---------------------------------------------------------------------------
# Stage 1: Adults and study years
# ---------------------------------------------------------------------------
# %%
def stage_1_base_population(hosp, crrt):
    """Count how many hospitalizations are eligible: adults, in the study years,
    with a CRRT record.

    A reporting step, not a data-producing one. Returns the STROBE count ladder
    as a list of (label, count) pairs; Stage 2 builds the actual cohort from the
    full tables, because stitching must see the whole population.
    """
    # Facts about the INPUTS, captured before any filter reassigns hosp.
    site_years = hosp["admission_dttm"].dt.year
    site_first, site_last = site_years.min(), site_years.max()
    all_hosp_ids = set(hosp["hospitalization_id"])

    strobe = []
    strobe.append(("all hospitalizations", len(hosp))) # Start with all hospitalizations

    hosp = hosp[hosp["age_at_admission"] >= 18]
    strobe.append(("adults (age >= 18)", len(hosp))) # Filter for only adults over 18yo

    year = hosp["admission_dttm"].dt.year
    hosp = hosp[(year >= STUDY_YEAR_START) & (year <= STUDY_YEAR_END)]
    strobe.append((f"admitted {STUDY_YEAR_START}-{STUDY_YEAR_END}", len(hosp))) # Filter for study years

    crrt_ids = set(crrt["hospitalization_id"])
    hosp = hosp[hosp["hospitalization_id"].isin(crrt_ids)]
    strobe.append(("received CRRT", len(hosp))) # Filter for only those hospitalizations with a CRRT record

    # Referential integrity: CRRT records whose hospitalization does not exist.
    # The .isin above drops these silently, so count them deliberately.
    orphans = crrt_ids - all_hosp_ids

    # Print the STROBE ladder
    print(f"\nsite data covers {site_first}-{site_last} "
          f"(protocol window {STUDY_YEAR_START}-{STUDY_YEAR_END})")
    if orphans:
        print(f"WARNING: {len(orphans)} CRRT hospitalization_id(s) have no row in "
              f"clif_hospitalization and are excluded")

    print("\nSTROBE ladder")
    for label, n in strobe:
        print(f"  {label:<34} {n:>9,}")

    # Only the ladder. This is a reporting step: its deliverable is the
    # hospitalization-level rows of the STROBE diagram, not a table. Stage 2
    # works from the full hosp/crrt tables so that stitching can see a CRRT
    # encounter's partners.
    return strobe

strobe_1 = stage_1_base_population(hosp, crrt)

# ---------------------------------------------------------------------------
# Stage 2: CRRT encounters and the encounter block
# ---------------------------------------------------------------------------
# %%
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

