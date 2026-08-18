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
from clifpy import Adt, CrrtTherapy, Hospitalization, HospitalDiagnosis, stitch_encounters

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

STITCH_HOURS = design["cohort"]["stitch_time_interval_hours"]

ESRD_CODES = {c["code"] for c in design["cohort"]["esrd_exclusion"]["codes"]}

_kw = dict(
    data_directory=config["data_directory"],
    filetype=config["filetype"],
    timezone=config["timezone"],
    output_directory=str(REPO_ROOT / "output"),
)

crrt = CrrtTherapy.from_file(**_kw).df
hosp = Hospitalization.from_file(**_kw).df
adt = Adt.from_file(**_kw).df
diag = HospitalDiagnosis.from_file(**_kw).df

print(f"crrt: {len(crrt):,} rows   hosp: {len(hosp):,} rows")
print(f"adt: {len(adt):,} rows")
print(f"diag: {len(diag):,} rows   esrd codes: {len(ESRD_CODES)}")

# ---------------------------------------------------------------------------
# Stage 0: What do we actually have?
# ---------------------------------------------------------------------------
# %%
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
def stage_2_encounter_blocks(hosp, adt, crrt):
    """Stitch hospitalizations into encounter blocks; keep blocks with CRRT.
    
    Takes the full hospitalization table (i.e. not Stage 1's output) and finds partners\n
    for CRRT encounters across the whole table. Returns one row per encounter block. """

    strobe = []
    # Creating DataFrames for a stitched hospitalization table, a stitched ADT table, and the mapping
    _, _, mapping = stitch_encounters(hosp, adt, time_interval=STITCH_HOURS)
    strobe.append(("hospitalizations before stitching", len(mapping)))
    strobe.append(("encounter blocks after stitching", mapping["encounter_block"].nunique()))

    # Attaching encounter block ID to hospitalization
    hb = hosp.merge(mapping, on="hospitalization_id", how="left", validate="one_to_one")
    assert hb["encounter_block"].notna().all(), "some hospitalizations assigned no block"

    # Find blocks that contain CRRT episodes
    crrt_ids = set(crrt["hospitalization_id"])
    crrt_blocks = set(hb.loc[hb["hospitalization_id"].isin(crrt_ids), "encounter_block"])
    hb = hb[hb["encounter_block"].isin(crrt_blocks)]
    strobe.append(("blocks with any CRRT record", len(crrt_blocks)))

    # Collapse from one row per hospitalization to one row per encounter block
    n_in = len(hb)
    hb = hb.sort_values(["encounter_block", "admission_dttm"])
    blocks = hb.groupby("encounter_block").agg(
        patient_id=("patient_id", "first"),
        n_hospitalizations=("hospitalization_id", "size"),
        block_admission_dttm=("admission_dttm", "min"),
        block_discharge_dttm=("discharge_dttm", "max"),
        age_at_admission=("age_at_admission", "first"),
    ).reset_index()
    assert len(blocks) == blocks["encounter_block"].nunique(), "duplicate blocks"
    assert blocks["n_hospitalizations"].sum() == n_in, "rows lost or duplicated"

    # Discharge dttm to end encounter block based on last hospitalization in the block
    last = (hb.sort_values(["encounter_block", "discharge_dttm"])
            .groupby("encounter_block")
            .last()[["discharge_category"]]
            .reset_index()
            )
    blocks = blocks.merge(last, on="encounter_block", how="left", validate="one_to_one")

    assert blocks["discharge_category"].notna().all(), "block lost its disposition"

    # STROBE Criteria by Encounter Block
    blocks = blocks[blocks["age_at_admission"] >= 18]
    strobe.append(("adult blocks", len(blocks)))

    y = blocks["block_admission_dttm"].dt.year
    blocks = blocks[(y >= STUDY_YEAR_START) & (y <= STUDY_YEAR_END)]
    strobe.append((f"admitted {STUDY_YEAR_START}-{STUDY_YEAR_END}", len(blocks)))

    print("\nSTROBE, stage 2")
    for label, n in strobe:
        print(f"  {label:<36} {n:>9,}")
    multi = (blocks["n_hospitalizations"] > 1).sum()
    print(f"\n  blocks built from >1 hospitalization: {multi:,}")

    return blocks, strobe

# %%
blocks, strobe_2 = stage_2_encounter_blocks(hosp, adt, crrt)

# ---------------------------------------------------------------------------
# Stage 3: The ESRD exclusion
# ---------------------------------------------------------------------------

# %%
def stage_3_exclude_esrd():
    """Drop encounter blocks with pre-existing end-stage renal disease.
    Applied at the encounter_block level as ESRD is patient-specific at that time scale."""
# %%
# Indent into the def function after testing
    strobe = []

    code = diag["diagnosis_code"].str.replace(".","",regex=False).str.lower()
    is_esrd = code.isin(ESRD_CODES)
    print(f"  diagnosis rows matching an ESRD code: {int(is_esrd.sum()):,}")
    # Present on admission, if uninformative then favor dropping ESRD to be conservative
    poa = diag["poa_present"]
    poa_informative = bool((poa == 1).any())
    if poa_informative:
        is_esrd = is_esrd & (poa == 1)
    else:
        print("  NOTE: poa_present has no positive values. Treating it as "
            "UNINFORMATIVE and matching on diagnosis code alone. With the POA "
            "condition applied, ZERO blocks would be excluded.")
    print(f"  after the POA rule:                   {int(is_esrd.sum()):,}")

    # Exclude ESRD at the block grain
    esrd_hosp_ids = set(diag.loc[is_esrd, "hospitalization_id"])
    esrd_blocks = set(mapping.loc[
        mapping["hospitalization_id"].isin(esrd_hosp_ids), "encounter_block"])

    n_before = len(blocks)
    blocks = blocks[~blocks["encounter_block"].isin(esrd_blocks)]
    n_excluded = n_before - len(blocks)
    assert n_excluded + len(blocks) == n_before, "blocks lost outside the filter"

    strobe.append(("blocks with pre-existing ESRD", n_excluded))
    strobe.append(("blocks without ESRD", len(blocks)))

    print("\nSTROBE, stage 3")
    for label, n in strobe:
        print(f"  {label:<36} {n:>9,}")

    return blocks, strobe

# ---------------------------------------------------------------------------
# Stage 4: CRRT initiation, the index event
# ---------------------------------------------------------------------------
# %%
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


# %%
