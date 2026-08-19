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

Built stage by stage. The walkthrough lives in the private repo at
crrt-manuscript-tools/lmtp-docs/clif_cohort_tutorial.md (symlinked here as docs/).

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
import subprocess
from datetime import datetime, timezone as _tz
from pathlib import Path

import pandas as pd
from clifpy import Adt, CrrtTherapy, Hospitalization, HospitalDiagnosis, Patient, Vitals, stitch_encounters

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

DOSE_FLOWS = ["dialysate_flow_rate",
              "pre_filter_replacement_fluid_rate",
              "post_filter_replacement_fluid_rate"]

_kw = dict(
    data_directory=config["data_directory"],
    filetype=config["filetype"],
    timezone=config["timezone"],
    output_directory=str(REPO_ROOT / "output"),
)

outliers = json.loads((REPO_ROOT / "config" / "outlier_config.json").read_text())
WEIGHT_LO, WEIGHT_HI = outliers["weight_kg"]

crrt = CrrtTherapy.from_file(**_kw).df
hosp = Hospitalization.from_file(**_kw).df
adt = Adt.from_file(**_kw).df
diag = HospitalDiagnosis.from_file(**_kw).df
vit = Vitals.from_file(
    **_kw,
    columns=["hospitalization_id", "recorded_dttm", "vital_category", "vital_value"],
    filters={"vital_category": ["weight_kg"]},
).df

pat = Patient.from_file(**_kw).df

# Defining the outcome of mortality
DEATH_DISPOSITIONS = {"expired", "hospice"}
OUTCOME_HORIZON_D = design["time"]["outcome_horizon_days"]


# Flow bounds are applied ONCE, here, before anything reads a flow. 
for _col in DOSE_FLOWS:
    _lo, _hi = outliers[_col]
    _bad = crrt[_col].notna() & ~crrt[_col].between(_lo, _hi)
    if _bad.any():
        print(f"  {_col}: {int(_bad.sum()):,} values outside [{_lo}, {_hi}] set to null")
    crrt.loc[_bad, _col] = pd.NA

# The arithmetic ceiling a dose can reach IF the bounds above were applied: every
# dose-eligible flow at its maximum, over the minimum plausible weight. Derived, not
# written down, so it cannot drift when a bound changes. This is not a clinical
# limit (a real dose is 20-35); it is the point past which a value is proof that
# the bounds did not run.
DOSE_CEILING = sum(outliers[c][1] for c in DOSE_FLOWS) / outliers["weight_kg"][0]
print(f"  bounded-dose arithmetic ceiling: {DOSE_CEILING:,.0f} mL/kg/hr")

print(f"crrt: {len(crrt):,} rows   hosp: {len(hosp):,} rows")
print(f"adt: {len(adt):,} rows")
print(f"diag: {len(diag):,} rows   esrd codes: {len(ESRD_CODES)}")
print(f"weight rows: {len(vit):,}")

# Defining Tau, from the design config
NODE_HOURS = design["time"]["exposure_node_hours"]
EXPOSURE_WINDOW_H = max(NODE_HOURS) + 24
print(f"nodes {NODE_HOURS}, window 0-{EXPOSURE_WINDOW_H}h")

# ---------------------------------------------------------------------------
# Utility: view the whole ladder so far
# ---------------------------------------------------------------------------
# %%
# Which stages chain from the one before, i.e. their first row should equal the
# previous stage's last row. Stage 1 is a standalone report and Stage 2 restarts
# from the full hospitalization table, so neither chains.
STROBE_STAGES = [
    (1, "Stage 1   hospitalizations", False),
    (2, "Stage 2   encounter blocks", False),
    (3, "Stage 3   ESRD exclusion",   True),
    (4, "Stage 4   CRRT initiation",  True),
    (5, "Stage 5   weight",           True),
    (6, "Stage 6   dose",             True),
    (7, "Stage 7   outcomes",         True),
]

def show_strobe():
    """Print every stage's STROBE rows as one ladder. Takes no arguments.

    Picks up whatever strobe_1 .. strobe_7 exist in the module namespace, so it
    can be run at any point and simply reports later stages as not built yet.

    Two things it checks, not just displays:
      - the change column, which is a real exclusion count because every row is a
        survivor count rather than a mixture of survivors and exclusions
      - continuity, since a chaining stage's first row must equal the previous
        stage's last row. A mismatch means a stage was handed the wrong frame
    """
    W, g = 44, globals()
    print("\n" + "=" * (W + 22))
    print("STROBE ladder, all stages")
    print("=" * (W + 22))
    last = None
    for num, label, chains in STROBE_STAGES:
        rows = g.get(f"strobe_{num}")
        if not rows:
            print(f"\n{label}   (not built yet)")
            continue
        print(f"\n{label}")
        if chains and last is not None and rows[0][1] != last:
            print(f"  !! expected to start at {last:,}, got {rows[0][1]:,}"
                  f"  <- this stage was handed the wrong frame")
        prev = None
        for name, n in rows:
            change = ""
            if prev is not None:
                d = n - prev
                change = "         same" if d == 0 else f"{d:>+14,}"
            print(f"  {name:<{W}} {n:>9,}{change}")
            prev = n
        last = rows[-1][1]
    print("=" * (W + 22))


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
    print(f"  stage 1 input: {len(hosp):,} hospitalizations")

    # Facts about the INPUTS, captured before any filtering.
    site_years = hosp["admission_dttm"].dt.year
    site_first, site_last = site_years.min(), site_years.max()
    all_hosp_ids = set(hosp["hospitalization_id"])

    strobe = []
    strobe.append(("all hospitalizations", len(hosp))) # Start with all hospitalizations

    # `eligible` is the narrowing frame; `hosp` stays untouched so the parameter
    # always means what the caller passed in.
    eligible = hosp[hosp["age_at_admission"] >= 18]
    strobe.append(("adults (age >= 18)", len(eligible))) # Filter for only adults over 18yo

    year = eligible["admission_dttm"].dt.year
    eligible = eligible[(year >= STUDY_YEAR_START) & (year <= STUDY_YEAR_END)]
    strobe.append((f"admitted {STUDY_YEAR_START}-{STUDY_YEAR_END}", len(eligible))) # Filter for study years

    crrt_ids = set(crrt["hospitalization_id"])
    eligible = eligible[eligible["hospitalization_id"].isin(crrt_ids)]
    strobe.append(("received CRRT", len(eligible))) # Filter for only those hospitalizations with a CRRT record

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

    print(f"  stage 2 input: {len(hosp):,} hospitalizations")

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

    return blocks, strobe, mapping

# %%
blocks, strobe_2, mapping = stage_2_encounter_blocks(hosp, adt, crrt)

# ---------------------------------------------------------------------------
# Stage 3: The ESRD exclusion
# ---------------------------------------------------------------------------

# %%
def stage_3_exclude_esrd(blocks, diag, mapping):
    """Drop encounter blocks with pre-existing end-stage renal disease.

    Applied at the encounter_block level, as ESRD is patient-specific at that
    time scale: an ESRD code on ANY hospitalization in a block excludes the whole
    block. Returns the surviving blocks and the strobe additions.
    """
    print(f"  stage 3 input: {len(blocks):,} blocks")

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

    # Output gets its own name, matching what the caller assigns it to. The
    # parameter `blocks` is never reassigned, so it always means what was passed in.
    n_before = len(blocks)
    blocks_no_esrd = blocks[~blocks["encounter_block"].isin(esrd_blocks)]
    n_excluded = n_before - len(blocks_no_esrd)
    assert n_excluded + len(blocks_no_esrd) == n_before, "blocks lost outside the filter"

    # Every strobe row is a count of what REMAINS; exclusions are the gaps
    # between rows. That keeps the ladder readable top to bottom and makes each
    # stage's first row the input check for the stage before it.
    strobe.append(("blocks entering ESRD exclusion", n_before))
    strobe.append(("blocks without ESRD", len(blocks_no_esrd)))

    print("\nSTROBE, stage 3")
    for label, n in strobe:
        print(f"  {label:<36} {n:>9,}")

    return blocks_no_esrd, strobe


# %%
blocks_no_esrd, strobe_3 = stage_3_exclude_esrd(blocks, diag, mapping)

# ---------------------------------------------------------------------------
# Stage 4: CRRT initiation, the index event
# ---------------------------------------------------------------------------
# %%
print(f"  stage 4 input: {len(blocks_no_esrd):,} blocks")
def stage_4_crrt_initiation(blocks_no_esrd, crrt, mapping):
    """Define t = 0 as the first crrt_therapy record per encounter block in which 
    patient actually receives a nonzero CRRT dose."""
    strobe = []

    c_all = crrt.merge(mapping, on="hospitalization_id", how="left", validate="many_to_one")
    c_all = c_all[c_all["encounter_block"].isin(set(blocks_no_esrd["encounter_block"]))]

    # Initiation of CRRT as the first record delivering a CRRT dose
    delivering = False
    for col in DOSE_FLOWS:
        delivering = delivering | (c_all[col].notna() & (c_all[col] > 0))

    init = (c_all[delivering].groupby("encounter_block")["recorded_dttm"].min()
            .rename("crrt_initiation_dttm").reset_index())


    n_before = len(blocks_no_esrd)
    blocks_with_init = blocks_no_esrd.merge(init, on="encounter_block", how="inner",
                                            validate="one_to_one")
    strobe.append(("blocks entering stage 4", n_before))
    strobe.append(("blocks with a CRRT initiation time", len(blocks_with_init)))
    assert blocks_with_init["crrt_initiation_dttm"].notna().all(), "block without t=0"

    # Determining the hours from admission to CRRT start
    # Includes warnings for CRRT start times before admission or after discharge, which suggest ETL issues
    late = blocks_with_init["crrt_initiation_dttm"] > blocks_with_init["block_discharge_dttm"]
    early = blocks_with_init["crrt_initiation_dttm"] < blocks_with_init["block_admission_dttm"]
    if late.any() or early.any():
        over = ((blocks_with_init.loc[late, "crrt_initiation_dttm"]
                    - blocks_with_init.loc[late, "block_discharge_dttm"]).dt.total_seconds()/3600)
        print(f"  WARNING: {int(early.sum())} blocks start CRRT before admission, "
                f"{int(late.sum())} after discharge "
                f"(median {over.median():.1f}h past discharge)")

    hrs = ((blocks_with_init["crrt_initiation_dttm"] - blocks_with_init["block_admission_dttm"])
            .dt.total_seconds() / 3600)
    print(f"  hours from admission to CRRT start: median {hrs.median():.1f}, "
            f"p5 {hrs.quantile(.05):.1f}, p95 {hrs.quantile(.95):.1f}")

    # What was actually charted for the blocks that got no t=0? 
    # Grouped by the actual dose, not crrt_mode_category
    no_init = set(blocks_no_esrd["encounter_block"]) - set(init["encounter_block"])
    d = c_all[c_all["encounter_block"].isin(no_init)]
    per_block = d.groupby("encounter_block").agg(
        n_rows=("recorded_dttm", "size"),
        any_uf=("ultrafiltration_out", lambda s: (s.notna() & (s > 0)).any()),
        any_bf=("blood_flow_rate",     lambda s: (s.notna() & (s > 0)).any()),
    )
    ran = per_block["any_uf"] | per_block["any_bf"]
    uf_only     = int((ran & (per_block["n_rows"] >= 6)).sum())
    brief_trace = int((ran & (per_block["n_rows"] <  6)).sum())
    no_evidence = int((~ran).sum())
    print(f"  no t=0 for {len(no_init)} blocks: {uf_only} ultrafiltration-only courses, "
            f"{brief_trace} brief traces (1-5 rows), {no_evidence} with no machine activity")

    print("\nSTROBE, stage 4")
    for label, n in strobe:
        print(f"  {label:<38} {n:>9,}")

    return blocks_with_init, strobe

# %%
blocks_with_init, strobe_4 = stage_4_crrt_initiation(blocks_no_esrd, crrt, mapping)

# ---------------------------------------------------------------------------
# Stage 5: Weight at initiation
# ---------------------------------------------------------------------------
# %%
def stage_5_weight(blocks_with_init, vit, mapping): 
    """Find the weight closest to initiation. This is the dose denominator."""
    print(f"  stage 5 input: {len(blocks_with_init):,} blocks")
    strobe = []

    w = vit.rename(columns={"vital_value": "weight_kg"})
    n_raw = len(w)
    w = w[w["weight_kg"].between(WEIGHT_LO, WEIGHT_HI)]
    print(f"  weight rows: {n_raw:,} -> {len(w):,} inside [{WEIGHT_LO}, {WEIGHT_HI}] kg")

    w = w.merge(mapping, on="hospitalization_id", how="inner", validate="many_to_one")
    w = w[w["encounter_block"].isin(set(blocks_with_init["encounter_block"]))]

    # Sort by dttm and align, then find the most recent weight recorded before CRRT initiation
    left = (blocks_with_init[["encounter_block", "crrt_initiation_dttm"]]
            .sort_values("crrt_initiation_dttm"))
    right = w[["encounter_block", "recorded_dttm", "weight_kg"]].sort_values("recorded_dttm")

    back = pd.merge_asof(left, right,
                            left_on="crrt_initiation_dttm", right_on="recorded_dttm",
                            by="encounter_block", direction="backward")
    # Some weights are first recorded after CRRT start, so use that as a fall back
    missing = back.loc[back["weight_kg"].isna(),
                        ["encounter_block", "crrt_initiation_dttm"]]
    fwd = pd.merge_asof(missing.sort_values("crrt_initiation_dttm"), right,
                        left_on="crrt_initiation_dttm", right_on="recorded_dttm",
                        by="encounter_block", direction="forward")

    got = pd.concat([back[back["weight_kg"].notna()], fwd[fwd["weight_kg"].notna()]])
    blocks_with_weight = blocks_with_init.merge(
        got[["encounter_block", "weight_kg", "recorded_dttm"]]
            .rename(columns={"recorded_dttm": "weight_dttm"}),
        on="encounter_block", how="inner", validate="one_to_one")

    strobe.append(("blocks entering stage 5", len(blocks_with_init)))
    strobe.append(("blocks with a weight", len(blocks_with_weight)))

    # How stale is that weight value? 
    lag = ((blocks_with_weight["crrt_initiation_dttm"]
            - blocks_with_weight["weight_dttm"]).dt.total_seconds() / 3600)
    print(f"  hours from weight to initiation: median {lag.median():.1f}, "
            f"p95 {lag.quantile(.95):.1f}, max {lag.max():.1f}")
    print(f"  weights taken AFTER initiation (the fallback): {int((lag < 0).sum())}")
    print(f"  weight kg: median {blocks_with_weight['weight_kg'].median():.1f}, "
            f"p1 {blocks_with_weight['weight_kg'].quantile(.01):.1f}, "
            f"p99 {blocks_with_weight['weight_kg'].quantile(.99):.1f}")

    print("\nSTROBE, stage 5")
    for label, n in strobe:
        print(f"  {label:<38} {n:>9,}")

    return blocks_with_weight, strobe

# %%
blocks_with_weight, strobe_5 = stage_5_weight(blocks_with_init, vit, mapping)

# ---------------------------------------------------------------------------
# Stage 6: The effluent dose
# ---------------------------------------------------------------------------
# %%
def stage_6_dose(blocks_with_weight, crrt, mapping):
    """Compute delivered effluent dose in mL/kg/hr, modality-agnostic.

    Sums dialysate + pre-filter + post-filter replacement for every dose-eligible
    mode in a modality-agnostic fashion.
    """
    print(f"  stage 6 input: {len(blocks_with_weight):,} blocks")
    strobe = []

    d = crrt.merge(mapping, on="hospitalization_id", how="left", validate="many_to_one")
    d = d.merge(blocks_with_weight[["encounter_block", "crrt_initiation_dttm",
                                        "weight_kg"]],
                    on="encounter_block", how="inner", validate="many_to_one")


    d["effluent_ml_hr"] = d[DOSE_FLOWS].sum(axis=1, min_count=1)
    n_no_flow = int(d["effluent_ml_hr"].isna().sum())
    d = d[d["effluent_ml_hr"].notna()]
    # Go from mL/hr to mL/kg/hr based on most recent weight
    d["dose_ml_kg_hr"] = d["effluent_ml_hr"] / d["weight_kg"]
    d["hours_from_init"] = ((d["recorded_dttm"] - d["crrt_initiation_dttm"])
                            .dt.total_seconds() / 3600)

    n_pre = len(d)
    dose_series = d.loc[d["hours_from_init"].between(0, EXPOSURE_WINDOW_H),
                        ["encounter_block", "recorded_dttm", "hours_from_init",
                            "effluent_ml_hr", "weight_kg", "dose_ml_kg_hr"]].copy()
    print(f"  records inside the 0-{EXPOSURE_WINDOW_H}h window: "
            f"{len(dose_series):,} of {n_pre:,}")

    with_dose = set(dose_series["encounter_block"])
    blocks_with_dose = blocks_with_weight[
        blocks_with_weight["encounter_block"].isin(with_dose)]
    assert len(blocks_with_dose) == len(blocks_with_weight), (
            f"{len(blocks_with_weight) - len(blocks_with_dose)} blocks lost a dose "
            "series they should structurally have; check the flow bounds")

    strobe.append(("blocks entering stage 6", len(blocks_with_weight)))
    strobe.append(("blocks with a dose series", len(blocks_with_dose)))

    # Dose diagnostics at each Tau node
    z = int((dose_series["dose_ml_kg_hr"] == 0).sum())
    print(f"  dose mL/kg/hr: median {dose_series['dose_ml_kg_hr'].median():.1f}, "
            f"p1 {dose_series['dose_ml_kg_hr'].quantile(.01):.1f}, "
            f"p99 {dose_series['dose_ml_kg_hr'].quantile(.99):.1f}, "
            f"max {dose_series['dose_ml_kg_hr'].max():.1f}")
    print(f"  exactly zero: {z:,} records ({100 * z / len(dose_series):.1f}%)")

    for h in NODE_HOURS:
        n = dose_series.loc[dose_series["hours_from_init"].between(h, h + 24),
                            "encounter_block"].nunique()
        print(f"    node {h:>2}-{h + 24:>2}h: {n:,} blocks "
                f"({100 * n / len(blocks_with_dose):.0f}%)")

    # Not a clinical limit (a real dose is 20-35). This is the point past which a
    # value proves the flow bounds did not run: every dose-eligible flow at its
    # ceiling over the minimum plausible weight.
    assert dose_series["dose_ml_kg_hr"].max() <= DOSE_CEILING, (
        f"max dose {dose_series['dose_ml_kg_hr'].max():,.1f} exceeds the arithmetic "
        f"ceiling of {DOSE_CEILING:,.0f} mL/kg/hr; the flow outlier bounds did not run "
        "(re-run the config cell if you are in an interactive session)")

    print("\nSTROBE, stage 6")
    for label, n in strobe:
        print(f"  {label:<38} {n:>9,}")

    return blocks_with_dose, dose_series, strobe

# %%
blocks_with_dose, dose_series, strobe_6 = stage_6_dose(blocks_with_weight, crrt, mapping)

# ---------------------------------------------------------------------------
# Stage 7: Outcomes
# ---------------------------------------------------------------------------
# %%
def stage_7_outcomes(blocks_with_dose, pat):
    """Death, event datetime, and 30-day mortality anchored to CRRT initiation.

    In-hospital by construction: death comes from discharge disposition due to unreliable death_dttm
    , hospice discharge counts as a death.
    """
    print(f"  stage 7 input: {len(blocks_with_dose):,} blocks")
    strobe = []

    b = blocks_with_dose.merge(pat[["patient_id", "death_dttm"]],
                                on="patient_id", how="left", validate="many_to_one")

    disp = b["discharge_category"].str.lower()
    b["died_in_hospital"] = disp.isin(DEATH_DISPOSITIONS)

    # Death as defined by discharge not death_dttm
    b["event_dttm"] = b["block_discharge_dttm"].where(b["died_in_hospital"])
    b["days_to_event"] = ((b["event_dttm"] - b["crrt_initiation_dttm"])
                            .dt.total_seconds() / 86400)

    # Mortality in the hospital at 30d post initiation of CRRT
    b["mortality_30d"] = (b["died_in_hospital"]
                            & (b["days_to_event"] <= OUTCOME_HORIZON_D)).fillna(False)
    impossible = int((b["died_in_hospital"] & (b["days_to_event"] < 0)).sum())
    if impossible:
        print(f"  WARNING: {impossible} deaths timed BEFORE CRRT initiation. Revisit ETL.")

        n_died = int(b["died_in_hospital"].sum())
    disagree = b["died_in_hospital"] & b["death_dttm"].notna()
    gap_h = ((b["death_dttm"] - b["block_discharge_dttm"]).dt.total_seconds() / 3600)[disagree]

    print(f"  died in hospital: {n_died:,} ({100*n_died/len(b):.1f}%), "
            f"of which hospice {int((disp == 'hospice').sum()):,}")
    print(f"  death_dttm present for {int(disagree.sum()):,} of them "
            f"({100*disagree.sum()/n_died:.0f}%); vs discharge it differs by median "
            f"{gap_h.median():.1f}h, |diff|>24h for {int((gap_h.abs() > 24).sum())}")
    print(f"  days from initiation to death: median {b['days_to_event'].median():.1f}, "
            f"p95 {b['days_to_event'].quantile(.95):.1f}, max {b['days_to_event'].max():.1f}")

    strobe.append(("blocks entering stage 7", len(blocks_with_dose)))
    strobe.append(("blocks with an outcome", len(b)))
    print("\nSTROBE, stage 7")
    for label, n in strobe:
        print(f"  {label:<38} {n:>9,}")
    print(f"  {OUTCOME_HORIZON_D}-day in-hospital mortality: "
            f"{int(b['mortality_30d'].sum()):,} ({100*b['mortality_30d'].mean():.1f}%)")
    print(f"  in-hospital mortality, any time  : {n_died:,} ({100*n_died/len(b):.1f}%)")

    return b, strobe

# %%
cohort, strobe_7 = stage_7_outcomes(blocks_with_dose, pat)


# ---------------------------------------------------------------------------
# Stage 8: Write output and STROBE counts
# ---------------------------------------------------------------------------

def _code_version(repo_root):
    """git describe, or "unknown" where this is not a git checkout.

    Sites may receive the code as a zip rather than a clone. "unknown" is a worse
    answer than a SHA and a much better one than crashing at the final step of a
    long run.
    """
    try:
        return subprocess.run(
            ["git", "-C", str(repo_root), "describe", "--always", "--dirty"],
            capture_output=True, text=True, timeout=10, check=True,
        ).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def stage_8_write(cohort, dose_series, mapping, strobes, config, design, repo_root):
    """Write patient-level artifacts and the shareable STROBE count table.

    Computes nothing. The work is splitting outputs by whether they can leave the
    site, and stamping the shareable one with enough provenance to trace it back to
    the code and the estimand that produced it.
    """
    print(f"  stage 8 input: {len(cohort):,} blocks, {len(dose_series):,} dose records")

    # Five fields, answering the questions someone pooling ten sites will ask.
    # --dirty matters: it says the tree had uncommitted changes, so the SHA does
    # not fully identify the code that ran.
    prov = {
        "site_id": config["site_name"],
        "code_version": _code_version(repo_root),
        "clif_version": config["clif_version"],
        "definition_version": design["definition_version"],
        "generated": datetime.now(_tz.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    for k, v in prov.items():
        print(f"    {k:<20} {v}")

    phi = repo_root / "output" / "intermediate_phi"
    share = repo_root / "output" / "final_no_phi"
    phi.mkdir(parents=True, exist_ok=True)
    share.mkdir(parents=True, exist_ok=True)

    # Sorted before writing so a re-run is byte-comparable. Row order out of a
    # groupby is not guaranteed stable across pandas versions, and an artifact
    # whose bytes change without its contents changing cannot be hash-checked.
    cohort_out = cohort.sort_values("encounter_block").reset_index(drop=True)
    dose_out = (dose_series.sort_values(["encounter_block", "recorded_dttm"])
                           .reset_index(drop=True))
    cohort_out.to_parquet(phi / "cohort.parquet", index=False)
    dose_out.to_parquet(phi / "dose_series.parquet", index=False)
    print(f"    wrote intermediate_phi/cohort.parquet      {len(cohort_out):,} rows")
    print(f"    wrote intermediate_phi/dose_series.parquet {len(dose_out):,} rows")

    # 02 needs hospitalization_id -> encounter_block to pull covariates, since every
    # raw CLIF table is keyed on hospitalization_id and nothing downstream of here
    # carries one. Persisted rather than re-derived: re-running stitch_encounters in
    # 02 would be a second, independent claim about which hospitalizations belong to
    # which block, and two scripts that disagree about that disagree about the cohort.
    # Restricted to the final blocks, so it inherits the same PHI scope as cohort.parquet.
    map_out = (mapping[mapping["encounter_block"].isin(set(cohort_out["encounter_block"]))]
               .sort_values(["encounter_block", "hospitalization_id"])
               .reset_index(drop=True))
    assert set(map_out["encounter_block"]) == set(cohort_out["encounter_block"]), (
        "block map does not cover the cohort exactly")
    map_out.to_parquet(phi / "block_map.parquet", index=False)
    print(f"    wrote intermediate_phi/block_map.parquet   {len(map_out):,} rows "
          f"({map_out['encounter_block'].nunique():,} blocks)")

    # Long format, provenance repeated on every row: pooled files get concatenated
    # and filtered, and a header-only provenance block is lost the moment they are.
    rows = [{"stage": s, "step": i, "label": lab, "n": int(n)}
            for s, entries in strobes for i, (lab, n) in enumerate(entries, 1)]
    strobe_df = pd.DataFrame(rows)
    for k, v in prov.items():
        strobe_df[k] = v

    # Stage 8 is where an accident would leave the site. This cannot fire today;
    # it exists so a future edit adding a helpful "example encounter_block" column
    # stops the run rather than shipping silently. A timestamp would be equally
    # disqualifying: an admission time plus a site is re-identifying, a count is not.
    ID_LIKE = ("encounter_block", "hospitalization_id", "patient_id", "_dttm")
    leaked = [c for c in strobe_df.columns if any(s in c for s in ID_LIKE)]
    assert not leaked, f"identifier-like columns in a shareable file: {leaked}"
    assert strobe_df["n"].notna().all(), "a strobe row lost its count"

    path = share / f"{config['site_name']}_strobe_counts.csv"
    strobe_df.to_csv(path, index=False)
    print(f"    wrote final_no_phi/{path.name}  {len(strobe_df)} rows, PHI-checked")
    return prov


strobes = [(1, strobe_1), (2, strobe_2), (3, strobe_3), (4, strobe_4),
           (5, strobe_5), (6, strobe_6), (7, strobe_7)]
cohort_prov = stage_8_write(cohort, dose_series, mapping, strobes, config, design, REPO_ROOT)


# %%

# %%
show_strobe()
