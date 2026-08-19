"""Build the LMTP analysis frame from 01's cohort, the dose series, and raw CLIF tables.

Produces the single wide data.frame that `lmtp` requires: one row per encounter
block, exposure columns ordered by time, time-varying covariates per node, and
outcome / competing-event / censoring indicators on the day 3/7/14/30 grid.

  A_t   delivered CRRT dose in node t, time-weighted mean, under two
        charted-zero-versus-uncharted-gap conventions (S1 primary, S2 sensitivity)
  L_0   baseline covariates over [-24h, 0h)
  L_t   time-varying covariates over the node PRECEDING the exposure node they feed
  Y     in-hospital death, D  discharge alive (competing), C  censoring

Every measurement rule applied here is READ from config/lmtp_design.json under
`covariates`. This script decides nothing. If a lookback or a summary rule needs to
change, change it there: it is protocol, it is hashed into definition_version, and a
constant hardcoded here would be invisible to everyone who did not write it.

THE ORDERING THAT MATTERS: L_t must be measured before A_t. Every covariate window
ends at or before the start of the exposure node it feeds, and the config block is
ASSERTED against the node schedule at import rather than trusted to a comment. The
sibling repo's 04_build_causal_df.py violates exactly this rule by pairing a 0-24h
mean dose with covariates measured at 24h.

The walkthrough lives in the private repo at
crrt-manuscript-tools/lmtp-docs/lmtp_df_build_notes.md (symlinked here as docs/).

DATA SAFETY: this script reads protected patient data. Print aggregates only, never
rows. Outputs split into output/intermediate_phi/ (patient-level, stays at the site)
and output/final_no_phi/ (aggregate, shareable).
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

import numpy as np
import pandas as pd
import polars as pl
import yaml

import clifpy
from clifpy import (HospitalDiagnosis, Labs, MedicationAdminContinuous, Patient,
                    RespiratorySupport, Vitals, calculate_cci, compute_sofa_polars)
from clifpy.utils.unit_converter import convert_dose_units_by_med_category

try:
    REPO_ROOT = Path(__file__).resolve().parent.parent
except NameError:               # no __file__ in an interactive session
    REPO_ROOT = Path.cwd()

# %%
# ---------------------------------------------------------------------------
# CONFIG BLOCK  —  everything below is READ, nothing is decided here
# ---------------------------------------------------------------------------

config = json.loads((REPO_ROOT / "config" / "config.json").read_text())
design = json.loads((REPO_ROOT / "config" / "lmtp_design.json").read_text())

COV = design["covariates"]
DEFS = COV["definitions"]

TAU = design["time"]["tau"]
NODE_HOURS = design["time"]["exposure_node_hours"]
assert len(NODE_HOURS) == TAU, "tau disagrees with the node schedule"

# Exposure node k spans [NODE_HOURS[k-1], +24). Covariate window for node k is the
# 24h immediately before it. Both come from config; neither is written down here.
_W = COV["windows"]
COV_WINDOWS = {1: (_W["baseline_L0"]["start_h"], _W["baseline_L0"]["end_h"]),
               2: (_W["node_2"]["start_h"], _W["node_2"]["end_h"]),
               3: (_W["node_3"]["start_h"], _W["node_3"]["end_h"])}
EXPO_WINDOWS = {i + 1: (h, h + 24) for i, h in enumerate(NODE_HOURS)}

# ---------------------------------------------------------------------------
# THE ORDERING GUARD.  L_t -> A_t, asserted, not commented.
#
# This is the one invariant whose violation would not announce itself: the frame
# would build, the model would fit, and the estimate would be conditioned on the
# future. It fires at import so a covariate window edited into the wrong place
# stops the run before any data is read.
# ---------------------------------------------------------------------------
for _k in sorted(COV_WINDOWS):
    _ws, _we = COV_WINDOWS[_k]
    _ns, _ne = EXPO_WINDOWS[_k]
    assert _we <= _ns, (
        f"L_t -> A_t VIOLATED at node {_k}: the covariate window [{_ws}, {_we}) ends "
        f"after the exposure node [{_ns}, {_ne}) begins. Covariates would be measured "
        f"after the treatment they are supposed to confound. Fix "
        f"config/lmtp_design.json covariates.windows.")
    print(f"  node {_k}: covariates [{_ws:+3d}h, {_we:+3d}h)  ->  exposure "
          f"[{_ns:+3d}h, {_ne:+3d}h)   ordering OK")

WINDOW_SPAN_H = (min(w[0] for w in COV_WINDOWS.values()),
                 max(w[1] for w in COV_WINDOWS.values()))
EXPOSURE_WINDOW_H = max(NODE_HOURS) + 24

# --- pinned measurement parameters -----------------------------------------
FIO2_LOOKBACK_H = DEFS["pf_ratio"]["fio2_lookback_hours"]
SPO2_CEILING = DEFS["sf_ratio"]["spo2_ceiling"]
CREAT_LOOKBACK_H = DEFS["creatinine"]["creatinine_lookback_hours"]
NEE_COEF = DEFS["nee"]["coefficients"]
NEE_HOLD_H = DEFS["nee"]["hold_hours"]
NEE_DRUGS = list(NEE_COEF)
INOTROPES = DEFS["inotrope"]["sources"][0]["categories"]
CCI_ENTERED = DEFS["cci_components"]["entered"]

# variable name -> lab_category. The keys are our column names, the values are CLIF's.
LAB_VARS = {"creatinine": "creatinine", "bun": "bun", "potassium": "potassium",
            "lactate": "lactate", "bicarbonate": "bicarbonate",
            "ph_arterial": "ph_arterial", "ph_venous": "ph_venous",
            "po2_arterial": "po2_arterial"}

BASELINE_L0 = COV["baseline_L0"]
TIME_VARYING = COV["time_varying_Lt"]
OUTCOME_GRID_D = design["time"]["outcome_grid_days"]

print(f"\nprotocol: tau={TAU}, nodes {NODE_HOURS}, outcome grid {OUTCOME_GRID_D}d")
print(f"  fio2 lookback {FIO2_LOOKBACK_H}h, creatinine lookback {CREAT_LOOKBACK_H}h, "
      f"spo2 ceiling {SPO2_CEILING}")
print(f"  L_0 has {len(BASELINE_L0)} terms, L_t has {len(TIME_VARYING)} per node")

_kw = dict(
    data_directory=config["data_directory"],
    filetype=config["filetype"],
    timezone=config["timezone"],
    output_directory=str(REPO_ROOT / "output"),
)

# Covariate outlier bounds come from clifpy's own shipped config, NOT from the
# vendored config/outlier_config.json. That file is byte-pinned at SHA ee4774b and
# has no bounds for po2_arterial, ph_*, bicarbonate, potassium, bun or creatinine;
# adding them would break the vendor contract and the dose reconciliation at once.
_clifpy_outliers = yaml.safe_load(
    (Path(clifpy.__file__).parent / "schemas" / "outlier_config.yaml").read_text())


def _bounds(table: str, category: str):
    """Look up (lo, hi) for one category in clifpy's shipped outlier config.

    The file nests bounds under the table's VALUE COLUMN, e.g.
    tables -> labs -> lab_value_numeric -> creatinine -> {min, max}, so the value
    column is searched rather than assumed. Returns (None, None) when clifpy does
    not bound this variable, so the caller reports the gap instead of applying no
    bound and believing it applied one.
    """
    node = _clifpy_outliers.get("tables", _clifpy_outliers).get(table, {})
    candidates = [node]
    candidates += [v for v in node.values() if isinstance(v, dict)]
    for cand in candidates:
        e = cand.get(category)
        if isinstance(e, dict):
            return e.get("min", e.get("lower")), e.get("max", e.get("upper"))
        if isinstance(e, (list, tuple)) and len(e) == 2:
            return e[0], e[1]
    return None, None


def _apply_bounds(df, value_col, table, category, label=None):
    """Null out-of-range values in place and report how many went. Never drops rows."""
    lo, hi = _bounds(table, category)
    label = label or category
    if lo is None and hi is None:
        print(f"    {label:<16} NO BOUND in clifpy config  (left unbounded)")
        return df
    v = df[value_col]
    bad = v.notna() & ~v.between(lo if lo is not None else -np.inf,
                                 hi if hi is not None else np.inf)
    if bad.any():
        print(f"    {label:<16} [{lo}, {hi}]  {int(bad.sum()):,} of {len(df):,} nulled")
    df.loc[bad, value_col] = np.nan
    return df


# %%
# ---------------------------------------------------------------------------
# Stage 0: Load 01's artifacts
# ---------------------------------------------------------------------------
def stage_0_load(repo_root):
    """Read the cohort, dose series, and block map written by 01.

    The block map is why 01 persists it: every raw CLIF table is keyed on
    hospitalization_id, and nothing downstream of 01 carries one. Re-deriving it
    here with a second stitch_encounters call would be an independent claim about
    which hospitalizations belong to which block, and two scripts that disagree
    about that disagree about the cohort.
    """
    phi = repo_root / "output" / "intermediate_phi"
    cohort = pd.read_parquet(phi / "cohort.parquet")
    dose_series = pd.read_parquet(phi / "dose_series.parquet")
    block_map = pd.read_parquet(phi / "block_map.parquet")

    print(f"  cohort       {len(cohort):,} blocks")
    print(f"  dose series  {len(dose_series):,} records")
    print(f"  block map    {len(block_map):,} hospitalizations -> "
          f"{block_map['encounter_block'].nunique():,} blocks")

    assert set(block_map["encounter_block"]) == set(cohort["encounter_block"]), (
        "block map does not cover the cohort exactly; re-run 01")
    assert cohort["encounter_block"].is_unique, "cohort is not one row per block"
    return cohort, dose_series, block_map


# %%
cohort, dose_series, block_map = stage_0_load(REPO_ROOT)

# t0 per block, and the hospitalization-level key both carrying it.
t0 = cohort[["encounter_block", "crrt_initiation_dttm", "weight_kg"]].copy()
hosp_t0 = block_map.merge(t0, on="encounter_block", how="left", validate="many_to_one")
HOSP_IDS = block_map["hospitalization_id"].tolist()


def _to_hours(df, dttm_col):
    """Attach hours_from_init and trim to the union of every covariate window.

    One helper because every observation table needs exactly this: map to the block,
    subtract t0, and discard everything outside [-24h, +48h). Trimming here rather
    than per variable keeps the frames small enough to hold several at once.
    """
    d = df.merge(hosp_t0, on="hospitalization_id", how="inner", validate="many_to_one")
    d["hours_from_init"] = ((d[dttm_col] - d["crrt_initiation_dttm"])
                            .dt.total_seconds() / 3600)
    lo, hi = WINDOW_SPAN_H
    return d[d["hours_from_init"].between(lo, hi, inclusive="left")].copy()


# %%
# ---------------------------------------------------------------------------
# Stage 1: The node skeleton
# ---------------------------------------------------------------------------
def stage_1_node_skeleton(cohort):
    """One row per (block, exposure node), with whether the block is still at risk.

    `node_status` is the point of this stage. Whether a node holds a real exposure,
    a liberation zero, or nothing at all is a three-way distinction that decides how
    A_t is filled, and it deserves to be a column that can be counted rather than an
    accident of which groupby happened to return an empty group.
    """
    print(f"  stage 1 input: {len(cohort):,} blocks")

    # block_discharge_dttm, i.e. when the block ENDED however it ended. NOT
    # death_dttm_by_discharge, which 01 masks to deaths and is null for everyone
    # discharged alive. Under a design whose competing event IS discharge alive, a
    # patient who walks out on day 9 stops being at risk on day 9.
    assert cohort["block_discharge_dttm"].notna().all(), "a block has no discharge time"
    base = cohort[["encounter_block", "crrt_initiation_dttm",
                   "block_discharge_dttm"]].copy()
    n_survivors = int(cohort["death_dttm_by_discharge"].isna().sum())
    print(f"  block end from block_discharge_dttm, complete for all {len(cohort):,} "
          f"blocks ({n_survivors:,} of them discharged alive, so death-timed columns "
          f"are null there by design)")

    rows = []
    for k, (s_h, e_h) in EXPO_WINDOWS.items():
        r = base.copy()
        r["node"] = k
        r["node_start_dttm"] = r["crrt_initiation_dttm"] + pd.Timedelta(hours=s_h)
        r["node_end_dttm"] = r["crrt_initiation_dttm"] + pd.Timedelta(hours=e_h)
        # At risk means the block had not already ended when the node opened. A block
        # ending exactly at the boundary is over, so the comparison is strict.
        r["at_risk"] = r["block_discharge_dttm"] > r["node_start_dttm"]
        rows.append(r.drop(columns=["crrt_initiation_dttm", "block_discharge_dttm"]))
    skel = pd.concat(rows, ignore_index=True)

    print("  at risk at node open:")
    for k in sorted(EXPO_WINDOWS):
        n = int(skel.loc[skel["node"] == k, "at_risk"].sum())
        print(f"    node {k} ({EXPO_WINDOWS[k][0]:>2}-{EXPO_WINDOWS[k][1]:>2}h): "
              f"{n:,} of {len(cohort):,} ({100 * n / len(cohort):.0f}%)")
    assert len(skel) == len(cohort) * TAU, "skeleton is not blocks x nodes"
    return skel


# %%
skeleton = stage_1_node_skeleton(cohort)


# %%
# ---------------------------------------------------------------------------
# Stage 2: Exposure  —  A_t under S1 and S2
# ---------------------------------------------------------------------------
def stage_2_exposure(skeleton, dose_series, cohort):
    """Time-weighted mean delivered dose per node, under the S1 and S2 conventions.

    Built on an HOURLY BIN reconstruction, which is what makes the ladder coherent:
    CRRT charting here is hourly and regular, so each hour of a node is either
    charted or not, and the three rungs differ only in what goes in the denominator.

        S1  sum(charted bin doses) / n charted bins          (primary, upper bound)
        S2  sum(charted bin doses) / bins from first to last (gaps inside the span
                                                              count as downtime)
        S3  sum(charted bin doses) / 24                      (blocked, see todo 2e)

    A node with NO charted record at all, in a block still at risk, is liberation
    and takes exposure 0 unshifted. That is the empty-node rule, and it is the seam
    where two separately-settled decisions meet: S1's denominator is undefined at
    zero charted hours, while the discontinuation decision says liberation is a
    dose of zero. The empty node is the liberation case.
    """
    print(f"  stage 2 input: {len(dose_series):,} dose records")

    ds = dose_series.copy()
    ds["hour_bin"] = np.floor(ds["hours_from_init"]).astype(int)
    # One value per charted hour. Several records in an hour is rare here but must
    # not let a densely charted hour outvote a sparsely charted one.
    hourly = (ds.groupby(["encounter_block", "hour_bin"], as_index=False)["dose_ml_kg_hr"]
                .mean())
    print(f"  {len(hourly):,} charted block-hours across "
          f"{hourly['encounter_block'].nunique():,} blocks")

    out = []
    for k, (s_h, e_h) in EXPO_WINDOWS.items():
        w = hourly[hourly["hour_bin"].between(s_h, e_h - 1)]
        g = w.groupby("encounter_block")["dose_ml_kg_hr"].agg(
            total="sum", n_charted_hours="count", first_bin="idxmin", last_bin="idxmax")
        # idxmin/idxmax give row labels; translate back to bin numbers.
        span = w.groupby("encounter_block")["hour_bin"].agg(["min", "max"])
        g = g.drop(columns=["first_bin", "last_bin"]).join(span)
        g["node"] = k
        out.append(g.reset_index())

    ex = pd.concat(out, ignore_index=True)
    ex["observed_span_h"] = ex["max"] - ex["min"] + 1
    ex["gap_hours"] = ex["observed_span_h"] - ex["n_charted_hours"]
    ex["a_s1"] = ex["total"] / ex["n_charted_hours"]
    ex["a_s2"] = ex["total"] / ex["observed_span_h"]
    ex = ex.drop(columns=["total", "min", "max"])

    e = skeleton.merge(ex, on=["encounter_block", "node"], how="left")
    e["n_charted_hours"] = e["n_charted_hours"].fillna(0).astype(int)
    e["gap_hours"] = e["gap_hours"].fillna(0).astype(int)

    empty = e["n_charted_hours"] == 0
    e["node_status"] = np.where(~e["at_risk"], "post_event",
                        np.where(empty, "liberated", "on_crrt"))
    # The empty-node rule. Liberation is an exposure of zero that the policy leaves
    # alone, because 0 - delta < floor for every delta in the ladder.
    lib = e["node_status"] == "liberated"
    e.loc[lib, ["a_s1", "a_s2"]] = 0.0

    print("\n  node status:")
    for k in sorted(EXPO_WINDOWS):
        sub = e[e["node"] == k]
        counts = sub["node_status"].value_counts()
        bits = "  ".join(f"{s} {counts.get(s, 0):,}"
                         for s in ("on_crrt", "liberated", "post_event"))
        print(f"    node {k}: {bits}")

    print("\n  exposure, S1 primary (blocks on CRRT only):")
    for k in sorted(EXPO_WINDOWS):
        s = e.loc[(e["node"] == k) & (e["node_status"] == "on_crrt"), "a_s1"]
        if len(s):
            print(f"    node {k}: n {len(s):,}  median {s.median():.1f}  "
                  f"p5 {s.quantile(.05):.1f}  p95 {s.quantile(.95):.1f}  "
                  f"max {s.max():.1f}")

    on = e["node_status"] == "on_crrt"
    agree = np.isclose(e.loc[on, "a_s1"], e.loc[on, "a_s2"])
    print(f"\n  S1 == S2 for {int(agree.sum()):,} of {int(on.sum()):,} charted "
          f"node-rows ({100 * agree.mean():.1f}%), i.e. no uncharted gap")
    n1 = e[(e["node"] == 1) & on]
    if len(n1):
        a1 = np.isclose(n1["a_s1"], n1["a_s2"])
        print(f"    node 1 alone: {int(a1.sum()):,} of {len(n1):,} ({100 * a1.mean():.1f}%)")

    assert e.loc[e["node_status"] != "post_event", ["a_s1", "a_s2"]].notna().all().all(), (
        "an at-risk node has no exposure value; the empty-node rule did not run")
    return e


# %%
exposure = stage_2_exposure(skeleton, dose_series, cohort)


# %%
# ---------------------------------------------------------------------------
# Stage 3: Observation tables  —  load once, summarise many times
# ---------------------------------------------------------------------------
def stage_3_load_observations(hosp_ids):
    """Load every raw table the covariates need, trimmed to the window and cohort.

    Loaded once and kept tidy, because baseline and both time-varying nodes are the
    same summary applied to three different windows. Nothing here is summarised yet.
    """
    print(f"  loading for {len(hosp_ids):,} hospitalizations, window "
          f"[{WINDOW_SPAN_H[0]:+d}h, {WINDOW_SPAN_H[1]:+d}h)")

    labs = Labs.from_file(
        **_kw,
        columns=["hospitalization_id", "lab_result_dttm", "lab_category",
                 "lab_value_numeric"],
        filters={"lab_category": list(LAB_VARS.values())},
    ).df
    labs = _to_hours(labs, "lab_result_dttm")
    print(f"  labs                {len(labs):,} results in window")

    vit = Vitals.from_file(
        **_kw,
        columns=["hospitalization_id", "recorded_dttm", "vital_category", "vital_value"],
        filters={"vital_category": ["spo2", "weight_kg"]},
    ).df
    spo2 = _to_hours(vit[vit["vital_category"] == "spo2"], "recorded_dttm")
    print(f"  spo2                {len(spo2):,} readings in window")

    resp = RespiratorySupport.from_file(
        **_kw,
        columns=["hospitalization_id", "recorded_dttm", "device_category", "fio2_set"],
    ).df
    resp = resp[resp["hospitalization_id"].isin(set(hosp_ids))]
    resp = _to_hours(resp, "recorded_dttm")
    print(f"  respiratory support {len(resp):,} rows in window")

    meds = MedicationAdminContinuous.from_file(
        **_kw,
        columns=["hospitalization_id", "admin_dttm", "med_category", "med_dose",
                 "med_dose_unit", "mar_action_category"],
        filters={"med_category": NEE_DRUGS + INOTROPES},
    ).df
    # Meds are converted BEFORE trimming, since the unit converter needs the weight
    # history and a dose charted just outside the window still sets the rate inside it.
    meds_raw = meds[meds["hospitalization_id"].isin(set(hosp_ids))].copy()
    print(f"  continuous meds     {len(meds_raw):,} rows for the cohort")

    diag = HospitalDiagnosis.from_file(**_kw).df
    diag = diag[diag["hospitalization_id"].isin(set(hosp_ids))]
    print(f"  diagnoses           {len(diag):,} rows")

    pat = Patient.from_file(**_kw, columns=["patient_id", "sex_category"]).df

    print("\n  outlier bounds (clifpy schemas/outlier_config.yaml):")
    for var, cat in LAB_VARS.items():
        m = labs["lab_category"] == cat
        sub = labs.loc[m, ["lab_value_numeric"]]
        if len(sub):
            labs.loc[m, "lab_value_numeric"] = _apply_bounds(
                sub.copy(), "lab_value_numeric", "labs", cat, var)["lab_value_numeric"]
    spo2 = _apply_bounds(spo2, "vital_value", "vitals", "spo2")
    resp = _apply_bounds(resp, "fio2_set", "respiratory_support", "fio2_set")

    return dict(labs=labs, spo2=spo2, resp=resp, meds=meds_raw, diag=diag,
                pat=pat, vitals_weight=vit[vit["vital_category"] == "weight_kg"])


# %%
obs = stage_3_load_observations(HOSP_IDS)


# %%
# ---------------------------------------------------------------------------
# Window summarisation  —  the one path every simple covariate takes
# ---------------------------------------------------------------------------
def _window_mask(hours, node):
    lo, hi = COV_WINDOWS[node]
    return hours.between(lo, hi, inclusive="left")


def summarise_windows(frame, value_col, how, name):
    """Apply one summary rule to one variable across all three covariate windows.

    Returns long: encounter_block, node, <name>. `how` is read from the config, so a
    rule change is a config edit and never a code edit.
    """
    out = []
    for node in sorted(COV_WINDOWS):
        w = frame[_window_mask(frame["hours_from_init"], node)]
        w = w[w[value_col].notna()]
        if not len(w):
            continue
        if how in ("min", "max", "mean", "median"):
            g = w.groupby("encounter_block")[value_col].agg(how)
        elif how.startswith("last"):
            g = (w.sort_values("hours_from_init")
                  .groupby("encounter_block")[value_col].last())
        elif how.startswith("any"):
            g = w.groupby("encounter_block")[value_col].max().astype(float)
        else:
            raise ValueError(f"unknown summary rule {how!r} for {name}")
        out.append(g.rename(name).reset_index().assign(node=node))
    if not out:
        return pd.DataFrame(columns=["encounter_block", "node", name])
    return pd.concat(out, ignore_index=True)


# %%
# ---------------------------------------------------------------------------
# Stage 4: Covariates
# ---------------------------------------------------------------------------
def _pf_and_sf(labs, spo2, resp):
    """P/F and S/F by asof-pairing each oxygenation reading to a preceding FiO2.

    The config's rule stated literally: for every PaO2 in the window, take the most
    recent non-null fio2_set at or before it, searching back at most FIO2_LOOKBACK_H,
    then take the minimum ratio in the window. `merge_asof` with a tolerance IS that
    rule, which is why it is used instead of a hand-rolled join.

    S/F is a SEPARATE column, never blended into P/F. The two are on different
    scales, and blending would make the variable's meaning depend on whether an
    arterial gas happened to be drawn, which is itself informative.
    """
    fio2 = (resp.loc[resp["fio2_set"].notna(),
                     ["encounter_block", "recorded_dttm", "fio2_set"]]
                .sort_values("recorded_dttm"))

    def pair(num_df, dttm_col, value_col, out_name):
        left = num_df.sort_values(dttm_col)
        if not len(left) or not len(fio2):
            return pd.DataFrame(columns=["encounter_block", "hours_from_init", out_name])
        m = pd.merge_asof(
            left, fio2, left_on=dttm_col, right_on="recorded_dttm",
            by="encounter_block", direction="backward",
            tolerance=pd.Timedelta(hours=FIO2_LOOKBACK_H))
        n_pairable = int(m["fio2_set"].notna().sum())
        print(f"    {out_name}: {n_pairable:,} of {len(m):,} readings paired to an "
              f"FiO2 within {FIO2_LOOKBACK_H}h ({100 * n_pairable / max(len(m), 1):.0f}%)")
        m = m[m["fio2_set"].notna() & (m["fio2_set"] > 0)]
        m[out_name] = m[value_col] / m["fio2_set"]
        return m[["encounter_block", "hours_from_init", out_name]]

    po2 = labs[labs["lab_category"] == "po2_arterial"].copy()
    pf = pair(po2, "lab_result_dttm", "lab_value_numeric", "pf_ratio")

    # S/F only where SpO2 is still informative about PaO2. Above the ceiling the
    # dissociation curve is flat and the ratio just reports the FiO2.
    s = spo2[spo2["vital_value"] <= SPO2_CEILING].copy()
    print(f"    spo2 <= {SPO2_CEILING}: {len(s):,} of {len(spo2):,} readings usable")
    sf = pair(s, "recorded_dttm", "vital_value", "sf_ratio")

    return (summarise_windows(pf, "pf_ratio", DEFS["pf_ratio"]["summary"], "pf_ratio"),
            summarise_windows(sf, "sf_ratio", DEFS["sf_ratio"]["summary"], "sf_ratio"))


def _nee(meds, vitals_weight):
    """Norepinephrine-equivalent dose as the max of a summed step function.

    Two things this does that a row-wise maximum does not:

    1. SUMS CONCURRENT DRUGS. A patient on moderate doses of three pressors is
       sicker than one on a slightly higher dose of a single agent, and a row-wise
       max cannot see the difference.
    2. HOLDS RATES FORWARD between records, so a drug charted at 01:00 still counts
       at 01:30 when a second drug is charted. Held at most NEE_HOLD_H hours: with
       a median inter-record interval near one hour, a longer silence means the
       infusion ended without a stop being charted.

    `mar_action_category == "stop"` sets the rate to zero, which is why the step
    function is honest rather than an assumption.
    """
    if not len(meds):
        return pd.DataFrame(columns=["encounter_block", "node", "nee", "inotrope"])

    conv, _ = convert_dose_units_by_med_category(
        meds, vitals_df=vitals_weight,
        preferred_units={**{d: "mcg/kg/min" for d in NEE_DRUGS + INOTROPES},
                         "vasopressin": "u/min"})
    conv = conv.rename(columns={"med_dose_converted": "dose_std"})
    bad_unit = conv["dose_std"].isna() & conv["med_dose"].notna()
    if bad_unit.any():
        print(f"    {int(bad_unit.sum()):,} of {len(conv):,} med rows could not be "
              f"converted to a standard unit and are dropped")
    conv = _to_hours(conv, "admin_dttm")

    # A stop is a rate of zero, not a missing value.
    conv.loc[conv["mar_action_category"].str.lower().eq("stop"), "dose_std"] = 0.0
    conv = conv[conv["dose_std"].notna()]
    conv["hour_bin"] = np.floor(conv["hours_from_init"]).astype(int)

    # Last charted rate per drug-hour, then hold forward across the hourly grid.
    per = (conv.sort_values("hours_from_init")
               .groupby(["encounter_block", "med_category", "hour_bin"], as_index=False)
               ["dose_std"].last())
    grid = pd.MultiIndex.from_product(
        [per["encounter_block"].unique(),
         per["med_category"].unique(),
         range(WINDOW_SPAN_H[0], WINDOW_SPAN_H[1])],
        names=["encounter_block", "med_category", "hour_bin"]).to_frame(index=False)
    per = grid.merge(per, on=["encounter_block", "med_category", "hour_bin"], how="left")
    per["dose_std"] = (per.sort_values("hour_bin")
                          .groupby(["encounter_block", "med_category"])["dose_std"]
                          .ffill(limit=NEE_HOLD_H))
    per = per[per["dose_std"].notna()]

    press = per[per["med_category"].isin(NEE_DRUGS)].copy()
    press["nee_contrib"] = press["dose_std"] * press["med_category"].map(NEE_COEF)
    nee_hourly = (press.groupby(["encounter_block", "hour_bin"], as_index=False)
                       ["nee_contrib"].sum().rename(columns={"nee_contrib": "nee"}))
    nee_hourly["hours_from_init"] = nee_hourly["hour_bin"]

    ino = per[per["med_category"].isin(INOTROPES)].copy()
    ino["inotrope"] = (ino["dose_std"] > 0).astype(float)
    ino_hourly = (ino.groupby(["encounter_block", "hour_bin"], as_index=False)
                     ["inotrope"].max())
    ino_hourly["hours_from_init"] = ino_hourly["hour_bin"]

    a = summarise_windows(nee_hourly, "nee", DEFS["nee"]["summary"].split()[0], "nee")
    b = summarise_windows(ino_hourly, "inotrope", "any", "inotrope")
    nz = a.loc[a["nee"] > 0, "nee"]
    print(f"    nee       {len(a):,} block-nodes, {len(nz):,} on any vasopressor; "
          f"median {nz.median():.3f}  p95 {nz.quantile(.95):.3f} mcg/kg/min-equivalent"
          if len(nz) else "    nee: no vasopressor exposure found")
    print(f"    inotrope  {len(b):,} block-nodes, "
          f"{int(b['inotrope'].sum()):,} with an inotrope running")
    return a.merge(b, on=["encounter_block", "node"], how="outer")


def _imv(resp):
    """IMV in effect at each covariate window's END boundary.

    A point-in-time state at the boundary, not "any IMV during the window": what
    precedes the next dose decision is the patient's state when that decision is
    made. Nulls in device_category are skipped rather than carried, so the answer is
    the last KNOWN device.
    """
    known = (resp.loc[resp["device_category"].notna(),
                      ["encounter_block", "hours_from_init", "device_category"]]
                 .sort_values("hours_from_init"))
    known["imv_status"] = (known["device_category"].str.upper() == "IMV").astype(float)
    out = []
    for node in sorted(COV_WINDOWS):
        boundary = COV_WINDOWS[node][1]
        w = known[known["hours_from_init"] < boundary]
        if not len(w):
            continue
        g = w.groupby("encounter_block")["imv_status"].last()
        out.append(g.reset_index().assign(node=node))
    return (pd.concat(out, ignore_index=True) if out
            else pd.DataFrame(columns=["encounter_block", "node", "imv_status"]))


def _ph(labs):
    """pH, arterial preferred, venous only where the window holds no arterial result.

    Indicated, never corrected. Venous pH runs roughly 0.03-0.05 below arterial, and
    a fixed offset would be worse than an indicator because the true offset varies
    with perfusion, which is exactly what is abnormal in this cohort.
    """
    art = labs[labs["lab_category"] == "ph_arterial"]
    ven = labs[labs["lab_category"] == "ph_venous"]
    how = DEFS["ph"]["summary"]
    a = summarise_windows(art, "lab_value_numeric", how, "ph").assign(ph_source="arterial")
    v = summarise_windows(ven, "lab_value_numeric", how, "ph").assign(ph_source="venous")
    # Arterial wins wherever it exists for that block-node.
    keys = set(map(tuple, a[["encounter_block", "node"]].values.tolist()))
    v = v[~v[["encounter_block", "node"]].apply(tuple, axis=1).isin(keys)] if len(v) else v
    both = pd.concat([a, v], ignore_index=True)
    if len(both):
        print(f"    pH: {int((both['ph_source'] == 'arterial').sum()):,} block-nodes "
              f"arterial, {int((both['ph_source'] == 'venous').sum()):,} venous fallback")
    return both


def _sofa_by_window(block_map, cohort, repo_root):
    """Non-renal SOFA per covariate window, for the pre-specified fallback set.

    Computed even though the primary parameterisation does not use it: the fallback
    is triggered by a positivity criterion at node 2, and having the columns already
    in the frame means 03 can switch without rebuilding the dataset.

    fill_na_scores_with_zero=False deliberately. clifpy's default scores a missing
    component as 0, i.e. normal, which at sparse nodes biases severity downward
    exactly where the data are thinnest.
    """
    base = block_map.merge(cohort[["encounter_block", "crrt_initiation_dttm"]],
                           on="encounter_block", how="left")
    frames = []
    for node in sorted(COV_WINDOWS):
        lo, hi = COV_WINDOWS[node]
        c = base.copy()
        c["start_dttm"] = c["crrt_initiation_dttm"] + pd.Timedelta(hours=lo)
        c["end_dttm"] = c["crrt_initiation_dttm"] + pd.Timedelta(hours=hi)
        cdf = pl.from_pandas(c[["hospitalization_id", "encounter_block",
                                "start_dttm", "end_dttm"]])
        try:
            s = compute_sofa_polars(
                data_directory=config["data_directory"], cohort_df=cdf,
                filetype=config["filetype"], id_name="encounter_block",
                fill_na_scores_with_zero=False, timezone=config["timezone"],
            ).to_pandas()
        except Exception as exc:                       # noqa: BLE001
            print(f"    SOFA node {node}: FAILED ({type(exc).__name__}: {exc})")
            continue
        keep = [c for c in s.columns if c.startswith("sofa_") or c == "encounter_block"]
        s = s[keep].copy()
        nonrenal = [c for c in s.columns
                    if c.startswith("sofa_") and c not in ("sofa_total", "sofa_renal")]
        # NULL unless every non-renal component is present. Summing what happens to be
        # available would let a missing component contribute 0 and read as normal,
        # which is the exact bias fill_na_scores_with_zero=False was set to avoid.
        s["sofa_nonrenal"] = s[nonrenal].sum(axis=1).where(s[nonrenal].notna().all(axis=1))
        s["node"] = node
        frames.append(s)
        print(f"    SOFA node {node}: {len(s):,} blocks, non-renal median "
              f"{s['sofa_nonrenal'].median():.1f}")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def stage_4_covariates(obs, block_map, cohort, repo_root):
    """Every covariate, in long (block, node) form.

    Deliberately one function per hard covariate and one shared path for the simple
    ones, because the hard ones are hard for different reasons: P/F needs a pairing,
    NEE needs a step function, IMV needs a boundary state, pH needs a preference
    order. Forcing those through a single generic summariser would have hidden four
    different judgements behind one signature.
    """
    labs, spo2, resp = obs["labs"], obs["spo2"], obs["resp"]
    frames = []

    print("  oxygenation:")
    pf, sf = _pf_and_sf(labs, spo2, resp)
    frames += [pf, sf]

    print("  simple labs:")
    for var, cat in LAB_VARS.items():
        if var in ("po2_arterial", "ph_arterial", "ph_venous"):
            continue                                   # handled by P/F and pH
        how = DEFS[var]["summary"]
        sub = labs[labs["lab_category"] == cat]
        f = summarise_windows(sub, "lab_value_numeric", how, var)
        print(f"    {var:<12} {how:<22} {len(f):,} block-nodes")
        frames.append(f)

    print("  pH:")
    frames.append(_ph(labs))

    print("  vasoactives:")
    frames.append(_nee(obs["meds"], obs["vitals_weight"]))

    print("  ventilation:")
    imv = _imv(resp)
    print(f"    imv_status  {len(imv):,} block-nodes")
    frames.append(imv)

    print("  SOFA (fallback parameterisation only):")
    sofa = _sofa_by_window(block_map, cohort, repo_root)
    if len(sofa):
        frames.append(sofa)

    long = None
    for f in frames:
        if not len(f):
            continue
        long = f if long is None else long.merge(
            f, on=["encounter_block", "node"], how="outer")
    return long


# %%
cov_long = stage_4_covariates(obs, block_map, cohort, REPO_ROOT)


# %%
# ---------------------------------------------------------------------------
# Stage 5: Static covariates  —  demographics and comorbidity
# ---------------------------------------------------------------------------
def stage_5_static(cohort, obs, block_map):
    """Age, sex, weight, and the six entered Charlson components.

    CCI is computed at BLOCK grain, matching the ESRD exclusion: a comorbidity coded
    on any hospitalization in a stitched block belongs to the patient for the whole
    block. All 17 components are computed because one call returns them; six are
    entered, and cci_score is retained for description but is NOT a model term.
    """
    print(f"  stage 5 input: {len(cohort):,} blocks")
    stat = cohort[["encounter_block", "patient_id", "age_at_admission",
                   "weight_kg"]].copy()
    stat = stat.merge(obs["pat"], on="patient_id", how="left", validate="many_to_one")
    stat = stat.rename(columns={"sex_category": "sex"})
    print(f"  sex missing: {int(stat['sex'].isna().sum()):,}")

    cci = calculate_cci(obs["diag"], hierarchy=True)
    cci = cci.to_pandas() if hasattr(cci, "to_pandas") else pd.DataFrame(cci)
    cci = cci.merge(block_map, on="hospitalization_id", how="inner")
    comp_cols = [c for c in cci.columns
                 if c not in ("hospitalization_id", "encounter_block")]
    # Block grain: any hospitalization in the block coding the condition sets it.
    cci_block = cci.groupby("encounter_block")[comp_cols].max().reset_index()

    # Diabetes is lumped from clifpy's two mutually exclusive columns (hierarchy=True),
    # so the maximum is simply "either". Built before the membership check below, so
    # the lumped name is a real column by the time it is verified.
    lump = DEFS["cci_components"].get("diabetes_lumping", [])
    if lump and all(c in cci_block.columns for c in lump):
        cci_block["diabetes"] = cci_block[lump].max(axis=1)

    # The config names components exactly as clifpy emits them, so this is an exact
    # membership test rather than a fuzzy match. A fuzzy match would keep passing while
    # silently entering the wrong column the day clifpy renames one.
    print(f"  CCI: {len(comp_cols)} columns computed, {len(CCI_ENTERED)} entered")
    unmatched = [c for c in CCI_ENTERED if c not in cci_block.columns]
    assert not unmatched, (
        f"config names components clifpy does not emit: {unmatched}. "
        f"Available: {sorted(cci_block.columns)}")

    stat = stat.merge(cci_block, on="encounter_block", how="left")
    for c in list(comp_cols) + ["diabetes"]:
        if c in stat.columns and c != "cci_score":
            stat[c] = stat[c].fillna(0)
    for c in CCI_ENTERED:
        print(f"    {c:<32} {100 * stat[c].mean():5.1f}% prevalence")
    return stat, comp_cols


# %%
static, cci_cols = stage_5_static(cohort, obs, block_map)


# %%
# ---------------------------------------------------------------------------
# Stage 6: Outcomes, competing event, censoring
# ---------------------------------------------------------------------------
def stage_6_outcomes(cohort, grid_days):
    """Y (in-hospital death), D (discharge alive), C (censoring) on the outcome grid.

    Monotone by construction: once an event has happened it stays happened, which is
    what `lmtp`'s event_locf convention expects.

    Censoring is structurally EMPTY here, and that is worth saying out loud rather
    than leaving as an all-ones column nobody examines. Mortality is in-hospital by
    construction and every block ends in exactly one of death or discharge alive, so
    the two indicators partition every outcome and nothing is lost to follow-up. C is
    emitted anyway because `lmtp` requires it for survival outcomes.
    """
    print(f"  stage 6 input: {len(cohort):,} blocks")
    o = cohort[["encounter_block", "days_to_block_end", "died_in_hospital",
                "mortality_30d"]].copy()
    # Time to the END of the block, whichever event ended it. 01 now supplies this
    # directly as days_to_block_end; days_to_death is the other, death-only column.
    o["days_to_end"] = o["days_to_block_end"]
    assert o["days_to_end"].notna().all(), "a block has no time to block end"
    n_impossible = int((o["days_to_end"] <= 0).sum())
    if n_impossible:
        print(f"  WARNING: {n_impossible} blocks end at or before CRRT initiation "
              f"(the known timestamp defect, todo 2e). They are post-event at every "
              f"node and contribute no exposure.")

    for d in grid_days:
        reached = o["days_to_end"] <= d
        o[f"y_d{d}"] = (reached & o["died_in_hospital"]).astype(int)
        o[f"d_d{d}"] = (reached & ~o["died_in_hospital"]).astype(int)
        o[f"c_d{d}"] = 1
    for d in grid_days:
        print(f"    day {d:>2}: died {o[f'y_d{d}'].sum():>5,}   "
              f"discharged alive {o[f'd_d{d}'].sum():>5,}   "
              f"still in hospital {int((~(o[f'y_d{d}'] | o[f'd_d{d}']).astype(bool)).sum()):>5,}")

    last = grid_days[-1]
    assert (o[f"y_d{last}"] == o["mortality_30d"].astype(int)).all(), (
        f"y_d{last} disagrees with the cohort's own mortality flag")
    for a, b in zip(grid_days, grid_days[1:]):
        assert (o[f"y_d{b}"] >= o[f"y_d{a}"]).all(), "outcome is not monotone"
        assert (o[f"d_d{b}"] >= o[f"d_d{a}"]).all(), "competing event is not monotone"
    assert not ((o[f"y_d{last}"] == 1) & (o[f"d_d{last}"] == 1)).any(), (
        "a block both died and was discharged alive")
    return o.drop(columns=["days_to_block_end", "days_to_end", "died_in_hospital",
                           "mortality_30d"])


# %%
outcomes = stage_6_outcomes(cohort, OUTCOME_GRID_D)


# %%
# ---------------------------------------------------------------------------
# Stage 7: Assemble the wide frame
# ---------------------------------------------------------------------------
def stage_7_assemble(exposure, cov_long, static, outcomes, cci_cols):
    """Pivot to one row per block, LOCF the time-varying covariates, and fill.

    Missingness is measured BEFORE the fill and reported, because the honest
    description of a covariate is not its value but how often that value was carried
    from somewhere else.
    """
    print(f"  stage 7 input: {len(static):,} blocks")
    cov_vars = [c for c in cov_long.columns
                if c not in ("encounter_block", "node")]

    # Complete (block x node) grid first, so "no row" and "row with NA" become the
    # same thing and every count below is over the same denominator.
    full = pd.MultiIndex.from_product(
        [static["encounter_block"], sorted(COV_WINDOWS)],
        names=["encounter_block", "node"]).to_frame(index=False)
    cl = full.merge(cov_long, on=["encounter_block", "node"], how="left")
    cl = cl.sort_values(["encounter_block", "node"]).reset_index(drop=True)

    # A medication with no record was not being given. Setting these to zero BEFORE
    # the missingness count is deliberate: counting them as missing would report
    # ignorance where the data are in fact informative, and LOCF-ing them would
    # invent vasopressor exposure for patients who had been weaned off it.
    zero_vars = [v for v in COV["missing_values"].get("absence_means_zero", [])
                 if v in cl.columns]
    if zero_vars:
        n_zeroed = int(cl[zero_vars].isna().sum().sum())
        cl[zero_vars] = cl[zero_vars].fillna(0.0)
        print(f"  absence-means-zero: {n_zeroed:,} cells set to 0 across {zero_vars}")

    # Missingness is counted among AT-RISK node-rows only. A post-event node has no
    # covariates because the block had ended, which is structure rather than data
    # quality, and mixing the two would make node 3 look far worse than it is.
    at_risk = (exposure[["encounter_block", "node", "at_risk"]]
               .rename(columns={"node": "node"}))
    clr = cl.merge(at_risk, on=["encounter_block", "node"], how="left")
    miss = []
    for node in sorted(COV_WINDOWS):
        sub = clr[(clr["node"] == node) & clr["at_risk"].fillna(False)]
        denom = max(len(sub), 1)
        for v in cov_vars:
            n_missing = int(sub[v].isna().sum())
            miss.append({"node": node, "variable": v, "n_missing": n_missing,
                         "n_at_risk": len(sub),
                         "pct_missing": round(100 * n_missing / denom, 1)})
    miss_df = pd.DataFrame(miss)
    print("\n  covariate missingness BEFORE LOCF (% of AT-RISK blocks at that node):")
    piv = miss_df.pivot(index="variable", columns="node", values="pct_missing")
    print(piv.to_string(float_format=lambda v: f"{v:5.1f}"))
    print("  at-risk denominators: " + ", ".join(
        f"node {n} = {miss_df.loc[miss_df['node'] == n, 'n_at_risk'].iloc[0]:,}"
        for n in sorted(COV_WINDOWS)))

    # LOCF forward through the nodes, per the Diaz convention for later time points.
    # zero_vars are excluded: they were already resolved by absence-means-zero.
    numeric = [v for v in cov_vars
               if cl[v].dtype.kind in "fiu" and v not in zero_vars]
    before = cl[numeric].isna()
    cl[numeric] = cl.groupby("encounter_block")[numeric].ffill()
    filled = before & cl[numeric].notna()
    print(f"\n  LOCF filled {int(filled.values.sum()):,} cells "
          f"({100 * filled.values.mean():.1f}% of LOCF-eligible covariate cells)")

    # Built as a dict and concatenated once. Inserting ~160 columns one at a time
    # fragments the frame and pandas warns about it.
    cols = {}
    for node in sorted(COV_WINDOWS):
        m = cl["node"] == node
        sub = cl[m].set_index("encounter_block")
        f = filled[m.values].set_index(cl.loc[m, "encounter_block"])
        suffix = "0" if node == 1 else str(node)
        for v in cov_vars:
            cols[f"{v}_{suffix}"] = sub[v].reindex(static["encounter_block"]).values
        for v in numeric:
            cols[f"{v}_{suffix}_locf"] = (
                f[v].reindex(static["encounter_block"]).fillna(False).values)

    for k in sorted(EXPO_WINDOWS):
        e = exposure[exposure["node"] == k].set_index("encounter_block")
        for src, dst in (("a_s1", f"a{k}_s1"), ("a_s2", f"a{k}_s2")):
            cols[dst] = e[src].reindex(static["encounter_block"]).values
        for col in ("n_charted_hours", "gap_hours", "node_status"):
            cols[f"{col}_{k}"] = e[col].reindex(static["encounter_block"]).values

    # Post-window exposure columns, one per outcome period beyond the exposure nodes.
    #
    # lmtp needs trt to be length 1 or exactly tau, and tau is the number of outcome
    # columns (Task.R:125-130), so a 6-period outcome grid needs 6 exposure columns even
    # though the policy only acts on the first three. These are a CONSTANT 0: there is no
    # therapy to intervene on after the exposure window, and 03's shift function returns
    # them untouched, which makes their density ratio exactly 1 and therefore free.
    #
    # Constant, not the last observed dose. Carrying dose forward would assert a therapy
    # that is not being delivered, and would do so most often in the patients who were
    # liberated earliest, i.e. the recovering ones.
    n_expo = len(EXPO_WINDOWS)
    n_periods = len(OUTCOME_GRID_D)
    for k in range(n_expo + 1, n_periods + 1):
        for arm in ("s1", "s2"):
            cols[f"a{k}_{arm}"] = 0.0
        cols[f"node_status_{k}"] = "post_exposure_window"
    if n_periods > n_expo:
        print(f"  emitted {n_periods - n_expo} post-window exposure periods "
              f"(a{n_expo + 1}..a{n_periods}, constant 0, never intervened)")

    wide = pd.concat([static.reset_index(drop=True),
                      pd.DataFrame(cols, index=static.index).reset_index(drop=True)],
                     axis=1)

    # Post-event nodes carry no exposure. lmtp forbids NA in trt, and values after
    # the event never enter the estimand, so they are carried forward from the last
    # observed node purely to satisfy the interface.
    for arm in ("s1", "s2"):
        cols = [f"a{k}_{arm}" for k in sorted(EXPO_WINDOWS)]
        n_before = int(wide[cols].isna().sum().sum())
        wide[cols] = wide[cols].ffill(axis=1)
        wide[cols] = wide[cols].fillna(0.0)
        print(f"  {arm}: filled {n_before:,} post-event exposure cells")

    wide = wide.merge(outcomes, on="encounter_block", how="left", validate="one_to_one")

    trt = [f"a{k}_{a}" for k in range(1, len(OUTCOME_GRID_D) + 1)
           for a in ("s1", "s2")]
    out_cols = [c for c in wide.columns if c[:2] in ("y_", "d_", "c_")]
    assert wide[trt].notna().all().all(), "an exposure column still holds NA"
    assert wide[out_cols].notna().all().all(), "an outcome column holds NA"
    assert len(wide) == len(static), "assembly changed the row count"
    return wide, miss_df


# %%
lmtp_df, missingness = stage_7_assemble(exposure, cov_long, static, outcomes, cci_cols)


# %%
# ---------------------------------------------------------------------------
# Stage 8: Diagnostics and write
# ---------------------------------------------------------------------------
def _code_version(repo_root):
    """git describe, or "unknown" where this is not a git checkout."""
    try:
        return subprocess.run(
            ["git", "-C", str(repo_root), "describe", "--always", "--dirty"],
            capture_output=True, text=True, timeout=10, check=True,
        ).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def stage_8_write(lmtp_df, exposure, missingness, config, design, repo_root):
    """Write the analysis frame and the shareable diagnostics.

    The diagnostics file carries what a coordinating centre needs to see per site
    BEFORE pooling: node availability, how much of each node was actually charted,
    how often a charted record was an explicit zero, and how wide the S1/S2 bracket
    is. A site whose downtime charting differs from everyone else's is then visible
    in the pooled comparison rather than silently shifting its own exposure.
    """
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

    out = lmtp_df.sort_values("encounter_block").reset_index(drop=True)
    out.to_parquet(phi / "lmtp_df.parquet", index=False)
    print(f"    wrote intermediate_phi/lmtp_df.parquet  {len(out):,} x {out.shape[1]}")

    rows = []
    for k in sorted(EXPO_WINDOWS):
        e = exposure[exposure["node"] == k]
        on = e[e["node_status"] == "on_crrt"]
        n = len(e)
        rows += [
            {"metric": "blocks_at_risk", "node": k, "value": float(e["at_risk"].sum())},
            {"metric": "blocks_on_crrt", "node": k, "value": float(len(on))},
            {"metric": "blocks_liberated", "node": k,
             "value": float((e["node_status"] == "liberated").sum())},
            {"metric": "blocks_post_event", "node": k,
             "value": float((e["node_status"] == "post_event").sum())},
            {"metric": "node_availability_pct", "node": k, "value": round(100 * len(on) / n, 1)},
            {"metric": "mean_charted_hours_of_24", "node": k,
             "value": round(float(on["n_charted_hours"].mean()), 2) if len(on) else np.nan},
            {"metric": "mean_gap_hours", "node": k,
             "value": round(float(on["gap_hours"].mean()), 2) if len(on) else np.nan},
            {"metric": "pct_no_gap", "node": k,
             "value": round(100 * float((on["gap_hours"] == 0).mean()), 1) if len(on) else np.nan},
            {"metric": "a_s1_median", "node": k,
             "value": round(float(on["a_s1"].median()), 2) if len(on) else np.nan},
            {"metric": "a_s2_median", "node": k,
             "value": round(float(on["a_s2"].median()), 2) if len(on) else np.nan},
            {"metric": "bracket_width_median", "node": k,
             "value": round(float((on["a_s1"] - on["a_s2"]).median()), 2) if len(on) else np.nan},
        ]
    diag = pd.DataFrame(rows)
    m = missingness.rename(columns={"pct_missing": "value"})[["node", "variable", "value"]]
    m["metric"] = "pct_missing_" + m["variable"]
    diag = pd.concat([diag, m[["metric", "node", "value"]]], ignore_index=True)
    for k, v in prov.items():
        diag[k] = v

    ID_LIKE = ("encounter_block", "hospitalization_id", "patient_id", "_dttm")
    leaked = [c for c in diag.columns if any(s in c for s in ID_LIKE)]
    assert not leaked, f"identifier-like columns in a shareable file: {leaked}"
    assert diag["value"].notna().any(), "diagnostics carry no values"

    path = share / f"{config['site_name']}_lmtp_df_diagnostics.csv"
    diag.to_csv(path, index=False)
    print(f"    wrote final_no_phi/{path.name}  {len(diag)} rows, PHI-checked")
    return prov


# %%
prov = stage_8_write(lmtp_df, exposure, missingness, config, design, REPO_ROOT)
print(f"\nlmtp_df: {len(lmtp_df):,} rows x {lmtp_df.shape[1]} columns")
