# CRRT-dose-lmtp

Longitudinal Modified Treatment Policies (LMTP) for time-varying continuous renal
replacement therapy (CRRT) dose across the CLIF consortium.

**CLIF Version:** 2.1.0

> ### Status: the Python half runs. The fit does not exist yet.
>
> `01_build_cohort.py` and `02_build_lmtp_df.py` are complete and run end to end,
> from raw CLIF 2.1.0 tables to the wide analysis frame `lmtp` consumes.
> `03_lmtp_fit.R` **does not exist yet** and no R packages are installed, so the
> estimation step will not run. Sites can build the analysis frame and the
> shareable diagnostics; they cannot yet fit anything.

---

## Objective

Estimate the effect of **reducing delivered CRRT dose** on 30-day in-hospital
mortality, treating dose as a **time-varying** exposure measured at 0, 24, and
48 hours after CRRT initiation, with time-varying confounders and discharge alive
as a competing event.

The intervention is a *modified treatment policy*, not a static dose assignment:
each patient's own observed dose is shifted down by a fixed amount, and only if
the result stays above a clinical floor.

```
d(a) = a - delta   if a - delta >= FLOOR
       a           otherwise
```

The policy is evaluated over a ladder of delta = {2.5, 5, 10} mL/kg/hr with
FLOOR = 15. The ladder is the dose-response result, so all three rungs are
reported.

**Why this framing.** A prior point-treatment analysis (high vs low dose at
30 mL/kg/hr) returned a null, and a time-varying marginal structural model was
abandoned for positivity violations: 29.5% of patients had a propensity score
above 0.95 at 12 hours, and truncation did not help. A feasibility audit across
ten sites found that under the shift framing the binned density ratio never
exceeds 6.0 anywhere at delta = 2.5. The gain is structural, not a matter of
tuning. Full audit: `.claude/lmtp_feasibility_findings.md`.

Scope is **causal estimation only**. Descriptive epidemiology of this cohort
lives in the sibling repository `CLIF-epidemiology-of-CRRT` and is not
reproduced here.

---

## Required CLIF tables and fields

The cohort is built by code vendored from `CLIF-epidemiology-of-CRRT`, so the
table requirements are identical to that project. The authoritative, complete
specification is `config/clif_data_requirements.yaml`; this is the summary.

### Shared by every step

| Table | Required columns | Required categories |
|---|---|---|
| **clif_patient** | `patient_id`, `race_category`, `ethnicity_category`, `sex_category` | - |
| **clif_hospitalization** | `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `discharge_category` | `discharge_category` drives the mortality outcome |
| **clif_adt** | `hospitalization_id`, `in_dttm`, `out_dttm`, `location_category` | ICU, Ward, ED, OR, Procedural, Other |
| **clif_hospital_diagnosis** | `hospitalization_id`, `diagnosis_code`, `diagnosis_code_format` | ICD-9/10, used for the ESRD exclusion and Charlson |

### Exposure

| Table | Required columns | Required categories |
|---|---|---|
| **clif_crrt_therapy** | `hospitalization_id`, `recorded_dttm`, `crrt_mode_category`, `blood_flow_rate`, `dialysate_flow_rate`, `pre_filter_replacement_fluid_rate`, `post_filter_replacement_fluid_rate`, `ultrafiltration_out` | CVVH, CVVHD, CVVHDF, SCUF |

`has_crrt_settings: true` is **mandatory**. Dose cannot be computed from the
minimal two-column CRRT table, and dose is the exposure, so a site without flow
rates cannot participate.

Effluent dose sums dialysate plus pre- and post-filter replacement for every
dose-eligible mode, counting whichever are charted, divided by body weight.
SCUF is excluded.

### Time-varying covariates at each node

| Table | Required columns | Required categories |
|---|---|---|
| **clif_vitals** | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | map, sbp, dbp, heart_rate, respiratory_rate, spo2, temp_c, weight_kg, height_cm |
| **clif_labs** | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value_numeric` | creatinine, bun, potassium, bicarbonate, lactate, platelet_count, bilirubin_total, po2_arterial, fio2_set, sodium, chloride, albumin |
| **clif_medication_admin_continuous** | `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit` | norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine, dobutamine, angiotensin (for norepinephrine-equivalent) |
| **clif_respiratory_support** | `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `fio2_set`, `peep_set` | IMV, NIPPV, High Flow NC, Nasal Cannula, Room Air |
| **clif_microbiology_culture** | `hospitalization_id`, `collected_dttm`, `result_category` | *(optional)* sepsis flag; absent means the flag is NA and the pipeline continues |

These feed SOFA-1 (Vincent 1996), P/F and S/F ratios, norepinephrine-equivalent,
and IMV status computed **at each exposure node**, so that the `L_t -> A_t`
ordering holds.

---

## Cohort identification

Identical by design to the sibling CJASN analysis. The cohort code is vendored at
a pinned commit rather than reimplemented, precisely so that two papers from the
same consortium cannot report two different Ns.

- **Population:** hospitalized adults receiving CRRT for acute kidney injury.
- **Unit of analysis:** the **encounter block**, not the hospitalization. CRRT
  initiation is the first `clif_crrt_therapy` record per encounter block.
- **Time zero:** CRRT initiation. The 30-day outcome window is anchored there.
- **Exclusions:** pre-existing ESRD, identified by ICD code.
- **Not excluded:** short CRRT courses. The sibling analysis dropped courses
  under 24 hours because a point-treatment design needs 24 hours to define the
  exposure. LMTP handles early death structurally, so that exclusion is dropped
  here and would now induce selection.
- **Outcome:** 30-day in-hospital mortality, with discharge alive as a competing
  event. In-hospital by construction: death is read from discharge disposition,
  so post-discharge deaths are unobservable, and hospice discharge counts as a
  death.

---

## Configuration

Copy the template and edit it for your site:

```bash
cp config/config_template.json config/config.json
```

```json
{
    "site_name": "Your_Site_Name",
    "data_directory": "/path/to/clif/tables/",
    "filetype": "parquet",
    "timezone": "America/Chicago",
    "project_root": "/path/to/CRRT-dose-lmtp",
    "output_dir": "output",
    "has_crrt_settings": true
}
```

`data_directory`, `filetype`, and `timezone` are named to match what `clifpy`
expects, so one file serves both the library and this project's own code with no
translation layer. The remaining keys are ours; clifpy ignores what it does not
recognise. Note this differs from the sibling `CLIF-epidemiology-of-CRRT`, which
calls the first two `tables_path` and `file_type`, so a config cannot be copied
between the two repos unchanged.

`config/config.json` is gitignored and never leaves your site.

`config/lmtp_design.json` holds the estimand itself: the delta ladder, the floor,
the node schedule, the shift form, the estimator, and the competing-event
handling. **Do not edit it to make a run work.** Changing a value there changes
what is being estimated, and it is the `definition_version` stamped onto every
shareable output.

### Where outputs land

| Path | Contents | Shareable |
|---|---|---|
| `output/final_no_phi/` | Aggregate estimates, diagnostics, figures, the federated export set | **Yes.** This is what the coordinating center receives |
| `output/intermediate_phi/` | Patient-level node datasets, fitted objects | **No. Never leaves the site** |

Nothing person-level is ever exported. The federated contract is per-site
point estimates, the influence-function covariance matrix, n, learner
coefficients, and diagnostics.

---

## Prerequisites

- **Python 3.11** (3.11.15 pinned via `.python-version`)
- **[uv](https://docs.astral.sh/uv/)** package manager
- **R 4.3+** with `renv`
- Read access to your site's CLIF 2.1.0 tables

```bash
uv sync                 # creates .venv and installs the pinned stack
```

The Python compute stack is pinned to exact versions for cross-site numeric
reproducibility. `clifpy` especially: its minor releases have changed CLIF
datetime timezone handling, which silently moves every windowed exposure node.

For sites without uv, `requirements.txt` carries the same pins.

R packages are **not yet installed**. See `code/R_PACKAGES.md`, which records the
required set and an R-version conflict to be aware of before installing.

---

## Running the pipeline

*(not built)* Neither runner script executes anything yet; both currently report
which steps are missing and exit.

### macOS / Linux

```bash
bash run_pipeline.sh
```

### Windows

```powershell
.\run_pipeline.ps1
```

---

## Pipeline steps

| Step | Language | Script | Status |
|---|---|---|---|
| 01 | Python | `code/01_build_cohort.py` | **Runs.** Cohort identification, ESRD exclusion, CRRT initiation, dose series, outcomes |
| 02 | Python | `code/02_build_lmtp_df.py` | **Runs.** Exposure and covariate nodes at 0/24/48h, `L_t -> A_t` ordering asserted |
| 03 | R | `code/03_lmtp_fit.R` | *(not built)* `lmtp_sdr` fit over the delta ladder, influence-function exports |

There is no step 00. An earlier plan vendored `00_cohort.py` and four dependencies
from `CLIF-epidemiology-of-CRRT`; that was abandoned on 2026-08-16 because the file
is a converted notebook with no callable API and because `clifpy` 0.4.9 already
supplies `stitch_encounters`, `compute_sofa_polars`, `create_wide_dataset` and
`apply_outlier_handling`. `code/vendor/` now pins two **config** files only.

### What the steps write

| Step | Patient-level (`intermediate_phi/`) | Shareable (`final_no_phi/`) |
|---|---|---|
| 01 | `cohort.parquet`, `dose_series.parquet`, `block_map.parquet` | `<SITE>_strobe_counts.csv` |
| 02 | `lmtp_df.parquet` (one row per encounter block) | `<SITE>_lmtp_df_diagnostics.csv` |

Step 02 reads its every measurement rule from `config/lmtp_design.json` under
`covariates`. It decides nothing itself: a change to a lookback window or a summary
rule is a protocol amendment made in that file, which bumps `definition_version`.

---

## Project structure

```
CRRT-dose-lmtp/
├── CLAUDE.md                     Project guidance, inherited decisions
├── README.md
├── pyproject.toml                Pinned Python stack
├── requirements.txt              Same pins, for sites without uv
├── .python-version               3.11.15
├── .Rprofile                     Activates renv
├── renv/                         R environment (scaffolded, empty)
├── run_pipeline.sh / .ps1        Runners (not yet functional)
│
├── .claude/
│   ├── lmtp_feasibility_findings.md   The GO decision and its evidence
│   ├── feasibility_results/           Aggregate audit CSVs, 10 sites
│   ├── feasibility_code/              The three audit scripts
│   ├── claude-todo.md                 Prioritized next steps
│   └── claude-progress.md             Session log and decisions
│
├── code/
│   ├── R_PACKAGES.md             R stack plan and version landmine
│   └── vendor/
│       ├── VENDOR_SHA            Upstream pin and manifest
│       └── sync_vendor.sh        Pull manifest files at the pin
│
├── config/
│   ├── config_template.json      Copy to config.json
│   ├── lmtp_design.json          The estimand. Treat as a protocol
│   ├── clif_data_requirements.yaml
│   └── outlier_config.json
│
├── output/
│   ├── final_no_phi/             Shareable aggregates
│   └── intermediate_phi/         Patient-level, never shared
└── tests/
    └── test_vendor_integrity.py  Byte-equality guard on vendored code
```

---

## Data safety

Raw CLIF tables contain protected patient data. They are never committed, never
copied into this tree, and never read into an analysis transcript. Only
aggregate outputs under `output/final_no_phi/` are shareable, and only
`output/final_no_phi/` is assembled for transfer to the coordinating center.

---

## Related repositories

| Repository | Relationship |
|---|---|
| `CLIF-epidemiology-of-CRRT` | Source of the vendored cohort code. Descriptive epidemiology and the point-treatment analysis |
| `crrt-manuscript-tools` | Private coordinating-center tooling for the CJASN manuscript |
| `fluid_ARDS` | Methodological source of truth for LMTP. Inherited by reference, never copied |
