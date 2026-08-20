# CRRT-dose-lmtp

Longitudinal Modified Treatment Policies (LMTP) for time-varying continuous renal
replacement therapy (CRRT) dose across the CLIF consortium.

**CLIF spec version:** 2.1.0 &nbsp;|&nbsp; **Definition version:** 0.2.0 (`config/lmtp_design.json`)

> ### Status: all three steps are built. One known blocker.
>
> `01_build_cohort.py`, `02_build_lmtp_df.py` and `03_lmtp_fit.R` all exist and
> run; the R environment is installed and locked (`renv.lock`, 74 packages,
> `lmtp` 1.5.4 on R 4.3.1), and `03`'s smoke stage passes.
>
> **Known blocker:** step 02 halts by design when a vasopressor cannot be unit
> converted. See [Known issues](#known-issues) before running. The halt is
> deliberate and dropping the rows is not a safe workaround.

---

## Quick start

```bash
# 1. Clone and install both stacks
git clone https://github.com/shanguleria/CLIF-CRRT-dose-LMTP.git
cd CLIF-CRRT-dose-LMTP
uv sync
Rscript -e 'renv::restore(prompt = FALSE)'

# 2. Create your site config
cp config/config_template.json config/config.json
#    set site_name, data_directory (your CLIF path), timezone, has_crrt_settings

# 3. Build the frame and run the smoke fit
./run_pipeline.sh                       # Windows: .\run_pipeline.ps1

# 4. Read the diagnostics, THEN estimate
Rscript code/03_lmtp_fit.R gate         # diagnostics only, no effect estimate
Rscript code/03_lmtp_fit.R expand       # the full delta ladder

# 5. Send output/final_no_phi/ to the coordinating center
#    (study Box folder; the link comes from the coordinating center)
```

**Five things to know before you start.**

1. **`has_crrt_settings` must be true.** Dose is the exposure and cannot be
   computed without CRRT flow rates, so both runners refuse to start without it.
2. **Step 3 stops on purpose.** The runner builds the frame and runs a cheap smoke
   fit, then hands you back the prompt. Step 4 is separate because stage 3 must not
   run until a human has read stage 2's diagnostics.
3. **Step 02 may halt on vasopressor units.** That is a deliberate guard, not a
   crash. See [Known issues](#known-issues) before working around it.
4. **Send `output/final_no_phi/` only.** `output/intermediate_phi/` is
   patient-level and never leaves your site. PHI-check before sending.
5. **Do not edit `config/lmtp_design.json`** to make a run work. It is the
   protocol, identical at every site, and changing it changes what is being
   estimated.

Fuller detail: [Configuration](#configuration), [Prerequisites](#prerequisites),
[Running the pipeline](#running-the-pipeline),
[Onboarding a new site](#onboarding-a-new-site).

---

## Objective

Estimate the effect of **reducing delivered CRRT dose** on 30-day in-hospital
mortality, treating dose as a **time-varying** exposure measured at 0, 24 and
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
FLOOR = 15 (primary delta = 5). The ladder **is** the dose-response result, so all
three rungs are reported. The shift is self-limiting and is **never clamped**:
clamping would break invertibility, and `lmtp` supplies no guard and will not warn.

**Why this framing.** A prior point-treatment analysis (high vs low dose at
30 mL/kg/hr) returned a null, and a time-varying marginal structural model was
abandoned for positivity violations: 29.5% of patients had a propensity score
above 0.95 at 12 hours, and truncation did not help. A feasibility audit across
ten sites found that under the shift framing the binned density ratio never
exceeds 6.0 at any site at delta = 2.5. The gain is structural, not a matter of
tuning.

Scope is **causal estimation only**. Descriptive epidemiology of this cohort lives
in the sibling repository `CLIF-epidemiology-of-CRRT` and is not reproduced here.

---

## Required CLIF tables and fields

The authoritative, complete specification is `config/clif_data_requirements.yaml`;
this is the summary. All tables are CLIF **2.1.0**.

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
rates cannot participate. Both runners refuse to start without it.

Effluent dose sums dialysate plus pre- and post-filter replacement for every
dose-eligible mode, counting whichever are charted, divided by **actual** body
weight. SCUF is excluded by definition rather than by name: it has neither
dialysate nor replacement, so it can never satisfy the "delivering" test.

The node statistic is a **time-weighted mean** over charted CRRT records in the
24-hour window. Charted zeros enter; uncharted gaps are not imputed.

### Time-varying covariates at each node

| Table | Required columns | Required categories |
|---|---|---|
| **clif_vitals** | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | map, sbp, dbp, heart_rate, respiratory_rate, spo2, temp_c, weight_kg, height_cm |
| **clif_labs** | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value_numeric` | creatinine, bun, potassium, bicarbonate, lactate, platelet_count, bilirubin_total, po2_arterial, fio2_set, sodium, chloride, albumin |
| **clif_medication_admin_continuous** | `hospitalization_id`, `admin_dttm`, `med_category`, `med_dose`, `med_dose_unit` | norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine, dobutamine, angiotensin (for norepinephrine-equivalent) |
| **clif_respiratory_support** | `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `fio2_set`, `peep_set` | IMV, NIPPV, High Flow NC, Nasal Cannula, Room Air |
| **clif_microbiology_culture** | `hospitalization_id`, `collected_dttm`, `result_category` | *(optional)* sepsis flag; absent means the flag is NA and the pipeline continues |

These feed P/F (with S/F fallback and a source indicator), norepinephrine
equivalent, inotrope status, lactate, potassium, pH, bicarbonate, BUN and IMV
status, computed **at each exposure node**, so the `L_t -> A_t` ordering holds.
That ordering is asserted at import against the node schedule, not left to a
comment.

**A stated limitation:** there is **no volume assessment**. Urine output and fluid
balance were dropped rather than proxied, because no `intake_output` table exists
in this CLIF version and `ultrafiltration_out` is unstable.

---

## Cohort identification

- **Population:** hospitalized adults receiving CRRT for acute kidney injury.
- **Unit of analysis:** the **encounter block**, not the hospitalization.
  Encounters are stitched first, then judged; blocks are the unit throughout.
- **Time zero:** CRRT initiation, the first *delivering* `clif_crrt_therapy`
  record per encounter block. The outcome window is anchored there.
- **Exclusions:** pre-existing ESRD, by ICD code at block grain.
- **Not excluded:** short CRRT courses. The sibling analysis dropped courses under
  24 hours because a point-treatment design needs 24 hours to define the exposure.
  LMTP handles early death structurally, so that exclusion is dropped here and
  would now induce selection.
- **Outcome:** 30-day in-hospital mortality, with **discharge alive** as a
  competing event. In-hospital by construction: death is read from discharge
  disposition, so post-discharge deaths are unobservable, and **hospice discharge
  counts as a death**.
- **Time grid:** one six-period axis, days `[1, 2, 3, 7, 14, 30]`. The policy
  intervenes in periods 1-3 only (0/24/48h); periods 4-6 carry the outcome curve
  to day 30 unintervened, at a density ratio of exactly 1.

**The cohort is written here, not vendored.** An earlier plan copied
`00_cohort.py` and four dependencies from `CLIF-epidemiology-of-CRRT`; that was
abandoned on 2026-08-16 because the file is a converted notebook with no callable
API, and because `clifpy` 0.4.9 already supplies `stitch_encounters`,
`compute_sofa_polars`, `create_wide_dataset` and `apply_outlier_handling`.
`code/vendor/` now pins **two config files and no code**.

The consequence is worth stating plainly: because the cohort logic is
reimplemented rather than copied, **matching the sibling analysis's N is a claim to
be tested at each site, not a fact**. `config/lmtp_design.json` records the
reconciliation under `cohort._must_reconcile`.

---

## Configuration

Two files, and putting a value in the wrong one is a scientific error rather than
a tidiness problem.

| | File | Contains | Differs by site? |
|---|---|---|---|
| **Site** | `config/config.json` | `data_directory`, `filetype`, `timezone`, `site_name`, `has_crrt_settings`, `n_workers` | **Yes** |
| **Protocol** | `config/lmtp_design.json` | delta ladder, floor, node schedule, study window, estimator, competing event, cohort rules | **No.** Identical everywhere |

The test: *if two sites set this differently, is the pooled result still
meaningful?* If no, it is protocol.

Copy the template and edit it for your site:

```bash
cp config/config_template.json config/config.json
```

```json
{
    "site_name": "Your_Site_Name",
    "clif_version": "2.1.0",
    "data_directory": "/path/to/clif/tables/",
    "filetype": "parquet",
    "timezone": "America/Chicago",
    "project_root": "/path/to/CRRT-dose-lmtp",
    "output_dir": "output",
    "has_crrt_settings": true,
    "n_workers": 6
}
```

`data_directory`, `filetype` and `timezone` are named to match what `clifpy`
expects, so one file serves both the library and this project's own code with no
translation layer. Note this differs from the sibling `CLIF-epidemiology-of-CRRT`,
which calls the first two `tables_path` and `file_type`, so a config **cannot** be
copied between the two repos unchanged.

`n_workers` is a **site** setting, not protocol: `lmtp` wraps each cross-fitting
fold in `future(..., seed = TRUE)`, so results are identical under any worker
count and only wall-clock changes. Set it to 1 to force sequential.

`config/config.json` is gitignored and never leaves your site.

`config/lmtp_design.json` holds the estimand itself. **Do not edit it to make a run
work.** Changing a value there changes what is being estimated, and it is the
`definition_version` stamped onto every shareable output.

### Where outputs land

| Path | Contents | Shareable |
|---|---|---|
| `output/final_no_phi/` | Aggregate estimates, diagnostics, the federated export set | **Yes.** This is what the coordinating center receives |
| `output/intermediate_phi/` | Patient-level node datasets, fitted objects | **No. Never leaves the site** |
| `output/logs/` | Run and clifpy validation logs | Review before sharing |

Nothing person-level is ever exported. The federated contract is per-site point
estimates, the T x T influence-function covariance matrix, n, learner
coefficients, and diagnostics.

---

## Prerequisites

- **Python 3.11** (3.11.15 pinned via `.python-version`)
- **[uv](https://docs.astral.sh/uv/)** package manager
- **R 4.3+** (developed and locked against 4.3.1) with `renv`
- Read access to your site's CLIF 2.1.0 tables, named `clif_<table>.parquet`

```bash
uv sync                      # creates .venv and installs the pinned Python stack
Rscript -e 'renv::restore()' # installs the pinned R stack from renv.lock
```

The Python compute stack is pinned to **exact** versions for cross-site numeric
reproducibility, notably [`clifpy`](https://pypi.org/project/clifpy/)`==0.4.9`:
its minor releases have changed CLIF datetime timezone handling, which silently
moves every windowed exposure node. Presentation libraries float. For sites
without uv, `requirements.txt` carries the same pins.

`renv.lock` pins **74** R packages including `lmtp` 1.5.4. This repo builds its own
lock against R 4.3.1 and deliberately does **not** reuse the `fluid_ARDS` lock,
which targets 4.5.2 and fails to restore here. CRAN's macOS binaries for R 4.3 are
frozen at `lmtp` 1.5.2, so 1.5.4 installs from source; this works because `lmtp` is
pure R. Detail and the version landmine: `code/R_PACKAGES.md`.

---

## Running the pipeline

The runner builds the analysis frame and runs the cheap smoke fit, then **stops**.

Stopping is deliberate. `03_lmtp_fit.R` runs in gated stages, and stage 3 must not
run until a human has read stage 2's diagnostics: once an effect estimate has been
seen, every later decision about trimming or covariates is contaminated. A runner
that drove all three stages end to end would defeat the gate.

### macOS / Linux

```bash
./run_pipeline.sh                      # 01 -> 02 -> 03 smoke, then stops
Rscript code/03_lmtp_fit.R gate        # diagnostics ONLY, no effect estimate
# ...read the diagnostics, then:
Rscript code/03_lmtp_fit.R expand      # the full delta ladder
```

### Windows

```powershell
.\run_pipeline.ps1
Rscript code\03_lmtp_fit.R gate
Rscript code\03_lmtp_fit.R expand
# if execution policy blocks it:
#   powershell -ExecutionPolicy Bypass -File .\run_pipeline.ps1
```

**Exit codes:** `1` preflight or step failure, `2` a pipeline step is missing,
`3` frame built but `Rscript` was not on PATH. A runner that exits 0 without doing
anything is how a site ends up believing it has results.

**Runtime.** The gate stage is roughly 2 hours and the full expand grid roughly
6 hours at `folds = 10`, hardware dependent. Useful arithmetic for judging whether
a long run is pathological or expected: one SuperLearner call on this design matrix
is ~40 s at V=10, and one `lmtp` fit is `tau x folds x 2 = 120` such calls.

---

## Pipeline steps

| Step | Language | Script | Description |
|---|---|---|---|
| 01 | Python | `code/01_build_cohort.py` | Cohort identification, ESRD exclusion, CRRT initiation, dose series, outcomes |
| 02 | Python | `code/02_build_lmtp_df.py` | Exposure and covariate nodes at 0/24/48h; `L_t -> A_t` ordering asserted at import |
| 03 | R | `code/03_lmtp_fit.R` | `lmtp_sdr` fit over the delta ladder, positivity diagnostics, influence-function exports |

**There is no step 00.** See [Cohort identification](#cohort-identification) for why
the vendoring plan was abandoned.

Step 03 takes a stage argument: `smoke` (one cheap fit, SL.glm, folds = 2),
`gate` (natural course plus the primary policy, full learner library, diagnostics
and **no** effect estimate), `expand` (the full delta ladder x {S1, S2} x
{SDR, TMLE}).

### What the steps write

| Step | Patient-level (`intermediate_phi/`) | Shareable (`final_no_phi/`) |
|---|---|---|
| 01 | `cohort.parquet`, `dose_series.parquet`, `block_map.parquet` | `<SITE>_strobe_counts.csv` |
| 02 | `lmtp_df.parquet` (one row per encounter block) | `<SITE>_lmtp_df_diagnostics.csv` |
| 03 | Fitted objects, per-observation influence functions | Estimates, diagnostics, the federated export set |

Steps 02 and 03 read every measurement rule from `config/lmtp_design.json`. They
decide nothing themselves: a change to a lookback window or a summary rule is a
protocol amendment made in that file, which bumps `definition_version`. Both steps
carry a coverage assertion that fails loudly if a config key is declared but never
consumed.

---

## Known issues

**Step 02 halts on vasopressor unit conversion.** `clifpy` 0.4.9 leaves the **raw**
value in `med_dose_converted` when it cannot convert a unit, and reports the failure
only in `med_dose_unit_converted` — so a check for NA sees nothing. An angiotensin
dose charted in ng/kg/min is then read as mcg/kg/min, and the x10 norepinephrine
equivalent coefficient turns it into a 200-fold error.

Step 02 therefore **raises** rather than continue. Dropping the unconvertible rows
is not a safe workaround: it removes a drug from the norepinephrine equivalent for
exactly the patients who received it, which is differential misclassification of a
confounder in a known direction. Across sites it would also let each site compute a
different NEE.

If your site charts every vasopressor in a unit `clifpy` converts, you will not see
this. Report it if you do; the fix belongs upstream.

---

## Project structure

```
CRRT-dose-lmtp/
├── CLAUDE.md                     Project guidance and inherited decisions
├── README.md
├── pyproject.toml                Pinned Python stack
├── requirements.txt              Same pins, for sites without uv
├── uv.lock                       Resolved Python lockfile
├── .python-version               3.11.15
├── .Rprofile                     Activates renv
├── renv.lock                     74 R packages, lmtp 1.5.4 on R 4.3.1
├── renv/                         R environment
├── run_pipeline.sh / .ps1        Runners (macOS/Linux and Windows)
│
├── code/
│   ├── 01_build_cohort.py        Step 01
│   ├── 02_build_lmtp_df.py       Step 02
│   ├── 03_lmtp_fit.R             Step 03, gated stages
│   ├── R_PACKAGES.md             R stack and the version landmine
│   ├── README.md
│   └── vendor/
│       ├── VENDOR_SHA            Upstream pin and manifest (two config files)
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
│   ├── intermediate_phi/         Patient-level, never shared
│   └── logs/
│
├── references/                   Papers (gitignored; README tracked)
└── tests/
    └── test_vendor_integrity.py  Byte-equality guard on the pinned config files
```

Two directories exist in the coordinating site's working copy but are **gitignored
and absent from a clone**: `.claude/` (internal planning and the feasibility audit,
which names individual sites) and `docs/` (build walkthroughs written around one
site's worked numbers). Nothing in the pipeline depends on either.

---

## Definitions & provenance

If someone asks *"which code and which definitions produced this number?"*, the
answer is auditable:

- **`config/lmtp_design.json` carries a `definition_version`** (currently `0.2.0`)
  and is the single source of truth for the estimand. It is machine-readable so it
  is diffable and hashable, and any change to it is a protocol amendment that
  leaves earlier results identifiable.
- **Every shareable output carries a provenance block**: `site_id`,
  `code_version` (git SHA), `clif_version`, `definition_version`, `generated`.
- **Config coverage is asserted, not assumed.** Steps 02 and 03 fail loudly if a
  key in `lmtp_design.json` is declared but never consumed, which is the failure
  mode where a config reads as policy while changing nothing.
- **Estimator:** `lmtp_sdr` primary, `lmtp_tmle` secondary. `lmtp_ipw()` and
  `lmtp_sub()` are defunct in `lmtp` 1.5.x and raise errors, so the estimation
  diagnostics are SDR-vs-TMLE agreement plus the density-ratio distribution and
  trimmed fraction.

---

## Onboarding a new site

The commands are in [Quick start](#quick-start); they are not repeated here, so
the two cannot drift apart. This section is what to check once they have run.

**Before trusting any number, check three things.**

1. **`has_crrt_settings` must be true.** Without flow rates there is no dose, and
   dose is the exposure. Both runners refuse to start otherwise.
2. **Vasopressor units.** See [Known issues](#known-issues). If step 02 raises on
   unit conversion, that is the guard working, not a bug in your data.
3. **`crrt_mode_category` is not trustworthy on its own.** At least one site
   labels every sustained ultrafiltration-only course as `cvvhd`. This pipeline
   does not branch on the mode label for dose eligibility, precisely for that
   reason, but do not build a modality analysis on that column without checking.

**Hard-failure gotchas.** `clifpy` needs files named `clif_<table>.parquet`. A
different `clif_version` can shift column names. Flow-rate outlier bounds are
applied once at load, so a site with wildly out-of-range charted flows will see
records nulled rather than silently used.

---

## Data safety

Raw CLIF tables contain protected patient data. They are never committed, never
copied into this tree, and never read into an analysis transcript. Scripts print
**aggregates only**, never rows.

- `output/intermediate_phi/` is the PHI working space and is gitignored. It holds
  one row per patient and never leaves the site.
- `output/final_no_phi/` is aggregate by construction and is what the coordinating
  center receives. PHI-check it before sending.
- `config/config.json` is gitignored, since it carries your data path.

**Site anonymization.** Anything audience-facing (manuscripts, abstracts, posters,
slides, public dashboards) refers to the *number* of participating sites and uses
anonymized "Site 1, Site 2, …" labels for per-site figures. Author affiliations are
the only place real institution names appear.

---

## Related repositories

| Repository | Relationship |
|---|---|
| `CLIF-epidemiology-of-CRRT` | Descriptive epidemiology and the point-treatment analysis of the same cohort. Source of the two pinned config files |
| `fluid_ARDS` | Methodological source of truth for LMTP. Inherited by reference, never copied |
| `crrt-manuscript-tools` | Private coordinating-center tooling for the CJASN manuscript |

---

## Acknowledgements

Built on the [Common Longitudinal ICU Format (CLIF)](https://clif-consortium.github.io/website/)
and the [`clifpy`](https://pypi.org/project/clifpy/) library, with estimation by
[`lmtp`](https://cran.r-project.org/package=lmtp) (Williams & Díaz).

Method references: Díaz, Williams, Hoffman & Schenck, *Journal of the American
Statistical Association* 2023 (modified treatment policies); Díaz, Hoffman &
Hejazi, *Lifetime Data Analysis* 2024;30:213-236 (competing risks).
