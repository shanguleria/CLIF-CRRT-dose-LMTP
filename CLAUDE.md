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

## Current status: scaffold only

Nothing is estimable yet. `02_build_lmtp_df.py` and `03_lmtp_fit.R` do not
exist, no code is vendored, and no R packages are installed. What exists is the
repo, the pinned Python environment, the design spec, and the rescued
feasibility evidence.

**The build order is in `.claude/lmtp_feasibility_findings.md` section 7.** Step 1
of it is a blocking decision, described below.

## The blocking decision, before any dataset code

**How to handle CRRT discontinuation.** The exposure is undefined once CRRT
stops, and 39% of encounter blocks at the coordinating site have no third node.
Three options, and they are **three different estimands, not three
implementations of one**:

1. treat liberation as a competing event;
2. carry dose forward as 0, which silently changes the question to "dose or no dose";
3. define the policy only while the patient is on therapy.

This determines the node schema, so it comes before `02_build_lmtp_df.py`. It is
recorded as `discontinuation_handling: null` in `config/lmtp_design.json`, and
that null is deliberate. Do not pick one implicitly while writing dataset code.

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
| Node statistic | **windowed median** | The archived MSM used means, which reintroduce the charting-lag artefact |
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

## Vendoring contract

`00_cohort.py` and its four dependencies are **copied verbatim** from
`CLIF-epidemiology-of-CRRT` at SHA `ee4774b`, never rewritten. That is where the
science lives: the ESRD gate, CRRT initiation per encounter block, the
modality-agnostic effluent formula, the first-3h dose median, the 30-day anchor
to CRRT initiation. Reimplementing it forks the cohort definition, and then two
papers from the same consortium report two different Ns.

- Manifest and pin: `code/vendor/VENDOR_SHA`
- Sync: `bash code/vendor/sync_vendor.sh`
- Guard: `uv run pytest tests/test_vendor_integrity.py` asserts byte-equality

**Never edit a vendored file in place.** To take an upstream change, move the SHA
and re-sync in one deliberate commit.

Do **not** bring `02_construct_crrt_tableone.py`, `03_crrt_epidemiology.py`,
`03b`, `02c`, `06`, `04_build_causal_df.py`, `05`, `05b`, or the sibling's renv.

**Do not inherit `04`'s schema.** It reads `tableone_analysis_df.parquet`, so its
baseline block is coupled to the Table 1 script. More importantly, **`04`
violates the `L_t -> A_t` ordering** by pairing a 0-24h mean dose with covariates
measured at 24h. Do not reproduce that.

## Repo conventions

Read `.gitignore`'s header before adding a rule to it. This repo is **private**,
so unlike the sibling consortium repo, **`.claude/`, `CLAUDE.md`, and `docs/` are
tracked**. `.claude/lmtp_feasibility_findings.md` is the entire evidentiary basis
for the project and previously existed as one untracked copy under a gitignored
path; ignoring `.claude/` here would recreate that.

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
