# References

Papers this project's design decisions rest on. **The PDFs are gitignored**: they
are publisher copyright and this repository is public, so they are read locally and
never redistributed. Obtain them from your institution's library.

| Citation | DOI | What it decides here |
|---|---|---|
| Quickfall D, La A, Bell E, Pisano J, Costello P, Gunning S, Koyner JL. Variability in Antibiotic Dosing and Resistance Development during Continuous Renal Replacement Therapy in Critically Ill Patients. *Blood Purif* 2026;55:197-209. | `10.1159/000550381` | The node statistic. Delivered CRRT dose is computed from hourly effluent flow rates averaged within patient; short interruptions such as filter changes are retained as real delivered dose; gaps of 48h or more are excluded rather than zero-filled. See `exposure._node_statistic_why` in `config/lmtp_design.json`. Also a same-institution external check on our dose distribution (todo 1b). |
| Diaz I, Williams N, Hoffman KL, Schenck EJ. Non-Parametric Causal Effects Based on Longitudinal Modified Treatment Policies. *JASA* 2023;118(542):846-857. | `10.1080/01621459.2021.1955691` | The estimand and the estimator. Vector-valued exposure with censoring intervened to 1 (p.850), Assumption 1 positivity, and the SDR estimator. |
| Williams N, Diaz I. `lmtp`: Non-Parametric Causal Effects of Feasible Interventions Based on Modified Treatment Policies. R package v1.5.4 (manual). | | Implementation mechanics: `cens` is mandatory for survival outcomes, `compete` is separate and applies only when `outcome_type = "survival"`, `mtp = TRUE`. |

## Not yet obtained

Both bear on the competing-event framing and should be read before it is finalised
in the manuscript (todo 1):

- Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernan MA. A causal framework for
  classical statistical estimands in failure-time settings with competing events.
  *Stat Med* 2020;39(8):1199-1236.
- Diaz I, Hoffman KL, Hejazi NS. Causal survival analysis under competing risks using
  longitudinal modified treatment policies. *Lifetime Data Anal* 2024;30:213-236.
