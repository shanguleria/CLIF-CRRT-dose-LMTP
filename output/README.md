# Output layout

Two trees, split by whether the contents can leave the site.

| Path | Contents | Shareable |
|---|---|---|
| `final_no_phi/` | Aggregate estimates, diagnostics, figures, the federated export set | **Yes** |
| `intermediate_phi/` | Patient-level node datasets, fitted model objects | **No. Never leaves the site** |

Both directories' contents are gitignored. This file, and the warning label in
`intermediate_phi/`, are tracked so the layout survives a fresh clone (git does
not track empty directories, hence the `.gitkeep` in `final_no_phi/`).

Every file written to `final_no_phi/` carries a provenance block (`site_id`,
`code_version`, `clif_version`, `definition_version`, `generated`) so a pooled
result can be traced back to the commit and estimand that produced it.
PHI-check the directory before sharing.

Nothing person-level is exported. The federated contract is per-site point
estimates, the T x T influence-function covariance, n, learner coefficients, and
diagnostics.
