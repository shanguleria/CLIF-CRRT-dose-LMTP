# Output layout

Two trees, split by whether the contents can leave the site.

| Path | Contents | Shareable |
|---|---|---|
| `final_no_phi/` | Aggregate estimates, diagnostics, figures, the federated export set | **Yes** |
| `intermediate_phi/` | Patient-level node datasets, fitted model objects | **No. Never leaves the site** |

Both directories' contents are gitignored. These README files are tracked so the
layout is documented in a fresh clone.

Nothing person-level is exported. The federated contract is per-site point
estimates, the T x T influence-function covariance, n, learner coefficients, and
diagnostics.
