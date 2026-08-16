# R environment for the LMTP fit

**Status: scaffolded, not populated.** `renv::init(bare = TRUE)` has created
`.Rprofile`, `renv/activate.R`, and `renv/settings.json`. No packages are
installed and there is no `renv.lock` yet. `Rscript` will print
"The project is out-of-sync" until the first `renv::snapshot()`; that is
expected at this stage, not a fault.

## The R version landmine, read before installing anything

| | Version |
|---|---|
| System R on this machine | **4.3.1** (2023-06-16) |
| Sibling `CLIF-epidemiology-of-CRRT` `renv.lock` targets | **4.5.2** |

`renv::restore()` against the sibling's lock **fails on this machine**: there are
no matching binaries for 4.3.1, so it compiles from source and dies on a batch
(`arrow`, `Matrix`, `RcppArmadillo`, `V8`, `systemfonts`, `fracdiff`), emptying
`renv/library` on the way out.

**Consequence for this repo: build our own lock against 4.3.1. Do not copy,
symlink, or restore the sibling's `renv.lock`.** The two projects have different
R stacks and only one of them needs `lmtp`.

If R is later standardised to 4.5.x across sites, regenerate this lock against
that version deliberately, in its own commit.

## Packages needed, and what is already available

System library: `/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library` (300 packages).

| Package | Purpose | System lib | Sibling lock |
|---|---|---|---|
| `lmtp` | the estimator; target v1.5.4 | **no** | no |
| `mlr3superlearner` | learner stack `lmtp` expects | **no** | no |
| `earth` | MARS learner | **no** | no |
| `nnls` | convex learner weights | **no** | no |
| `SuperLearner` | legacy learner interface | yes | yes (2.0-40) |
| `ranger` | random forest learner | yes | no |
| `glmnet` | penalised regression learner | yes | yes (5.0) |
| `xgboost` | boosted-tree learner | check | yes |
| `survival` | outcome handling | yes | yes |
| `data.table` | wide-format assembly | yes | yes |
| `progressr`, `future` | progress + parallelism for `lmtp` | check | check |

So the four that must be installed fresh are **`lmtp`, `mlr3superlearner`,
`earth`, `nnls`**.

## Install, when the time comes

```r
# From the repo root, so .Rprofile activates renv.
renv::install(c("lmtp", "mlr3superlearner", "earth", "nnls",
                "SuperLearner", "ranger", "glmnet", "xgboost",
                "survival", "data.table", "progressr", "future"))
renv::snapshot()
```

Then commit `renv.lock` together with whatever script triggered the need.

**Verify `lmtp` actually builds on R 4.3.1 before designing around it.** This is
the first real risk to retire on the R side. If it needs a newer R, that changes
the deployment story for every participating site, and it is much cheaper to
learn now than after the dataset is built.

## Open mechanical question to settle at install time

`lmtp` v1.5.4 behaviour when `tau_trt < tau_outcome`, i.e. exposure nodes at
0/24/48h but outcome nodes out to day 30. The outcome grid is coarsened to
day 3/7/14/30 (`config/lmtp_design.json`) precisely so the T x T
influence-function covariance stays conditionable at small sites, but the
package mechanics for the mismatched horizon need verifying against the manual
rather than assumed. See `fluid_ARDS/plans/lmtp_discussion.md` section 6 and
open question 5 in `.claude/lmtp_feasibility_findings.md`.
