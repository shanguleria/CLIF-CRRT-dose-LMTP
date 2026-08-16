# config

| File | Purpose |
|---|---|
| `config_template.json` | Copy to `config.json` and edit for your site |
| `config.json` | Your site's settings. **Gitignored**, never leaves the site |
| `lmtp_design.json` | The estimand: delta ladder, floor, node schedule, shift form, estimator. Treat it as a protocol, not a tuning file |
| `clif_data_requirements.yaml` | Full CLIF table, column, and category specification. Pinned copy, see `code/vendor/VENDOR_SHA` |
| `outlier_config.json` | Physiological ranges for vitals, labs, vasopressors, CRRT flow rates. Pinned copy |

Changing a value in `lmtp_design.json` changes what is being estimated, and it is
the `definition_version` stamped onto every shareable output. It is not the place
to make a failing run pass.
