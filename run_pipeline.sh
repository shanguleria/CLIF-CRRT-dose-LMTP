#!/usr/bin/env bash
# CRRT-dose-lmtp pipeline runner (macOS / Linux).
#
# Builds the analysis frame and runs the cheap smoke fit, then STOPS.
#
# Stopping is deliberate. 03_lmtp_fit.R runs in gated stages and stage 3 must not
# run until a human has read stage 2's diagnostics: once an effect estimate has
# been seen, every later decision about trimming or covariates is contaminated.
# A runner that drove all three stages end to end would defeat the gate, so the
# gate and expand stages are invoked by hand. This script prints the commands.
#
# Exits non-zero when it cannot complete: a runner that exits 0 without doing
# anything is how a site ends up believing it has results.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONIOENCODING=utf-8
cd "$SCRIPT_DIR"

CFG="${CLIF_CONFIG:-$SCRIPT_DIR/config/config.json}"

echo "=============================================================="
echo " CRRT-dose-lmtp"
echo "=============================================================="

# --- Preflight ---------------------------------------------------------------
if [ ! -f "$CFG" ]; then
    echo "ERROR: config not found at $CFG" >&2
    echo "       Fix: cp config/config_template.json config/config.json" >&2
    echo "       then edit it for your site." >&2
    exit 1
fi

if ! command -v uv >/dev/null 2>&1; then
    echo "ERROR: uv is not installed. See https://docs.astral.sh/uv/" >&2
    exit 1
fi

SITE=$(uv run --quiet python -c "import json,sys;print(json.load(open(sys.argv[1]))['site_name'])" "$CFG")
HAS_SETTINGS=$(uv run --quiet python -c "import json,sys;print(json.load(open(sys.argv[1])).get('has_crrt_settings',False))" "$CFG")

echo "Site   : $SITE"
echo "Config : $CFG"
echo

if [ "$HAS_SETTINGS" != "True" ]; then
    echo "ERROR: has_crrt_settings is not true." >&2
    echo "       Dose is the exposure for this study and cannot be computed" >&2
    echo "       without CRRT flow rates. A site without them cannot participate." >&2
    exit 1
fi

# Lockfile completeness. Cheap, needs no R, and fails BEFORE step 01's long table
# read rather than at step 03's library() call an hour later. nanoparquet shipped
# unrecorded in renv.lock and only the coordinating machine survived it, because the
# package happened to be installed there. See code/check_r_deps.py.
if ! uv run --quiet python "$SCRIPT_DIR/code/check_r_deps.py"; then
    echo "ERROR: R dependency preflight failed. Nothing was run." >&2
    exit 1
fi
echo

# --- Step inventory ----------------------------------------------------------
# There is no step 00. The plan to vendor 00_cohort.py from
# CLIF-epidemiology-of-CRRT was abandoned on 2026-08-16; code/vendor/ now pins two
# config files and no code. See README, "Pipeline steps".
PYTHON_STEPS=(
    "code/01_build_cohort.py"
    "code/02_build_lmtp_df.py"
)
R_SCRIPT="code/03_lmtp_fit.R"

MISSING=()
for step in "${PYTHON_STEPS[@]}" "$R_SCRIPT"; do
    [ -f "$SCRIPT_DIR/$step" ] || MISSING+=("$step")
done

if [ ${#MISSING[@]} -gt 0 ]; then
    echo "--------------------------------------------------------------"
    echo " PIPELINE INCOMPLETE. Nothing was run."
    echo "--------------------------------------------------------------"
    echo "Missing step(s):"
    for m in "${MISSING[@]}"; do echo "  - $m"; done
    echo
    echo "Preflight itself passed: config is valid and the environment resolves."
    exit 2
fi

# --- Execution ---------------------------------------------------------------
START=$SECONDS

for step in "${PYTHON_STEPS[@]}"; do
    echo ">>> $step"
    STEP_START=$SECONDS
    if ! uv run python "$SCRIPT_DIR/$step"; then
        echo >&2
        echo "--------------------------------------------------------------" >&2
        echo " FAILED: $step" >&2
        echo "--------------------------------------------------------------" >&2
        if [ "$step" = "code/02_build_lmtp_df.py" ]; then
            echo "If this failed on medication unit conversion, that is the KNOWN" >&2
            echo "NEE blocker and the halt is deliberate, not a crash." >&2
            echo "clifpy leaves the RAW value in med_dose_converted when it cannot" >&2
            echo "convert, so an unconverted vasopressor silently inflates the" >&2
            echo "norepinephrine equivalent. Step 02 raises rather than continue." >&2
            echo "Dropping the rows is NOT a safe workaround: it removes a drug" >&2
            echo "from NEE for exactly the patients who received it." >&2
        fi
        exit 1
    fi
    echo "    done in $((SECONDS - STEP_START))s"
done

if ! command -v Rscript >/dev/null 2>&1; then
    echo >&2
    echo "WARNING: Rscript not on PATH, so the smoke fit was skipped." >&2
    echo "         The analysis frame was still built." >&2
    echo "         Install R 4.3+ and run: Rscript $R_SCRIPT smoke" >&2
    exit 3
fi

echo ">>> $R_SCRIPT smoke"
STEP_START=$SECONDS
Rscript "$SCRIPT_DIR/$R_SCRIPT" smoke
echo "    done in $((SECONDS - STEP_START))s"

ELAPSED=$((SECONDS - START))
echo
echo "--------------------------------------------------------------"
echo " Frame built and smoke fit passed in $((ELAPSED / 60))m $((ELAPSED % 60))s."
echo "--------------------------------------------------------------"
echo "The remaining stages are run by hand, on purpose:"
echo
echo "  Rscript $R_SCRIPT gate     # diagnostics ONLY, no effect estimate"
echo "  # ...read the diagnostics, then:"
echo "  Rscript $R_SCRIPT expand   # the full delta ladder"
echo
echo "Shareable outputs: output/final_no_phi/  (PHI-check before sending)"
