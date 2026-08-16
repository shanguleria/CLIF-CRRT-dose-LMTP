#!/usr/bin/env bash
# CRRT-dose-lmtp pipeline runner (macOS / Linux).
#
# SCAFFOLD. No analysis step exists yet, so this script performs the preflight
# checks and then reports exactly what is missing rather than pretending to run.
# It exits non-zero while the pipeline is incomplete: a runner that exits 0
# without doing anything is how a site ends up believing it has results.
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

# --- Step inventory ----------------------------------------------------------
PYTHON_STEPS=(
    "code/vendor/00_cohort.py"
    "code/vendor/01_create_wide_df.py"
    "code/02_build_lmtp_df.py"
)
R_STEPS=(
    "code/03_lmtp_fit.R"
)

MISSING=()
for step in "${PYTHON_STEPS[@]}" "${R_STEPS[@]}"; do
    [ -f "$SCRIPT_DIR/$step" ] || MISSING+=("$step")
done

if [ ${#MISSING[@]} -gt 0 ]; then
    echo "--------------------------------------------------------------"
    echo " PIPELINE INCOMPLETE. Nothing was run."
    echo "--------------------------------------------------------------"
    echo "Missing step(s):"
    for m in "${MISSING[@]}"; do echo "  - $m"; done
    echo
    echo "Steps 00 and 01 are vendored, not written here:"
    echo "  bash code/vendor/sync_vendor.sh"
    echo "Steps 02 and 03 are not yet implemented. See .claude/claude-todo.md."
    echo
    echo "Preflight itself passed: config is valid and the environment resolves."
    exit 2
fi

# --- Execution (reached only once every step exists) -------------------------
START=$SECONDS
for step in "${PYTHON_STEPS[@]}"; do
    echo ">>> $step"
    STEP_START=$SECONDS
    uv run python "$SCRIPT_DIR/$step"
    echo "    done in $((SECONDS - STEP_START))s"
done

if command -v Rscript >/dev/null 2>&1; then
    for step in "${R_STEPS[@]}"; do
        echo ">>> $step"
        STEP_START=$SECONDS
        Rscript "$SCRIPT_DIR/$step"
        echo "    done in $((SECONDS - STEP_START))s"
    done
else
    echo "WARNING: Rscript not on PATH; R steps skipped." >&2
fi

echo "Total: $(( (SECONDS - START) / 60 ))m $(( (SECONDS - START) % 60 ))s"
echo "Shareable outputs: output/final_no_phi/  (PHI-check before sending)"
