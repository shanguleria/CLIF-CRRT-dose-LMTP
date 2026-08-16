#!/usr/bin/env bash
# Pull every file in the VENDOR_SHA manifest from the pinned upstream commit.
#
# Vendoring, not rewriting, is a deliberate choice. 00_cohort.py is where the
# science lives and it carries most of the correctness commits in the upstream
# fix/pipeline-corrections branch: the ESRD gate, CRRT initiation per encounter
# block, the modality-agnostic effluent formula, the first-3h dose median, the
# 30-day anchor to CRRT initiation, bit-reproducibility. Reimplementing it forks
# the cohort definition, and then two papers from the same consortium report two
# different Ns. See .claude/lmtp_feasibility_findings.md S6.2.
#
# Usage: bash code/vendor/sync_vendor.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
VENDOR_SHA_FILE="$REPO_ROOT/code/vendor/VENDOR_SHA"

[ -f "$VENDOR_SHA_FILE" ] || { echo "ERROR: $VENDOR_SHA_FILE not found" >&2; exit 1; }

SHA=$(grep -E '^SHA=' "$VENDOR_SHA_FILE" | head -1 | cut -d= -f2)
UPSTREAM_REL=$(grep -E '^UPSTREAM_REPO=' "$VENDOR_SHA_FILE" | head -1 | cut -d= -f2)
UPSTREAM=$(cd "$REPO_ROOT/code/vendor" && cd "$UPSTREAM_REL" 2>/dev/null && pwd) || {
    echo "ERROR: upstream repo not found at $UPSTREAM_REL (relative to code/vendor/)" >&2
    echo "       Clone it beside this repo, or edit UPSTREAM_REPO in VENDOR_SHA." >&2
    exit 1
}

echo "Upstream : $UPSTREAM"
echo "Pinned   : $SHA"

git -C "$UPSTREAM" cat-file -e "${SHA}^{commit}" 2>/dev/null || {
    echo "ERROR: commit $SHA is not present in $UPSTREAM. Fetch it first." >&2
    exit 1
}

n=0
while read -r _ up local; do
    [ -n "${up:-}" ] || continue
    mkdir -p "$REPO_ROOT/$(dirname "$local")"
    git -C "$UPSTREAM" show "${SHA}:${up}" > "$REPO_ROOT/$local"
    echo "  vendored  $up  ->  $local"
    n=$((n + 1))
done < <(grep -E '^MANIFEST[[:space:]]' "$VENDOR_SHA_FILE")

echo "$n file(s) synced. Now run: uv run pytest tests/test_vendor_integrity.py"
