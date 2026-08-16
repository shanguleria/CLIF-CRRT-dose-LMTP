"""Assert every vendored file is byte-identical to its pinned upstream commit.

Vendored code is copied, not imported, so nothing stops someone from "just
fixing" a line in code/vendor/ and silently forking the cohort definition from
the sibling CJASN analysis. This test is the guard: drift becomes a failing
check instead of two papers reporting two different Ns.

The manifest lives in code/vendor/VENDOR_SHA. While it lists nothing but the
config files, the code-vendoring test skips with an explicit reason rather than
passing vacuously.
"""

from __future__ import annotations

import hashlib
import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
VENDOR_SHA_FILE = REPO_ROOT / "code" / "vendor" / "VENDOR_SHA"

# Paths in the manifest that are config rather than the ~6,100 lines of pipeline
# code. Kept separate so an empty *code* manifest still reports as "not yet
# vendored" instead of hiding behind two green config checks.
CONFIG_PREFIXES = ("config/",)


def _parse_vendor_sha() -> tuple[str, Path, list[tuple[str, str]]]:
    if not VENDOR_SHA_FILE.exists():
        pytest.fail(f"{VENDOR_SHA_FILE} is missing; the vendor pin is the contract.")

    text = VENDOR_SHA_FILE.read_text()

    sha_match = re.search(r"^SHA=([0-9a-f]{40})$", text, re.MULTILINE)
    if not sha_match:
        pytest.fail("VENDOR_SHA has no full 40-character SHA= line.")
    sha = sha_match.group(1)

    up_match = re.search(r"^UPSTREAM_REPO=(.+)$", text, re.MULTILINE)
    if not up_match:
        pytest.fail("VENDOR_SHA has no UPSTREAM_REPO= line.")
    upstream = (VENDOR_SHA_FILE.parent / up_match.group(1).strip()).resolve()

    manifest = [
        (m.group(1), m.group(2))
        for m in re.finditer(r"^MANIFEST\s+(\S+)\s+(\S+)\s*$", text, re.MULTILINE)
    ]
    return sha, upstream, manifest


def _upstream_blob(upstream: Path, sha: str, path: str) -> bytes:
    result = subprocess.run(
        ["git", "-C", str(upstream), "show", f"{sha}:{path}"],
        capture_output=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"cannot read {path} at {sha[:12]} from {upstream}: "
            f"{result.stderr.decode(errors='replace').strip()}"
        )
    return result.stdout


SHA, UPSTREAM, MANIFEST = _parse_vendor_sha()

CODE_ENTRIES = [e for e in MANIFEST if not e[1].startswith(CONFIG_PREFIXES)]


def _upstream_available() -> bool:
    if not (UPSTREAM / ".git").exists():
        return False
    return (
        subprocess.run(
            ["git", "-C", str(UPSTREAM), "cat-file", "-e", f"{SHA}^{{commit}}"],
            capture_output=True,
        ).returncode
        == 0
    )


requires_upstream = pytest.mark.skipif(
    not _upstream_available(),
    reason=(
        f"upstream repo unavailable at {UPSTREAM} or commit {SHA[:12]} not fetched. "
        "This check only runs where the sibling repo is checked out (i.e. the "
        "coordinating site), not at participating sites."
    ),
)


@requires_upstream
@pytest.mark.parametrize("upstream_path,local_path", MANIFEST, ids=lambda p: str(p))
def test_vendored_file_matches_pin(upstream_path: str, local_path: str) -> None:
    local = REPO_ROOT / local_path
    assert local.exists(), (
        f"{local_path} is in the manifest but missing on disk. "
        "Run: bash code/vendor/sync_vendor.sh"
    )

    expected = _upstream_blob(UPSTREAM, SHA, upstream_path)
    actual = local.read_bytes()

    if expected != actual:
        pytest.fail(
            f"{local_path} has DRIFTED from {upstream_path} @ {SHA[:12]}.\n"
            f"  expected sha256 {hashlib.sha256(expected).hexdigest()[:16]} "
            f"({len(expected)} bytes)\n"
            f"  actual   sha256 {hashlib.sha256(actual).hexdigest()[:16]} "
            f"({len(actual)} bytes)\n"
            "Vendored files are never edited in place. Either revert the edit, or "
            "move the pin in VENDOR_SHA and re-sync deliberately."
        )


@requires_upstream
def test_pipeline_code_is_vendored() -> None:
    """Fails once code vendoring is expected; skips clearly until then."""
    if not CODE_ENTRIES:
        pytest.skip(
            "No pipeline code vendored yet. The manifest in code/vendor/VENDOR_SHA "
            "lists the five files to bring across (00_cohort.py, sofa_calculator.py, "
            "01_create_wide_df.py, pipeline_helpers.py, utils.py) as commented-out "
            "MANIFEST lines. Uncomment them and run sync_vendor.sh; this test then "
            "activates automatically."
        )

    missing = [lp for _, lp in CODE_ENTRIES if not (REPO_ROOT / lp).exists()]
    assert not missing, f"manifest lists code not on disk: {missing}"
