"""Assert renv.lock records every R package this repo's code loads.

The bug this guards against shipped silently and was invisible on the machine
that wrote it: `03_lmtp_fit.R:33` has called `library(nanoparquet)` since the day
it was written, and nanoparquet appeared zero times in renv.lock. The
coordinating site ran fine because the package happened to be installed there.
Every consortium site cloning the repo would have hit a fatal error at that line,
one long step-01 table read and one step-02 build later.

A lockfile is only a contract if something enforces it. This is that something.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "code"))

from check_r_deps import (  # noqa: E402
    LOCKFILE,
    R_SOURCE_DIR,
    locked_packages,
    scan_r_dependencies,
    strip_r_comment,
)


def test_every_used_r_package_is_locked() -> None:
    used = scan_r_dependencies(R_SOURCE_DIR)
    locked = locked_packages(LOCKFILE)

    missing = {name: sites for name, sites in used.items() if name not in locked}
    assert not missing, (
        "R package(s) used by code but absent from renv.lock: "
        + "; ".join(f"{n} (used at {s[0]})" for n, s in sorted(missing.items()))
        + ". A fresh clone would restore an incomplete environment. Fix with "
        "renv::record(), not renv::snapshot(), so the diff is one block."
    )


def test_the_scan_actually_finds_something() -> None:
    """A check that can only ever pass vacuously is not a check.

    If the R sources move or the scan regexes break, the test above goes green by
    finding nothing at all. This pins the floor.
    """
    used = scan_r_dependencies(R_SOURCE_DIR)
    assert "lmtp" in used, f"scan found {sorted(used)}; expected at least lmtp"
    assert "nanoparquet" in used, "the regression case itself must still be detected"


def test_missing_package_is_detected(tmp_path: Path) -> None:
    """Prove the check goes RED, by reconstructing the exact bug."""
    truncated = json.loads(LOCKFILE.read_text())
    removed = truncated["Packages"].pop("nanoparquet", None)
    assert removed is not None, "nanoparquet must be in renv.lock for this test to mean anything"

    fake_lock = tmp_path / "renv.lock"
    fake_lock.write_text(json.dumps(truncated))

    used = scan_r_dependencies(R_SOURCE_DIR)
    locked = locked_packages(fake_lock)
    missing = [name for name in used if name not in locked]

    assert missing == ["nanoparquet"], f"expected exactly nanoparquet to be flagged, got {missing}"


def test_check_refuses_to_pass_vacuously(tmp_path: Path, monkeypatch) -> None:
    """An empty source tree must FAIL, not report "all recorded".

    Found by accident, not by design: with 03_lmtp_fit.R temporarily renamed, an
    earlier version of this check printed "0 package(s) used ... all recorded" and
    exited 0. Renaming the R script would have silently disarmed the guard.
    """
    import check_r_deps

    empty = tmp_path / "code"
    empty.mkdir()
    monkeypatch.setattr(check_r_deps, "R_SOURCE_DIR", empty)

    assert check_r_deps.main([]) == 1, "an empty R source tree must exit non-zero"


@pytest.mark.parametrize(
    "line,expected",
    [
        ("library(real)", "library(real)"),
        ("# library(ghost)", ""),
        ("x <- 1  # library(ghost)", "x <- 1  "),
        ('f("# not a comment")  # real', 'f("# not a comment")  '),
        ("g('a#b') # tail", "g('a#b') "),
    ],
)
def test_comment_stripping_is_quote_aware(line: str, expected: str) -> None:
    """A package named only in prose must not block a site's run."""
    assert strip_r_comment(line) == expected
