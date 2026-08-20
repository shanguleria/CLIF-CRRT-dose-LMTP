#!/usr/bin/env python3
"""Assert that every R package the code loads is recorded in renv.lock.

This is a LOCKFILE-COMPLETENESS check, not an environment check. It answers one
question: if a site clones this repo and runs `renv::restore()`, will every
package that `code/*.R` calls `library()` on actually be installed?

It needs neither R nor an installed library, because the failure it catches is a
maintainer's, not a site's. `nanoparquet` was used at `03_lmtp_fit.R:33` from the
day that script was written and was never snapshotted; the coordinating machine
ran fine because the package happened to be installed there already. Every fresh
clone would have died at line 33, after steps 01 and 02 had already burned an
hour. That is the class of bug this exists to make impossible.

The check is deliberately one-directional. Used-but-not-locked is fatal. The
reverse, locked-but-not-used, is ordinary transitive-dependency noise: 74 of the
75 entries in this lockfile are dependencies of dependencies and no R file names
them. Failing on those would make the check unusable.

Usage:
    python code/check_r_deps.py           # exit 1 on any missing package
    python code/check_r_deps.py --list    # print the resolved dependency set
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
R_SOURCE_DIR = REPO_ROOT / "code"
LOCKFILE = REPO_ROOT / "renv.lock"

# Base and "priority: base" packages ship with R itself and are never recorded in
# a lockfile. Recommended packages (MASS, Matrix, survival, ...) ARE recorded, so
# they must not appear here.
BASE_PACKAGES = frozenset(
    {
        "base",
        "compiler",
        "datasets",
        "grDevices",
        "graphics",
        "grid",
        "methods",
        "parallel",
        "splines",
        "stats",
        "stats4",
        "tcltk",
        "tools",
        "utils",
    }
)

LIBRARY_CALL = re.compile(
    r"""\b(?:library|require)\s*\(\s*["']?([A-Za-z][A-Za-z0-9._]*)["']?\s*[,)]"""
)
REQUIRE_NAMESPACE = re.compile(r"""\brequireNamespace\s*\(\s*["']([A-Za-z][A-Za-z0-9._]*)["']""")
NAMESPACED_CALL = re.compile(r"""\b([A-Za-z][A-Za-z0-9._]*)::""")


def strip_r_comment(line: str) -> str:
    """Drop everything from the first unquoted '#' to end of line.

    Quote-aware rather than a plain split, because a false positive here would
    block a site's run over a package name that appears only in prose.
    """
    quote = None
    escaped = False
    for i, ch in enumerate(line):
        if escaped:
            escaped = False
            continue
        if ch == "\\" and quote is not None:
            escaped = True
        elif quote is not None:
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
        elif ch == "#":
            return line[:i]
    return line


def r_source_files(source_dir: Path) -> list[Path]:
    return sorted(source_dir.rglob("*.R"))


def scan_r_dependencies(source_dir: Path) -> dict[str, list[str]]:
    """Map each package named by R source to the `file:line` sites that name it."""
    found: dict[str, list[str]] = {}
    for path in r_source_files(source_dir):
        rel = path.relative_to(REPO_ROOT)
        for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            code = strip_r_comment(raw)
            names = set()
            for pattern in (LIBRARY_CALL, REQUIRE_NAMESPACE, NAMESPACED_CALL):
                names.update(pattern.findall(code))
            for name in names - BASE_PACKAGES:
                found.setdefault(name, []).append(f"{rel}:{lineno}")
    return found


def locked_packages(lockfile: Path) -> set[str]:
    with lockfile.open(encoding="utf-8") as fh:
        return set(json.load(fh).get("Packages", {}))


def main(argv: list[str]) -> int:
    if not LOCKFILE.exists():
        print(f"ERROR: no lockfile at {LOCKFILE}", file=sys.stderr)
        return 1

    r_files = r_source_files(R_SOURCE_DIR)
    used = scan_r_dependencies(R_SOURCE_DIR)
    locked = locked_packages(LOCKFILE)

    # A check that can pass by finding nothing is not a check. Renaming or moving
    # the R script would otherwise turn this green forever, which is the exact
    # shape of the bug it exists to catch. Seen for real: with 03_lmtp_fit.R
    # temporarily renamed, an earlier version of this script reported
    # "0 package(s) used ... all recorded" and exited 0.
    if not r_files:
        print(
            f"ERROR: no .R files found under {R_SOURCE_DIR}. This check would pass\n"
            "       vacuously, so it fails instead. Has the R source moved?",
            file=sys.stderr,
        )
        return 1
    if not used:
        print(
            f"ERROR: scanned {len(r_files)} .R file(s) and found no package "
            "dependencies at all.\n"
            "       Expected at least library(lmtp). The scan is probably broken.",
            file=sys.stderr,
        )
        return 1

    if "--list" in argv:
        for name in sorted(used):
            mark = "ok     " if name in locked else "MISSING"
            print(f"  {mark}  {name:<16} {used[name][0]}")
        print(f"\n{len(used)} package(s) used, {len(locked)} recorded in renv.lock")

    missing = {name: sites for name, sites in used.items() if name not in locked}
    if not missing:
        if "--list" not in argv:
            print(
                f"R dependency check: {len(used)} package(s) across {len(r_files)} "
                f".R file(s), all recorded in renv.lock ({len(locked)} entries)."
            )
        return 0

    print(
        "\nERROR: R package(s) used by this repo's code but ABSENT from renv.lock.\n"
        "       A fresh clone will restore an incomplete environment and step 03\n"
        "       will fail at the library() call, after steps 01 and 02 have run.\n",
        file=sys.stderr,
    )
    for name, sites in sorted(missing.items()):
        print(f"  {name}", file=sys.stderr)
        for site in sites:
            print(f"      used at {site}", file=sys.stderr)
    print(
        "\n       Fix, from a machine where the package IS installed:\n"
        "         Rscript -e 'renv::record(list(<pkg> = list("
        'Package = "<pkg>",\n'
        '           Version = "<version>", Source = "Repository", '
        'Repository = "CRAN")))\'\n'
        "       then commit renv.lock. Prefer record() over snapshot() so the diff\n"
        "       is exactly one block.\n",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
