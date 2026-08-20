"""Which site a run belongs to, and where that site's outputs go.

ONE resolution rule, in ONE place, used by every step.

Before this module existed the rule was written down three times and disagreed with
itself. `run_pipeline.sh` resolved `CLIF_CONFIG` and used it for preflight, while
`01_build_cohort.py`, `02_build_lmtp_df.py` and `03_lmtp_fit.R` each hardcoded
`config/config.json`. So `CLIF_CONFIG=config/config_SiteA.json ./run_pipeline.sh` printed
`Site : SiteA`, checked SiteA's `has_crrt_settings`, and then ran the other site's config
against the other site's data. The run exited 0 and its output looked fine.

That is why the equality check below is not optional politeness. A wrong-site run does
not announce itself: it produces a full set of plausible artifacts under the wrong name.

`code/03_lmtp_fit.R` carries an R transcription of `resolve_config`'s precedence and of
the site_name equality check, because R cannot import this file. The two are a matched
pair; changing the rule here means changing it there. `tests/test_site_resolution.py`
drives the R script and asserts the two agree, so the pair is checked by execution
rather than by hope.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

# A site_name becomes a directory name under output/. Validate it BEFORE it is used to
# build a path: "../../.." is a perfectly good JSON string, and os.path.join would
# cheerfully honour it.
SITE_NAME_RE = re.compile(r"^[A-Za-z0-9_-]+$")

# Written into every site's intermediate_phi/ when the directory is created. This used to
# be a tracked file at output/intermediate_phi/README.md, which meant it survived only
# until someone ran `git clean -fdx`. Site directories are created at runtime and cannot
# be tracked, so the label is now written at runtime too, which is strictly more reliable.
PHI_LABEL = """# intermediate_phi

Patient-level intermediates: node datasets, fitted objects, diagnostics keyed to
individuals.

**This directory never leaves the site.** It is not part of any export bundle and
must not be committed, copied into a shared drive, or read into an analysis
transcript. Empty until the pipeline is built.
"""


def available_sites(repo_root):
    """Site names this repo has a config for, from config/config_<NAME>.json.

    The name comes from the FILENAME, not from site_name inside the file, because the
    filename is what --site takes. resolve_config then requires the two to agree, so a
    name listed here is a name that will actually work.
    """
    cfg_dir = Path(repo_root) / "config"
    if not cfg_dir.is_dir():
        return []
    names = []
    for p in sorted(cfg_dir.glob("config_*.json")):
        name = p.stem[len("config_"):]
        if name == "template":
            continue
        if SITE_NAME_RE.match(name):
            names.append(name)
    return names


def _invocation_hint():
    """How the caller should be re-run, e.g. 'uv run python code/01_build_cohort.py --site'.

    The listing appends a site name to whatever this returns, so it must end at the
    flag and not carry a placeholder of its own.
    """
    import sys
    prog = Path(sys.argv[0]).name if sys.argv and sys.argv[0] else ""
    if prog.endswith(".py") and prog != "_site.py":
        return f"uv run python code/{prog} --site"
    return "--site"


def missing_default_message(repo_root, how=None):
    """What to say when config/config.json is absent.

    Two different readers hit this. A consortium site setting up for the first time
    needs "create the file". A machine that deliberately has no default, because it
    holds several sites, needs "name one of these" -- telling that reader to create a
    default is advice against their own setup. Which one is speaking is decided by
    whether any per-site config exists.
    """
    how = _invocation_hint() if how is None else how
    sites = available_sites(repo_root)
    if not sites:
        return ("ERROR: config not found at "
                f"{Path(repo_root) / 'config' / 'config.json'}\n"
                "       Fix: cp config/config_template.json config/config.json\n"
                "       then edit it for your site.")
    listed = "\n".join(f"         {how} {s}" for s in sites)
    return ("ERROR: no default config at config/config.json\n"
            "       This machine has per-site configs. Name one:\n"
            f"{listed}\n"
            "       (or create config/config.json for a single-site setup)")


def parse_site(argv):
    """Pull `--site NAME` out of argv. Returns (site_or_None, remaining_args).

    Kept separate from resolve_config so 03's R equivalent has an obvious counterpart,
    and so a caller that has its own positional arguments (the R script's STAGE) can
    strip the pair before reading them.
    """
    site, rest = None, []
    argv = list(argv)
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok == "--site":
            if i + 1 >= len(argv):
                raise SystemExit("--site requires a site name, e.g. --site SiteA")
            site = argv[i + 1]
            i += 2
        elif tok.startswith("--site="):
            site = tok.split("=", 1)[1]
            i += 1
        else:
            rest.append(tok)
            i += 1
    return site, rest


def resolve_config(repo_root, argv=None):
    """Resolve which config this run uses. Returns (path, parsed_config).

    Precedence, matched exactly by run_pipeline.sh, run_pipeline.ps1 and 03_lmtp_fit.R:

      1. --site NAME   ->  config/config_<NAME>.json
      2. $CLIF_CONFIG  ->  that path verbatim
      3. default       ->  config/config.json

    A single consortium site sets none of these and lands on (3), which is the workflow
    the README documents and this change deliberately leaves alone.
    """
    repo_root = Path(repo_root)
    argv = [] if argv is None else argv
    site, _ = parse_site(argv)

    if site is not None:
        _validate_site_name(site, source="--site")
        path = repo_root / "config" / f"config_{site}.json"
        hint = (f"Fix: cp config/config_template.json config/config_{site}.json\n"
                f"     then edit it for {site}.")
    elif os.environ.get("CLIF_CONFIG"):
        path = Path(os.environ["CLIF_CONFIG"])
        if not path.is_absolute():
            path = repo_root / path
        hint = "Fix: point CLIF_CONFIG at a config that exists, or unset it."
    else:
        path = repo_root / "config" / "config.json"
        # Site-aware, so the message suits whichever kind of machine this is. On one
        # holding several sites, "create config.json" is advice against the operator's
        # own setup; the real fix there is to name a site.
        if not path.is_file():
            raise SystemExit(missing_default_message(repo_root))
        hint = None

    if hint is not None and not path.is_file():
        raise SystemExit(f"ERROR: config not found at {path}\n       {hint}")

    config = json.loads(path.read_text())

    if "site_name" not in config:
        raise SystemExit(f"ERROR: {path} has no 'site_name'. Every output is stamped and "
                         f"filed under it, so a run without one cannot be attributed.")
    _validate_site_name(config["site_name"], source=f"site_name in {path.name}")

    # The check this module exists for. Without it, config_SiteA.json could declare
    # site_name "SiteB" and SiteA's results would be written into SiteB's directory,
    # stamped SiteB, with nothing anywhere reporting a problem.
    if site is not None and config["site_name"] != site:
        raise SystemExit(
            f"ERROR: site mismatch.\n"
            f"       You asked for --site {site}, which resolved to {path.name},\n"
            f"       but that file declares site_name = {config['site_name']!r}.\n"
            f"       Refusing to run: outputs would be filed and stamped under "
            f"{config['site_name']!r}, not {site!r}.\n"
            f"       Fix whichever is wrong, the filename or the site_name.")

    return path, config


def site_output_root(repo_root, config):
    """output/<site_name>/ — the root of everything this run may write.

    Always nested, including for a site running alone. One rule beats two, and
    output/SiteB/ is self-describing in a way output/ is not.
    """
    return Path(repo_root) / "output" / config["site_name"]


def ensure_site_dirs(site_out):
    """Create this site's output tree and return (intermediate_phi, final_no_phi)."""
    site_out = Path(site_out)
    phi = site_out / "intermediate_phi"
    share = site_out / "final_no_phi"
    phi.mkdir(parents=True, exist_ok=True)
    share.mkdir(parents=True, exist_ok=True)

    label = phi / "README.md"
    if not label.exists():
        label.write_text(PHI_LABEL)

    return phi, share


def _validate_site_name(name, source):
    if not isinstance(name, str) or not SITE_NAME_RE.match(name):
        raise SystemExit(
            f"ERROR: invalid site name {name!r} (from {source}).\n"
            f"       A site name becomes a directory under output/, so it is restricted "
            f"to letters, digits, underscore and hyphen.")


# ---------------------------------------------------------------------------
# A tiny CLI, so the runners do not carry a third copy of the rule.
#
#   python code/_site.py --explain-missing "./run_pipeline.sh"
#
# Both runners already shell out to python for site_name, so delegating this costs
# them nothing and keeps available_sites() the single place that knows how a per-site
# config is named.
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    _repo = Path(__file__).resolve().parent.parent
    if len(sys.argv) >= 2 and sys.argv[1] == "--explain-missing":
        how = sys.argv[2] if len(sys.argv) > 2 else "--site NAME"
        print(missing_default_message(_repo, how))
        raise SystemExit(0)
    if len(sys.argv) >= 2 and sys.argv[1] == "--list-sites":
        print("\n".join(available_sites(_repo)))
        raise SystemExit(0)
    if len(sys.argv) >= 2 and sys.argv[1] == "--check":
        # Answer "is this config actually usable?" in a second, rather than after a
        # table read. Reports only paths and settings; it opens no CLIF table and
        # prints no file listing, so it stays clear of patient data entirely.
        _argv = ["--site", sys.argv[2]] if len(sys.argv) > 2 else []
        _path, _cfg = resolve_config(_repo, _argv)
        _out = site_output_root(_repo, _cfg)
        _data = Path(_cfg.get("data_directory", ""))
        _problems = []
        if not _cfg.get("has_crrt_settings"):
            _problems.append("has_crrt_settings is not true: dose is the exposure and "
                             "cannot be computed without CRRT flow rates.")
        if not str(_data) or "/path/to/" in str(_data):
            _problems.append(f"data_directory is still a placeholder ({_data}).")
        elif not _data.is_dir():
            _problems.append(f"data_directory does not exist: {_data}")
        print(f"config          {_path}")
        print(f"site_name       {_cfg['site_name']}")
        print(f"clif_version    {_cfg.get('clif_version')}")
        print(f"filetype        {_cfg.get('filetype')}")
        print(f"timezone        {_cfg.get('timezone')}")
        print(f"n_workers       {_cfg.get('n_workers')}")
        print(f"data_directory  {_data}")
        print(f"output would go {_out}")
        if _problems:
            print("\nNOT READY:")
            for _p in _problems:
                print(f"  - {_p}")
            raise SystemExit(1)
        print("\nREADY.")
        raise SystemExit(0)
    raise SystemExit("usage: _site.py [--explain-missing HOW | --list-sites | --check SITE]")
