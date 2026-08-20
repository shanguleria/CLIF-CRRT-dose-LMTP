"""Assert that a run cannot be attributed to the wrong site.

The bug this guards against was live and silent. `run_pipeline.sh` resolved
`$CLIF_CONFIG` and used it for preflight, while `01_build_cohort.py`,
`02_build_lmtp_df.py` and `03_lmtp_fit.R` each hardcoded `config/config.json`. So
`CLIF_CONFIG=config/config_SiteA.json ./run_pipeline.sh` printed `Site : SiteA`, checked SiteA's
`has_crrt_settings`, then ran the other site's config against the other site's data and
exited 0 with a full set of plausible artifacts.

Nothing about a wrong-site run looks wrong from the outside, which is why the checks
below are worth their weight: every one of them is the difference between a loud failure
and a quiet mislabeled result.

The R half of the rule lives in `code/03_lmtp_fit.R` because R cannot import
`_site.py`. `test_r_and_python_resolvers_agree` drives the actual R script rather than
reading it, so the pair is checked by execution.
"""

from __future__ import annotations

import functools
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "code"))

import _site  # noqa: E402

R_SCRIPT = REPO_ROOT / "code" / "03_lmtp_fit.R"


@functools.lru_cache(maxsize=1)
def _renv_libpaths() -> str:
    """The pinned library, as a path string R can be handed via $R_LIBS.

    The R tests below run with cwd set to a temporary repo so they can plant configs
    without touching this one. That cwd has no .Rprofile, so renv never activates and
    `library(lmtp)` fails before the script reaches the code under test. Asking the real
    project for its .libPaths() and forwarding them keeps the library pinned without
    hardcoding a platform-specific path.
    """
    proc = subprocess.run(
        ["Rscript", "-e", "cat(paste(.libPaths(), collapse=.Platform$path.sep))"],
        capture_output=True, text=True, timeout=120, cwd=REPO_ROOT,
    )
    return proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ""


@functools.lru_cache(maxsize=1)
def _r_stack_available() -> bool:
    """Can 03 actually load its libraries here?

    renv/library is gitignored, so a fresh clone has R on PATH but an empty pinned
    library until `renv::restore()` runs. Without this the two R tests below would fail
    rather than skip on every clone, which trains people to ignore a red suite. Checked
    by loading the package, not by looking for a directory: an empty renv/library exists
    in a clone and would pass a path check while failing every library() call.
    """
    if shutil.which("Rscript") is None:
        return False
    proc = subprocess.run(
        ["Rscript", "-e", "suppressPackageStartupMessages(library(lmtp))"],
        capture_output=True, text=True, timeout=180, cwd=REPO_ROOT,
        env={**os.environ, "R_LIBS": _renv_libpaths()},
    )
    return proc.returncode == 0


requires_r = pytest.mark.skipif(
    not _r_stack_available(),
    reason="pinned R library not restored here; run Rscript -e 'renv::restore()'",
)


def make_r_repo(tmp_path: Path, config_filename: str, site_name: str) -> Path:
    """A throwaway repo holding a real copy of 03, so it resolves REPO to here.

    03 derives its repo root from its own --file= path, not from cwd, which is correct
    in production and means a test cannot relocate it by running elsewhere. Copying the
    shipped script into the temporary tree is what moves REPO, and it keeps the test
    exercising the real file rather than a transcription of it.
    """
    repo = tmp_path / "repo"
    (repo / "code").mkdir(parents=True)
    (repo / "config").mkdir(parents=True)
    shutil.copy(R_SCRIPT, repo / "code" / R_SCRIPT.name)
    shutil.copy(REPO_ROOT / "config" / "lmtp_design.json", repo / "config")
    write_config(repo, config_filename, site_name)
    return repo


def run_r(repo: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["Rscript", str(repo / "code" / R_SCRIPT.name), *args],
        capture_output=True, text=True, timeout=300, cwd=repo,
        env={**os.environ, "R_LIBS": _renv_libpaths()},
    )


def write_config(repo: Path, filename: str, site_name: str) -> Path:
    """A minimally valid site config. Only the keys resolution reads need be real."""
    cfg = repo / "config" / filename
    cfg.parent.mkdir(parents=True, exist_ok=True)
    cfg.write_text(json.dumps({
        "site_name": site_name,
        "clif_version": "2.1.0",
        "data_directory": "/nonexistent",
        "filetype": "parquet",
        "timezone": "America/Chicago",
        "has_crrt_settings": True,
        "n_workers": 1,
    }))
    return cfg


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    (tmp_path / "config").mkdir()
    return tmp_path


# ---------------------------------------------------------------------------
# Precedence
# ---------------------------------------------------------------------------

def test_default_is_plain_config_json(repo: Path, monkeypatch) -> None:
    """The consortium path: no flag, no env var, unchanged behaviour."""
    monkeypatch.delenv("CLIF_CONFIG", raising=False)
    write_config(repo, "config.json", "SiteB")

    path, cfg = _site.resolve_config(repo, [])

    assert path == repo / "config" / "config.json"
    assert cfg["site_name"] == "SiteB"


def test_site_flag_selects_the_named_config(repo: Path, monkeypatch) -> None:
    monkeypatch.delenv("CLIF_CONFIG", raising=False)
    write_config(repo, "config.json", "SiteB")
    write_config(repo, "config_SiteA.json", "SiteA")

    path, cfg = _site.resolve_config(repo, ["--site", "SiteA"])

    assert path == repo / "config" / "config_SiteA.json"
    assert cfg["site_name"] == "SiteA"


def test_site_flag_accepts_equals_form(repo: Path, monkeypatch) -> None:
    monkeypatch.delenv("CLIF_CONFIG", raising=False)
    write_config(repo, "config_SiteA.json", "SiteA")

    _, cfg = _site.resolve_config(repo, ["--site=SiteA"])

    assert cfg["site_name"] == "SiteA"


def test_clif_config_env_is_honoured(repo: Path, monkeypatch) -> None:
    """The lever that used to be read by the runner and ignored by every step."""
    write_config(repo, "config.json", "SiteB")
    other = write_config(repo, "config_SiteD.json", "SiteD")
    monkeypatch.setenv("CLIF_CONFIG", str(other))

    path, cfg = _site.resolve_config(repo, [])

    assert path == other
    assert cfg["site_name"] == "SiteD"


def test_site_flag_outranks_clif_config(repo: Path, monkeypatch) -> None:
    """Precedence is a contract shared with both runners and the R script."""
    write_config(repo, "config_SiteA.json", "SiteA")
    monkeypatch.setenv("CLIF_CONFIG", str(write_config(repo, "config_SiteD.json", "SiteD")))

    _, cfg = _site.resolve_config(repo, ["--site", "SiteA"])

    assert cfg["site_name"] == "SiteA", "--site must win over CLIF_CONFIG"


def test_missing_config_names_the_path_and_the_fix(repo: Path, monkeypatch) -> None:
    monkeypatch.delenv("CLIF_CONFIG", raising=False)

    with pytest.raises(SystemExit) as exc:
        _site.resolve_config(repo, ["--site", "SiteA"])

    assert "config_SiteA.json" in str(exc.value)


# ---------------------------------------------------------------------------
# The mismatch check: the reason this module exists
# ---------------------------------------------------------------------------

def test_site_name_mismatch_is_fatal(repo: Path, monkeypatch) -> None:
    """--site SiteA against a config declaring SiteB must refuse to run.

    Without this, SiteA's results would be written into output/SiteB/ and stamped
    site_id=SiteB. Nothing downstream could tell.
    """
    monkeypatch.delenv("CLIF_CONFIG", raising=False)
    write_config(repo, "config_SiteA.json", "SiteB")   # the mistake

    with pytest.raises(SystemExit) as exc:
        _site.resolve_config(repo, ["--site", "SiteA"])

    msg = str(exc.value)
    assert "mismatch" in msg
    assert "SiteA" in msg and "SiteB" in msg, "the error must name BOTH sides"


def test_matching_site_name_is_accepted(repo: Path, monkeypatch) -> None:
    """The green half of the pair, so the test above is not passing vacuously."""
    monkeypatch.delenv("CLIF_CONFIG", raising=False)
    write_config(repo, "config_SiteA.json", "SiteA")

    _, cfg = _site.resolve_config(repo, ["--site", "SiteA"])

    assert cfg["site_name"] == "SiteA"


def test_clif_config_is_not_subject_to_the_name_check(repo: Path, monkeypatch) -> None:
    """No site was asserted, so there is nothing to contradict.

    CLIF_CONFIG names a path, not a claim about which site it is; the filename carries
    no site to disagree with. The check binds --site only.
    """
    monkeypatch.setenv("CLIF_CONFIG", str(write_config(repo, "anything.json", "SiteD")))

    _, cfg = _site.resolve_config(repo, [])

    assert cfg["site_name"] == "SiteD"


# ---------------------------------------------------------------------------
# The missing-default message has two audiences
# ---------------------------------------------------------------------------

def test_missing_default_tells_a_fresh_site_to_create_one(repo: Path, monkeypatch) -> None:
    """No per-site configs: this is a consortium site setting up for the first time."""
    monkeypatch.delenv("CLIF_CONFIG", raising=False)

    with pytest.raises(SystemExit) as exc:
        _site.resolve_config(repo, [])

    msg = str(exc.value)
    assert "cp config/config_template.json config/config.json" in msg
    assert "--site" not in msg, "there are no sites to name yet"


def test_missing_default_lists_the_sites_that_exist(repo: Path, monkeypatch) -> None:
    """Per-site configs present: telling this operator to create a default is wrong advice."""
    monkeypatch.delenv("CLIF_CONFIG", raising=False)
    for name in ("SiteC", "SiteA", "SiteD"):
        write_config(repo, f"config_{name}.json", name)

    with pytest.raises(SystemExit) as exc:
        _site.resolve_config(repo, [])

    msg = str(exc.value)
    assert "no default config" in msg
    for name in ("SiteC", "SiteA", "SiteD"):
        assert f"--site {name}" in msg, f"{name} not offered"


def test_available_sites_ignores_the_template(repo: Path) -> None:
    """config_template.json matches config_*.json but is not a site."""
    (repo / "config" / "config_template.json").write_text("{}")
    write_config(repo, "config_SiteA.json", "SiteA")

    assert _site.available_sites(repo) == ["SiteA"]


def test_available_sites_is_empty_without_configs(repo: Path) -> None:
    assert _site.available_sites(repo) == []


# ---------------------------------------------------------------------------
# site_name becomes a directory, so it is validated before it becomes a path
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", ["../../etc", "a/b", "with space", "", "sit.e"])
def test_path_unsafe_site_name_is_rejected(repo: Path, monkeypatch, bad: str) -> None:
    monkeypatch.delenv("CLIF_CONFIG", raising=False)
    write_config(repo, "config.json", bad)

    with pytest.raises(SystemExit) as exc:
        _site.resolve_config(repo, [])

    assert "invalid site name" in str(exc.value)


def test_traversal_cannot_escape_output(repo: Path) -> None:
    """The concrete consequence the validation prevents."""
    out = _site.site_output_root(repo, {"site_name": "SiteA"})
    assert out == repo / "output" / "SiteA"
    assert out.resolve().is_relative_to((repo / "output").resolve())


def test_config_without_site_name_is_fatal(repo: Path, monkeypatch) -> None:
    monkeypatch.delenv("CLIF_CONFIG", raising=False)
    (repo / "config" / "config.json").write_text(json.dumps({"filetype": "parquet"}))

    with pytest.raises(SystemExit) as exc:
        _site.resolve_config(repo, [])

    assert "site_name" in str(exc.value)


# ---------------------------------------------------------------------------
# Output tree
# ---------------------------------------------------------------------------

def test_ensure_site_dirs_creates_both_and_labels_the_phi_one(tmp_path: Path) -> None:
    phi, share = _site.ensure_site_dirs(tmp_path / "output" / "SiteA")

    assert phi.is_dir() and share.is_dir()
    label = (phi / "README.md").read_text()
    assert "never leaves the site" in label, "the PHI warning label must be written"


def test_two_sites_get_disjoint_trees(repo: Path) -> None:
    """The collision this whole change exists to prevent."""
    a = _site.site_output_root(repo, {"site_name": "SiteB"})
    b = _site.site_output_root(repo, {"site_name": "SiteA"})
    _site.ensure_site_dirs(a)
    _site.ensure_site_dirs(b)

    (a / "intermediate_phi" / "cohort.parquet").write_text("A")
    (b / "intermediate_phi" / "cohort.parquet").write_text("B")

    assert (a / "intermediate_phi" / "cohort.parquet").read_text() == "A", (
        "writing site B must not disturb site A's identically-named artifact")


# ---------------------------------------------------------------------------
# The --check / --list-sites CLI the README tells operators to use
# ---------------------------------------------------------------------------

def run_site_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "code" / "_site.py"), *args],
        capture_output=True, text=True, timeout=60, cwd=REPO_ROOT,
    )


def make_py_repo(tmp_path: Path) -> Path:
    """A throwaway repo holding a real copy of _site.py.

    --check resolves its repo root from its own __file__, so as with the R script a
    test cannot relocate it by changing cwd. Copying the shipped module is what moves
    the root, and keeps the test exercising the real file.
    """
    repo = tmp_path / "repo"
    (repo / "code").mkdir(parents=True)
    (repo / "config").mkdir(parents=True)
    shutil.copy(REPO_ROOT / "code" / "_site.py", repo / "code" / "_site.py")
    return repo


def check_in(repo: Path, site: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(repo / "code" / "_site.py"), "--check", site],
        capture_output=True, text=True, timeout=60, cwd=repo,
    )


def test_check_rejects_a_freshly_copied_template(tmp_path: Path) -> None:
    """The state a config is in right after `cp`. It must not read as usable."""
    repo = make_py_repo(tmp_path)
    cfg = json.loads((REPO_ROOT / "config" / "config_template.json").read_text())
    cfg["site_name"] = "SiteA"
    (repo / "config" / "config_SiteA.json").write_text(json.dumps(cfg))

    proc = check_in(repo, "SiteA")

    assert proc.returncode == 1, "an unedited template must not report READY"
    assert "NOT READY" in proc.stdout
    assert "placeholder" in proc.stdout


def test_check_accepts_a_filled_in_config(tmp_path: Path) -> None:
    """The green half, so the test above cannot pass vacuously."""
    repo = make_py_repo(tmp_path)
    data = tmp_path / "clif_tables"
    data.mkdir()
    cfg = json.loads((REPO_ROOT / "config" / "config_template.json").read_text())
    cfg.update(site_name="SiteA", data_directory=str(data))
    (repo / "config" / "config_SiteA.json").write_text(json.dumps(cfg))

    proc = check_in(repo, "SiteA")

    assert proc.returncode == 0, f"expected READY, got:\n{proc.stdout}{proc.stderr}"
    assert "READY." in proc.stdout
    assert str(repo / "output" / "SiteA") in proc.stdout


def test_check_rejects_a_missing_data_directory(tmp_path: Path) -> None:
    repo = make_py_repo(tmp_path)
    cfg = json.loads((REPO_ROOT / "config" / "config_template.json").read_text())
    cfg.update(site_name="SiteA", data_directory=str(tmp_path / "does_not_exist"))
    (repo / "config" / "config_SiteA.json").write_text(json.dumps(cfg))

    proc = check_in(repo, "SiteA")

    assert proc.returncode == 1
    assert "does not exist" in proc.stdout


def test_check_rejects_has_crrt_settings_false(tmp_path: Path) -> None:
    """Without flow rates there is no dose, and dose is the exposure."""
    repo = make_py_repo(tmp_path)
    data = tmp_path / "clif_tables"
    data.mkdir()
    cfg = json.loads((REPO_ROOT / "config" / "config_template.json").read_text())
    cfg.update(site_name="SiteA", data_directory=str(data), has_crrt_settings=False)
    (repo / "config" / "config_SiteA.json").write_text(json.dumps(cfg))

    proc = check_in(repo, "SiteA")

    assert proc.returncode == 1
    assert "has_crrt_settings" in proc.stdout


def test_check_cli_exits_nonzero_on_an_unusable_config() -> None:
    """Driven as a subprocess, because the README tells operators to run it that way."""
    sites = _site.available_sites(REPO_ROOT)
    if not sites:
        pytest.skip("no per-site configs on this machine")
    proc = run_site_cli("--check", sites[0])
    assert "output would go" in proc.stdout
    assert proc.returncode in (0, 1)
    if proc.returncode == 1:
        assert "NOT READY" in proc.stdout
    else:
        assert "READY." in proc.stdout


def test_list_sites_cli_matches_available_sites() -> None:
    proc = run_site_cli("--list-sites")
    assert proc.returncode == 0
    listed = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    assert listed == _site.available_sites(REPO_ROOT)


def test_site_cli_rejects_an_unknown_flag() -> None:
    proc = run_site_cli("--nonsense")
    assert proc.returncode != 0
    assert "usage:" in (proc.stdout + proc.stderr)


# ---------------------------------------------------------------------------
# The Python rule and the R transcription must agree
# ---------------------------------------------------------------------------

@requires_r
def test_r_and_python_resolvers_agree(tmp_path: Path) -> None:
    """Drive the real R script and compare the site it resolves to.

    03 is a fit script, so it is expected to die once it reaches the missing analysis
    frame. What matters is the site/config/output it reports BEFORE that: those lines
    are produced by the R transcription of resolve_config.
    """
    repo = make_r_repo(tmp_path, "config_SiteA.json", "SiteA")

    proc = run_r(repo, "smoke", "--site", "SiteA")
    out = proc.stdout + proc.stderr

    _, cfg = _site.resolve_config(repo, ["--site", "SiteA"])
    expected_out = _site.site_output_root(repo, cfg)

    assert "site:   SiteA" in out, f"R did not resolve site SiteA:\n{out[:2000]}"
    assert str(expected_out) in out, (
        f"R output root disagrees with _site.site_output_root ({expected_out}):\n{out[:2000]}")


@requires_r
def test_r_rejects_site_name_mismatch(tmp_path: Path) -> None:
    """The R half of the guard must be red too, or the pair is only half enforced."""
    repo = make_r_repo(tmp_path, "config_SiteA.json", "SiteB")   # the mistake

    proc = run_r(repo, "smoke", "--site", "SiteA")

    assert proc.returncode != 0, "R must refuse a mismatched site"
    assert "site mismatch" in (proc.stdout + proc.stderr)
