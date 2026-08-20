# CRRT-dose-lmtp pipeline runner (Windows PowerShell).
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
#
# Usage:  .\run_pipeline.ps1              (config\config.json)
#         .\run_pipeline.ps1 SiteA           (config\config_SiteA.json)
# If execution policy blocks it:
#   powershell -ExecutionPolicy Bypass -File .\run_pipeline.ps1

param([string]$Site = "")

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir
$env:PYTHONIOENCODING = "utf-8"

# Which site. Precedence must match resolve_config() in code\_site.py and the
# transcription in code\03_lmtp_fit.R:
#
#   1. positional arg   .\run_pipeline.ps1 SiteA  ->  config\config_SiteA.json
#   2. $env:CLIF_CONFIG                        ->  that path
#   3. default                                 ->  config\config.json
#
# THE ARGUMENT IS ALSO PASSED TO EVERY STEP. It did not used to be: this script
# resolved CLIF_CONFIG for preflight while 01/02/03 each hardcoded config\config.json,
# so preflight validated one site and execution ran another, exiting 0.
$SiteArgs = @()
$SiteSuffix = ""
if ($Site) {
    $Cfg = Join-Path $ScriptDir "config\config_$Site.json"
    $SiteArgs = @("--site", $Site)
    $SiteSuffix = " --site $Site"
    $CfgHint = "copy config\config_template.json config\config_$Site.json"
} elseif ($env:CLIF_CONFIG) {
    $Cfg = $env:CLIF_CONFIG
    $CfgHint = "point CLIF_CONFIG at a config that exists, or clear it"
} else {
    $Cfg = Join-Path $ScriptDir "config\config.json"
    $CfgHint = "copy config\config_template.json config\config.json"
}

Write-Host "=============================================================="
Write-Host " CRRT-dose-lmtp"
Write-Host "=============================================================="

# --- Preflight ---------------------------------------------------------------
# uv is checked FIRST because the missing-config branch below shells out to it.
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Error "uv is not installed. See https://docs.astral.sh/uv/"
    exit 1
}

if (-not (Test-Path $Cfg)) {
    if ((-not $Site) -and (-not $env:CLIF_CONFIG)) {
        # No default config. What to advise depends on whether this machine holds
        # per-site configs, so ask _site.py rather than guessing here: it is the one
        # place that knows how a per-site config is named.
        & uv run --quiet python (Join-Path $ScriptDir "code\_site.py") --explain-missing ".\run_pipeline.ps1" | Write-Host
    } else {
        Write-Error "config not found at $Cfg`n       Fix: $CfgHint`n       then edit it for your site."
    }
    exit 1
}

$SiteName = uv run --quiet python -c "import json,sys;print(json.load(open(sys.argv[1]))['site_name'])" $Cfg
$HasSettings = uv run --quiet python -c "import json,sys;print(json.load(open(sys.argv[1])).get('has_crrt_settings',False))" $Cfg

# The same equality check _site.py enforces, applied here so the run dies during
# preflight rather than a table-read later.
if ($Site -and ($SiteName -ne $Site)) {
    Write-Error "site mismatch.`n       You asked for '$Site', which resolved to $(Split-Path -Leaf $Cfg),`n       but that file declares site_name = '$SiteName'.`n       Refusing to run: outputs would be filed and stamped under '$SiteName', not '$Site'."
    exit 1
}

Write-Host "Site   : $SiteName"
Write-Host "Config : $Cfg"
Write-Host "Output : output\$SiteName\"
Write-Host ""

if ($HasSettings -ne "True") {
    Write-Error "has_crrt_settings is not true.`n       Dose is the exposure for this study and cannot be computed`n       without CRRT flow rates. A site without them cannot participate."
    exit 1
}

# Lockfile completeness. Cheap, needs no R, and fails BEFORE step 01's long table
# read rather than at step 03's library() call an hour later. nanoparquet shipped
# unrecorded in renv.lock and only the coordinating machine survived it, because the
# package happened to be installed there. See code\check_r_deps.py.
& uv run --quiet python (Join-Path $ScriptDir "code\check_r_deps.py")
if ($LASTEXITCODE -ne 0) {
    Write-Error "R dependency preflight failed. Nothing was run."
    exit 1
}
Write-Host ""

# --- Step inventory ----------------------------------------------------------
# There is no step 00. The plan to vendor 00_cohort.py from
# CLIF-epidemiology-of-CRRT was abandoned on 2026-08-16; code\vendor\ now pins two
# config files and no code. See README, "Pipeline steps".
$PythonSteps = @(
    "code\01_build_cohort.py",
    "code\02_build_lmtp_df.py"
)
$RScript = "code\03_lmtp_fit.R"

$Missing = @()
foreach ($step in ($PythonSteps + $RScript)) {
    if (-not (Test-Path (Join-Path $ScriptDir $step))) { $Missing += $step }
}

if ($Missing.Count -gt 0) {
    Write-Host "--------------------------------------------------------------"
    Write-Host " PIPELINE INCOMPLETE. Nothing was run."
    Write-Host "--------------------------------------------------------------"
    Write-Host "Missing step(s):"
    foreach ($m in $Missing) { Write-Host "  - $m" }
    Write-Host ""
    Write-Host "Preflight itself passed: config is valid and the environment resolves."
    exit 2
}

# --- Execution ---------------------------------------------------------------
$Start = Get-Date

foreach ($step in $PythonSteps) {
    Write-Host ">>> $step"
    $StepStart = Get-Date
    & uv run python (Join-Path $ScriptDir $step) @SiteArgs
    if ($LASTEXITCODE -ne 0) {
        Write-Host ""
        Write-Host "--------------------------------------------------------------"
        Write-Host " FAILED: $step"
        Write-Host "--------------------------------------------------------------"
        if ($step -eq "code\02_build_lmtp_df.py") {
            Write-Host "If this failed on medication unit conversion, that is the KNOWN"
            Write-Host "NEE blocker and the halt is deliberate, not a crash."
            Write-Host "clifpy leaves the RAW value in med_dose_converted when it cannot"
            Write-Host "convert, so an unconverted vasopressor silently inflates the"
            Write-Host "norepinephrine equivalent. Step 02 raises rather than continue."
            Write-Host "Dropping the rows is NOT a safe workaround: it removes a drug"
            Write-Host "from NEE for exactly the patients who received it."
        }
        exit 1
    }
    Write-Host ("    done in {0:N0}s" -f ((Get-Date) - $StepStart).TotalSeconds)
}

if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) {
    Write-Host ""
    Write-Host "WARNING: Rscript not on PATH, so the smoke fit was skipped."
    Write-Host "         The analysis frame was still built."
    Write-Host "         Install R 4.3+ and run: Rscript $RScript smoke$SiteSuffix"
    exit 3
}

Write-Host ">>> $RScript smoke"
$StepStart = Get-Date
& Rscript (Join-Path $ScriptDir $RScript) smoke @SiteArgs
if ($LASTEXITCODE -ne 0) { Write-Error "smoke fit failed"; exit 1 }
Write-Host ("    done in {0:N0}s" -f ((Get-Date) - $StepStart).TotalSeconds)

$Elapsed = (Get-Date) - $Start
Write-Host ""
Write-Host "--------------------------------------------------------------"
Write-Host (" Frame built and smoke fit passed in {0:N0}m {1:N0}s." -f [math]::Floor($Elapsed.TotalMinutes), $Elapsed.Seconds)
Write-Host "--------------------------------------------------------------"
Write-Host "The remaining stages are run by hand, on purpose:"
Write-Host ""
Write-Host "  Rscript $RScript gate$SiteSuffix     # diagnostics ONLY, no effect estimate"
Write-Host "  # ...read the diagnostics, then:"
Write-Host "  Rscript $RScript expand$SiteSuffix   # the full delta ladder"
Write-Host ""
Write-Host "Shareable outputs: output\$SiteName\final_no_phi\  (PHI-check before sending)"
