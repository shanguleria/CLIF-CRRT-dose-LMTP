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
# Usage:  .\run_pipeline.ps1
# If execution policy blocks it:
#   powershell -ExecutionPolicy Bypass -File .\run_pipeline.ps1

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ScriptDir
$env:PYTHONIOENCODING = "utf-8"

$Cfg = if ($env:CLIF_CONFIG) { $env:CLIF_CONFIG } else { Join-Path $ScriptDir "config\config.json" }

Write-Host "=============================================================="
Write-Host " CRRT-dose-lmtp"
Write-Host "=============================================================="

# --- Preflight ---------------------------------------------------------------
if (-not (Test-Path $Cfg)) {
    Write-Error "config not found at $Cfg`n       Fix: copy config\config_template.json config\config.json`n       then edit it for your site."
    exit 1
}

if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Error "uv is not installed. See https://docs.astral.sh/uv/"
    exit 1
}

$Site = uv run --quiet python -c "import json,sys;print(json.load(open(sys.argv[1]))['site_name'])" $Cfg
$HasSettings = uv run --quiet python -c "import json,sys;print(json.load(open(sys.argv[1])).get('has_crrt_settings',False))" $Cfg

Write-Host "Site   : $Site"
Write-Host "Config : $Cfg"
Write-Host ""

if ($HasSettings -ne "True") {
    Write-Error "has_crrt_settings is not true.`n       Dose is the exposure for this study and cannot be computed`n       without CRRT flow rates. A site without them cannot participate."
    exit 1
}

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
    & uv run python (Join-Path $ScriptDir $step)
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
    Write-Host "         Install R 4.3+ and run: Rscript $RScript smoke"
    exit 3
}

Write-Host ">>> $RScript smoke"
$StepStart = Get-Date
& Rscript (Join-Path $ScriptDir $RScript) smoke
if ($LASTEXITCODE -ne 0) { Write-Error "smoke fit failed"; exit 1 }
Write-Host ("    done in {0:N0}s" -f ((Get-Date) - $StepStart).TotalSeconds)

$Elapsed = (Get-Date) - $Start
Write-Host ""
Write-Host "--------------------------------------------------------------"
Write-Host (" Frame built and smoke fit passed in {0:N0}m {1:N0}s." -f [math]::Floor($Elapsed.TotalMinutes), $Elapsed.Seconds)
Write-Host "--------------------------------------------------------------"
Write-Host "The remaining stages are run by hand, on purpose:"
Write-Host ""
Write-Host "  Rscript $RScript gate     # diagnostics ONLY, no effect estimate"
Write-Host "  # ...read the diagnostics, then:"
Write-Host "  Rscript $RScript expand   # the full delta ladder"
Write-Host ""
Write-Host "Shareable outputs: output\final_no_phi\  (PHI-check before sending)"
