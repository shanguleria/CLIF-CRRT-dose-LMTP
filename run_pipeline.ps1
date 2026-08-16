# CRRT-dose-lmtp pipeline runner (Windows PowerShell).
#
# SCAFFOLD. No analysis step exists yet, so this script performs the preflight
# checks and then reports exactly what is missing rather than pretending to run.
# It exits non-zero while the pipeline is incomplete: a runner that exits 0
# without doing anything is how a site ends up believing it has results.
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
$PythonSteps = @(
    "code\vendor\00_cohort.py",
    "code\vendor\01_create_wide_df.py",
    "code\02_build_lmtp_df.py"
)
$RSteps = @(
    "code\03_lmtp_fit.R"
)

$Missing = @()
foreach ($step in ($PythonSteps + $RSteps)) {
    if (-not (Test-Path (Join-Path $ScriptDir $step))) { $Missing += $step }
}

if ($Missing.Count -gt 0) {
    Write-Host "--------------------------------------------------------------"
    Write-Host " PIPELINE INCOMPLETE. Nothing was run."
    Write-Host "--------------------------------------------------------------"
    Write-Host "Missing step(s):"
    foreach ($m in $Missing) { Write-Host "  - $m" }
    Write-Host ""
    Write-Host "Steps 00 and 01 are vendored, not written here:"
    Write-Host "  bash code/vendor/sync_vendor.sh"
    Write-Host "Steps 02 and 03 are not yet implemented. See .claude\claude-todo.md."
    Write-Host ""
    Write-Host "Preflight itself passed: config is valid and the environment resolves."
    exit 2
}

# --- Execution (reached only once every step exists) -------------------------
$Total = [System.Diagnostics.Stopwatch]::StartNew()
foreach ($step in $PythonSteps) {
    Write-Host ">>> $step"
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    uv run python (Join-Path $ScriptDir $step)
    if ($LASTEXITCODE -ne 0) { Write-Error "step failed: $step"; exit 1 }
    Write-Host ("    done in {0:N0}s" -f $sw.Elapsed.TotalSeconds)
}

if (Get-Command Rscript -ErrorAction SilentlyContinue) {
    foreach ($step in $RSteps) {
        Write-Host ">>> $step"
        $sw = [System.Diagnostics.Stopwatch]::StartNew()
        Rscript (Join-Path $ScriptDir $step)
        if ($LASTEXITCODE -ne 0) { Write-Error "step failed: $step"; exit 1 }
        Write-Host ("    done in {0:N0}s" -f $sw.Elapsed.TotalSeconds)
    }
} else {
    Write-Warning "Rscript not on PATH; R steps skipped."
}

Write-Host ("Total: {0:N0}m {1:N0}s" -f [math]::Floor($Total.Elapsed.TotalMinutes), $Total.Elapsed.Seconds)
Write-Host "Shareable outputs: output\final_no_phi\  (PHI-check before sending)"
