$ErrorActionPreference = 'Continue'
# src/scripts/night_phase.ps1 -> 上两级 = 仓库根 (与 paths.py 的 ROOT 一致)
$Root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$OutRoot = Join-Path $Root 'outputs'
$LibRoot = Join-Path $Root 'results'
$TargetsDir = Join-Path $Root 'data\targets'
# 解释器: 优先环境变量 EONMOL_PYTHON (与 paths.PYTHON 同源), 否则用 PATH 上的 python
$py = if ($env:EONMOL_PYTHON) { $env:EONMOL_PYTHON } else { 'python' }
Set-Location $Root
# 等待 Top200 A1 对接完成
while (-not (Test-Path "$OutRoot\docking_a1_top200\summary.csv")) { Start-Sleep -Seconds 120 }
Add-Content "$OutRoot\night_progress.txt" -Encoding UTF8 -Value ("{0} TOP200_DONE" -f (Get-Date -Format 'HH:mm'))
# 导出剩余 2731 个 (非 Top200)
$rows = [System.IO.File]::ReadAllLines("$OutRoot\docking_a1_top200\summary.csv") | Select-Object -Skip 1 | ForEach-Object { ($_ -split ',')[0] }
$top200 = @{}
foreach ($r in $rows) { $top200[$r] = $true }
$all = [System.IO.File]::ReadAllLines("$LibRoot\all_smiles.txt") | Where-Object { $_ -and -not $top200.ContainsKey($_) }
[System.IO.File]::WriteAllLines("$LibRoot\rest2731_smiles.txt", $all)
Add-Content "$OutRoot\night_progress.txt" -Encoding UTF8 -Value ("{0} REST_COUNT={1}" -f (Get-Date -Format 'HH:mm'), $all.Count)
# 全库 A1 对接 (剩余部分)
& $py (Join-Path $Root 'src\scripts\docking_pipeline.py') --receptor (Join-Path $TargetsDir '5UEN.pdb') --center="55.969,58.897,143.624" --smiles-file "$LibRoot\rest2731_smiles.txt" --out "$OutRoot\docking_a1_rest" *>&1 | Out-File "$env:TEMP\dock_a1_rest.log" -Encoding UTF8
Add-Content "$OutRoot\night_progress.txt" -Encoding UTF8 -Value ("{0} REST_EXIT={1}" -f (Get-Date -Format 'HH:mm'), $LASTEXITCODE)
"NIGHT_PHASE_DONE"
