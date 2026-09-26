$ErrorActionPreference = 'Continue'
# src/scripts/dock_3targets.ps1 -> 上两级 = 仓库根 (与 paths.py 的 ROOT 一致)
$Root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$TargetsDir = Join-Path $Root 'data\targets'
$LibRoot = Join-Path $Root 'results'
# 解释器: 优先环境变量 EONMOL_PYTHON (与 paths.PYTHON 同源), 否则用 PATH 上的 python
$py = if ($env:EONMOL_PYTHON) { $env:EONMOL_PYTHON } else { 'python' }
Set-Location $Root
$targets = @(
  @{ key='b2ar';  pdb=(Get-ChildItem $TargetsDir -Filter "2RH1*.pdb" | Select-Object -First 1).FullName; center='-29.5,9.2,6.9' },
  @{ key='d3';    pdb=(Get-ChildItem $TargetsDir -Filter "3PBL*.pdb" | Select-Object -First 1).FullName; center='0.085,-14.828,10.432' },
  @{ key='5ht2b'; pdb=(Get-ChildItem $TargetsDir -Filter "4IB4*.pdb"  | Select-Object -First 1).FullName; center='22.448,18.284,11.726' }
)
foreach ($t in $targets) {
  $lib = Join-Path $LibRoot "$($t.key)\compounds.csv"
  if (-not (Test-Path $lib)) { continue }
  $smiFile = Join-Path $LibRoot "$($t.key)\all_smiles.txt"
  Import-Csv $lib | Select-Object -ExpandProperty smiles | Set-Content $smiFile -Encoding UTF8
  & $py (Join-Path $Root 'src\scripts\docking_pipeline.py') --receptor "$($t.pdb)" --center="$($t.center)" --smiles-file $smiFile --out (Join-Path $Root "outputs\docking_$($t.key)_library") *>&1 | Out-File "$env:TEMP\dock_$($t.key).log" -Encoding UTF8
  Add-Content (Join-Path $Root 'outputs\dock_progress.txt') -Encoding UTF8 -Value ("{0} {1} EXIT={2}" -f (Get-Date -Format 'HH:mm'), $t.key, $LASTEXITCODE)
}
"ALL_DOCKING_DONE"
