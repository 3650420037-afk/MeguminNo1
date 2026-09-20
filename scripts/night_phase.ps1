$ErrorActionPreference = 'Continue'
$repo = "D:\MMModel\Pocket2Mol"
$py = "D:\Miniconda3\envs\Pocket2Mol\python.exe"
Set-Location $repo
# 等待 Top200 A1 对接完成
while (-not (Test-Path "$repo\outputs\docking_a1_top200\summary.csv")) { Start-Sleep -Seconds 120 }
Add-Content "$repo\outputs\night_progress.txt" -Encoding UTF8 -Value ("{0} TOP200_DONE" -f (Get-Date -Format 'HH:mm'))
# 导出剩余 2731 个 (非 Top200)
$rows = [System.IO.File]::ReadAllLines("$repo\outputs\docking_a1_top200\summary.csv") | Select-Object -Skip 1 | ForEach-Object { ($_ -split ',')[0] }
$top200 = @{}
foreach ($r in $rows) { $top200[$r] = $true }
$all = [System.IO.File]::ReadAllLines("D:\MMModel\化合物库\all_smiles.txt") | Where-Object { $_ -and -not $top200.ContainsKey($_) }
[System.IO.File]::WriteAllLines("D:\MMModel\化合物库\rest2731_smiles.txt", $all)
Add-Content "$repo\outputs\night_progress.txt" -Encoding UTF8 -Value ("{0} REST_COUNT={1}" -f (Get-Date -Format 'HH:mm'), $all.Count)
# 全库 A1 对接 (剩余部分)
& $py scripts\docking_pipeline.py --receptor "D:\MMModel\靶点结构\5UEN.pdb" --center="55.969,58.897,143.624" --smiles-file "D:\MMModel\化合物库\rest2731_smiles.txt" --out ".\outputs\docking_a1_rest" *>&1 | Out-File "$env:TEMP\dock_a1_rest.log" -Encoding UTF8
Add-Content "$repo\outputs\night_progress.txt" -Encoding UTF8 -Value ("{0} REST_EXIT={1}" -f (Get-Date -Format 'HH:mm'), $LASTEXITCODE)
"NIGHT_PHASE_DONE"