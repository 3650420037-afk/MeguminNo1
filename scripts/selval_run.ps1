$ErrorActionPreference = 'Continue'
$repo = "D:\MMModel\Pocket2Mol"
$py = "D:\Miniconda3\envs\Pocket2Mol\python.exe"
Set-Location $repo
$a2a = (Get-ChildItem "D:\MMModel\靶点结构" -Filter "4EIY*.pdb" | Select-Object -First 1).FullName
$a1  = "D:\MMModel\靶点结构\5UEN.pdb"
$jobs = @(
  @{g='sel_a2a'; r=$a2a; c='-0.4,8.5,17.1'},
  @{g='sel_a2a'; r=$a1;  c='55.969,58.897,143.624'},
  @{g='sel_a1';  r=$a2a; c='-0.4,8.5,17.1'},
  @{g='sel_a1';  r=$a1;  c='55.969,58.897,143.624'}
)
foreach ($j in $jobs) {
  $smifile = "D:\MMModel\已知药物库\selectivity_a2a_selective.txt"
  if ($j.g -eq 'sel_a1') { $smifile = "D:\MMModel\已知药物库\selectivity_a1_selective.txt" }
  & $py scripts\docking_pipeline.py --receptor $j.r --center="$($j.c)" --smiles-file $smifile --out ".\outputs\selval_$($j.g)_$(Split-Path $j.r -Leaf)" *>&1 | Out-File "$env:TEMP\selval_$($j.g).log" -Encoding UTF8
  Add-Content "$repo\outputs\selval_progress.txt" -Encoding UTF8 -Value ("{0} {1} {2} EXIT={3}" -f (Get-Date -Format 'HH:mm'), $j.g, (Split-Path $j.r -Leaf), $LASTEXITCODE)
}
"SELVAL_DONE"