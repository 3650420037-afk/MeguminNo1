$ErrorActionPreference = 'Continue'
$repo = "D:\MMModel\Pocket2Mol"
$py = "D:\Miniconda3\envs\Pocket2Mol\python.exe"
Set-Location $repo
$targets = @(
  @{ key='b2ar';  pdb=(Get-ChildItem "D:\MMModel\靶点结构" -Filter "2RH1*.pdb" | Select-Object -First 1).FullName; center='-29.5,9.2,6.9' },
  @{ key='d3';    pdb=(Get-ChildItem "D:\MMModel\靶点结构" -Filter "3PBL*.pdb" | Select-Object -First 1).FullName; center='0.085,-14.828,10.432' },
  @{ key='5ht2b'; pdb=(Get-ChildItem "D:\MMModel\靶点结构" -Filter "4IB4*.pdb"  | Select-Object -First 1).FullName; center='22.448,18.284,11.726' }
)
foreach ($t in $targets) {
  $lib = "D:\MMModel\化合物库\$($t.key)\compounds.csv"
  if (-not (Test-Path $lib)) { continue }
  $smiFile = "D:\MMModel\化合物库\$($t.key)\all_smiles.txt"
  Import-Csv $lib | Select-Object -ExpandProperty smiles | Set-Content $smiFile -Encoding UTF8
  & $py scripts\docking_pipeline.py --receptor "$($t.pdb)" --center="$($t.center)" --smiles-file $smiFile --out ".\outputs\docking_$($t.key)_library" *>&1 | Out-File "$env:TEMP\dock_$($t.key).log" -Encoding UTF8
  Add-Content "$repo\outputs\dock_progress.txt" -Encoding UTF8 -Value ("{0} {1} EXIT={2}" -f (Get-Date -Format 'HH:mm'), $t.key, $LASTEXITCODE)
}
"ALL_DOCKING_DONE"