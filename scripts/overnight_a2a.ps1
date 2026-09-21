$ErrorActionPreference = 'Continue'
$repo = "D:\MMModel\Pocket2Mol"
$py = "D:\Miniconda3\envs\Pocket2Mol\python.exe"
$gen = "$repo\scripts\gen_sample_config.py"
$tmpl = "$repo\configs\sample_for_pdb_guided_l3.yml"
$pdb = (Get-ChildItem "D:\MMModel\靶点结构" -Filter "4EIY*.pdb" | Select-Object -First 1).FullName
$end = (Get-Date).AddHours(12)
$i = 0
$total = 0
while ((Get-Date) -lt $end) {
  $i++
  $seed = 7000 + $i
  $cfg = "$repo\configs\overnight_batch.yml"
  & $py $gen --template $tmpl --out $cfg --seed $seed --diversity-w 0.5 --guided 1 | Out-Null
  if ($LASTEXITCODE -ne 0 -or -not (Test-Path $cfg)) {
    Add-Content "$repo\outputs\overnight_a2a\progress.txt" -Encoding UTF8 -Value ("{0} CONFIG_GEN_FAIL batch={1} exit={2} -> 中止" -f (Get-Date -Format 'HH:mm'), $i, $LASTEXITCODE) -ErrorAction SilentlyContinue
    break
  }
  $t0 = Get-Date
  Set-Location $repo
  & $py sample_for_pdb.py --pdb_path $pdb --center " -0.4,8.5,17.1" --config $cfg --outdir ".\outputs\overnight_a2a\batch_$i" *>&1 | Out-File "$env:TEMP\overnight_b$i.log" -Encoding UTF8
  $rc = $LASTEXITCODE
  $d = Get-ChildItem "$repo\outputs\overnight_a2a\batch_$i" -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
  $n = 0
  if ($d) {
    Get-ChildItem $d.FullName -File -Filter "samples_*.pt" -ErrorAction SilentlyContinue | Where-Object { $_.Name -ne 'samples_all.pt' } | Remove-Item -Force -ErrorAction SilentlyContinue
    $smi = Join-Path $d.FullName "SMILES.txt"
    if (Test-Path $smi) {
      $lines = @(Get-Content $smi -Encoding UTF8 | Where-Object { $_ })
      $n = $lines.Count
      if ($n -gt 0) { Add-Content "$repo\outputs\overnight_a2a\all_smiles.txt" -Value $lines -Encoding UTF8 -ErrorAction SilentlyContinue }
    }
  }
  $total += $n
  Add-Content "$repo\outputs\overnight_a2a\progress.txt" -Encoding UTF8 -Value ("{0} batch_{1} seed={2} finished={3} total={4} elapsed={5}min exit={6}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm'), $i, $seed, $n, $total, [math]::Round(((Get-Date)-$t0).TotalMinutes,1), $rc) -ErrorAction SilentlyContinue
}
"OVERNIGHT_DONE batches=$i total=$total"