$ErrorActionPreference = 'Continue'
# src/scripts/overnight_a2a.ps1 -> 上两级 = 仓库根 (与 paths.py 的 ROOT 一致)
$Root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$TargetsDir = Join-Path $Root 'data\targets'
$OutRoot = Join-Path $Root 'outputs'
$ConfigsDir = Join-Path $Root 'configs'
$SrcDir = Join-Path $Root 'src'
# 解释器: 优先环境变量 EONMOL_PYTHON (与 paths.PYTHON 同源), 否则用 PATH 上的 python
$py = if ($env:EONMOL_PYTHON) { $env:EONMOL_PYTHON } else { 'python' }
$gen = Join-Path $SrcDir 'scripts\gen_sample_config.py'
$tmpl = Join-Path $ConfigsDir 'sample_for_pdb_guided_l3.yml'
$samplePy = Join-Path $SrcDir 'sample_for_pdb.py'
$pdb = (Get-ChildItem $TargetsDir -Filter "4EIY*.pdb" | Select-Object -First 1).FullName
$end = (Get-Date).AddHours(12)
$i = 0
$total = 0
while ((Get-Date) -lt $end) {
  $i++
  $seed = 7000 + $i
  $cfg = Join-Path $ConfigsDir 'overnight_batch.yml'
  & $py $gen --template $tmpl --out $cfg --seed $seed --diversity-w 0.5 --guided 1 | Out-Null
  if ($LASTEXITCODE -ne 0 -or -not (Test-Path $cfg)) {
    Add-Content "$OutRoot\overnight_a2a\progress.txt" -Encoding UTF8 -Value ("{0} CONFIG_GEN_FAIL batch={1} exit={2} -> 中止" -f (Get-Date -Format 'HH:mm'), $i, $LASTEXITCODE) -ErrorAction SilentlyContinue
    break
  }
  $t0 = Get-Date
  Set-Location $Root
  & $py $samplePy --pdb_path $pdb --center " -0.4,8.5,17.1" --config $cfg --outdir "$OutRoot\overnight_a2a\batch_$i" *>&1 | Out-File "$env:TEMP\overnight_b$i.log" -Encoding UTF8
  $rc = $LASTEXITCODE
  $d = Get-ChildItem "$OutRoot\overnight_a2a\batch_$i" -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
  $n = 0
  if ($d) {
    Get-ChildItem $d.FullName -File -Filter "samples_*.pt" -ErrorAction SilentlyContinue | Where-Object { $_.Name -ne 'samples_all.pt' } | Remove-Item -Force -ErrorAction SilentlyContinue
    $smi = Join-Path $d.FullName "SMILES.txt"
    if (Test-Path $smi) {
      $lines = @(Get-Content $smi -Encoding UTF8 | Where-Object { $_ })
      $n = $lines.Count
      if ($n -gt 0) { Add-Content "$OutRoot\overnight_a2a\all_smiles.txt" -Value $lines -Encoding UTF8 -ErrorAction SilentlyContinue }
    }
  }
  $total += $n
  Add-Content "$OutRoot\overnight_a2a\progress.txt" -Encoding UTF8 -Value ("{0} batch_{1} seed={2} finished={3} total={4} elapsed={5}min exit={6}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm'), $i, $seed, $n, $total, [math]::Round(((Get-Date)-$t0).TotalMinutes,1), $rc) -ErrorAction SilentlyContinue
}
"OVERNIGHT_DONE batches=$i total=$total"
