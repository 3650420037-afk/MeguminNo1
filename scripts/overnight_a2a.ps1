$ErrorActionPreference = 'Continue'
$repo = "D:\MMModel\Pocket2Mol"
$pdb = (Get-ChildItem "D:\MMModel\靶点结构" -Filter "4EIY*.pdb" | Select-Object -First 1).FullName
$end = (Get-Date).AddHours(12)
$i = 50
$total = 564
while ((Get-Date) -lt $end) {
  if ($total -ge 4200) { break }
  $free = (Get-PSDrive D).Free / 1GB
  if ($free -lt 15) {
    Get-ChildItem "$repo\outputs" -Recurse -File -Filter "samples_*.pt" -ErrorAction SilentlyContinue | Where-Object { $_.Name -ne 'samples_all.pt' } | Remove-Item -Force -ErrorAction SilentlyContinue
    $free = (Get-PSDrive D).Free / 1GB
    if ($free -lt 8) { Add-Content "$repo\outputs\overnight_a2a\progress.txt" -Encoding UTF8 -Value ("{0} DISK_LOW {1}GB paused" -f (Get-Date -Format 'HH:mm'), [math]::Round($free,1)) -ErrorAction SilentlyContinue; Start-Sleep -Seconds 600; continue }
  }
  $i++
  $seed = 7000 + $i
  $src = [IO.File]::ReadAllText("$repo\configs\sample_for_pdb_guided_l3.yml").TrimStart([char]0xFEFF)
  $cfg = ($src -replace 'seed: \d+', "seed: $seed") + "    diversity_w: 0.5`n"
  [IO.File]::WriteAllText("$repo\configs\overnight_batch.yml", $cfg, (New-Object System.Text.UTF8Encoding($false)))
  $t0 = Get-Date
  Set-Location $repo
  & "D:\Miniconda3\envs\Pocket2Mol\python.exe" sample_for_pdb.py --pdb_path $pdb --center " -0.4,8.5,17.1" --config ".\configs\overnight_batch.yml" --outdir ".\outputs\overnight_a2a\batch_$i" *>&1 | Out-File "$env:TEMP\overnight_b$i.log" -Encoding UTF8
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
  Add-Content "$repo\outputs\overnight_a2a\progress.txt" -Encoding UTF8 -Value ("{0} batch_{1} seed={2} finished={3} total={4} elapsed={5}min" -f (Get-Date -Format 'yyyy-MM-dd HH:mm'), $i, $seed, $n, $total, [math]::Round(((Get-Date)-$t0).TotalMinutes,1)) -ErrorAction SilentlyContinue
}
"OVERNIGHT_DONE batches=$i total=$total $(Get-Date -Format HH:mm)"