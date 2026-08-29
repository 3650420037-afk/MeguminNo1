$ErrorActionPreference = 'Continue'
$repo = "D:\MMModel\Pocket2Mol"
$py = "D:\Miniconda3\envs\Pocket2Mol\python.exe"
$out = "$repo\outputs\overnight_phase2"
New-Item $out -ItemType Directory -Force | Out-Null
Set-Location $repo

# ===== 阶段1: 四靶点正式微调 2000 步 =====
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE1_TRAIN_START" -f (Get-Date -Format 'HH:mm'))
& $py train.py --config ".\configs\train_multitarget_v2.yml" --logdir ".\logs" *>&1 | Out-File "$env:TEMP\phase2_train.log" -Encoding UTF8
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE1_TRAIN_EXIT={1}" -f (Get-Date -Format 'HH:mm'), $LASTEXITCODE)
$ckptDir = Get-ChildItem "$repo\logs" -Directory | Where-Object Name -like 'train_multitarget_v2*' | Sort-Object Name -Descending | Select-Object -First 1
$ckpt = Join-Path $ckptDir.FullName "checkpoints\2000.pt"
if (-not (Test-Path $ckpt)) {
  $all = Get-ChildItem (Join-Path $ckptDir.FullName "checkpoints") -Filter "*.pt" -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
  $ckpt = if ($all) { $all.FullName } else { "" }
}
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("CKPT={0}" -f $ckpt)

# ===== 阶段2: A/B 快评 (A2A, 30/50/20, seed 2024) =====
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE2_AB_START" -f (Get-Date -Format 'HH:mm'))
$pdb = (Get-ChildItem "D:\MMModel\靶点结构" -Filter "4EIY*.pdb" | Select-Object -First 1).FullName
& $py sample_for_pdb.py --pdb_path $pdb --center " -0.4,8.5,17.1" --config ".\configs\sample_for_pdb_ab_baseline_305020.yml" --outdir "$out\abft_baseline" *>&1 | Out-File "$env:TEMP\abft_b.log" -Encoding UTF8
$src = [IO.File]::ReadAllText("$repo\configs\sample_for_pdb_ab_local_305020.yml").TrimStart([char]0xFEFF)
$ft = ($src -replace 'checkpoint: [^\r\n]+', ("checkpoint: " + ($ckpt -replace '\\','/'))) -replace 'seed: \d+', 'seed: 2024'
[IO.File]::WriteAllText("$repo\configs\abft_finetuned.yml", $ft, (New-Object System.Text.UTF8Encoding($false)))
& $py sample_for_pdb.py --pdb_path $pdb --center " -0.4,8.5,17.1" --config ".\configs\abft_finetuned.yml" --outdir "$out\abft_finetuned" *>&1 | Out-File "$env:TEMP\abft_f.log" -Encoding UTF8
& $py scripts\judge_promotion.py "$out\abft_baseline" "$out\abft_finetuned" "$out\decision.txt"
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE2_AB_DONE decision={1}" -f (Get-Date -Format 'HH:mm'), ((Get-Content "$out\decision.txt" | Select-String 'DECISION').Line))

# ===== 阶段3: 三靶点轮转生成到 07:30 =====
$decision = Get-Content "$out\decision.txt" -Encoding UTF8 | Select-String 'DECISION' | ForEach-Object { $_.Line }
$useFt = $decision -match 'PROMOTED'
$ckptTemplate = [IO.File]::ReadAllText("$repo\configs\sample_for_pdb_guided_l3.yml").TrimStart([char]0xFEFF)
if ($useFt -and $ckpt) { $ckptTemplate = $ckptTemplate -replace 'checkpoint: [^\r\n]+', ("checkpoint: " + ($ckpt -replace '\\','/')) }
$targets = @(
  @{ key='b2ar';  pdb=(Get-ChildItem "D:\MMModel\靶点结构" -Filter "2RH1*.pdb" | Select-Object -First 1).FullName; center=' -29.5,9.2,6.9';      seedBase=8100 },
  @{ key='d3';    pdb=(Get-ChildItem "D:\MMModel\靶点结构" -Filter "3PBL*.pdb" | Select-Object -First 1).FullName; center=' 0.085,-14.828,10.432'; seedBase=8200 },
  @{ key='5ht2b'; pdb=(Get-ChildItem "D:\MMModel\靶点结构" -Filter "4IB4*.pdb"  | Select-Object -First 1).FullName; center=' 22.448,18.284,11.726'; seedBase=8300 }
)
$end = Get-Date '07:30'
if ((Get-Date) -gt $end) { $end = (Get-Date).AddHours(11) }
$counts = @{ b2ar=0; d3=0; '5ht2b'=0 }
$totals = @{ b2ar=0; d3=0; '5ht2b'=0 }
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE3_GEN_START ckpt={(1)}" -f (Get-Date -Format 'HH:mm'), $(if ($useFt) { 'finetuned' } else { 'pretrained' }))
while ((Get-Date) -lt $end) {
  $free = (Get-PSDrive D).Free / 1GB
  if ($free -lt 12) {
    Get-ChildItem "$repo\outputs" -Recurse -File -Filter "samples_*.pt" -ErrorAction SilentlyContinue | Where-Object { $_.Name -ne 'samples_all.pt' } | Remove-Item -Force -ErrorAction SilentlyContinue
    $free = (Get-PSDrive D).Free / 1GB
    if ($free -lt 6) { Start-Sleep -Seconds 900; continue }
  }
  foreach ($t in $targets) {
    if ((Get-Date) -ge $end) { break }
    $key = $t.key; $counts[$key]++
    $seed = $t.seedBase + $counts[$key]
    $cfg = ($ckptTemplate -replace 'seed: \d+', "seed: $seed") + "    diversity_w: 0.5`n"
    [IO.File]::WriteAllText("$repo\configs\gen_batch.yml", $cfg, (New-Object System.Text.UTF8Encoding($false)))
    $t0 = Get-Date
    & $py sample_for_pdb.py --pdb_path $t.pdb --center $t.center --config ".\configs\gen_batch.yml" --outdir "$out\gen_$key\batch_$($counts[$key])" *>&1 | Out-File "$env:TEMP\gen_${key}_$($counts[$key]).log" -Encoding UTF8
    $d = Get-ChildItem "$out\gen_$key\batch_$($counts[$key])" -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
    $n = 0
    if ($d) {
      Get-ChildItem $d.FullName -File -Filter "samples_*.pt" -ErrorAction SilentlyContinue | Where-Object { $_.Name -ne 'samples_all.pt' } | Remove-Item -Force -ErrorAction SilentlyContinue
      $smi = Join-Path $d.FullName "SMILES.txt"
      if (Test-Path $smi) {
        $lines = @(Get-Content $smi -Encoding UTF8 | Where-Object { $_ })
        $n = $lines.Count
        if ($n -gt 0) { Add-Content "$out\gen_$key\all_smiles.txt" -Value $lines -Encoding UTF8 -ErrorAction SilentlyContinue }
      }
    }
    $totals[$key] += $n
    Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} {1} batch_{2} seed={3} finished={4} total={5}" -f (Get-Date -Format 'HH:mm'), $key, $counts[$key], $seed, $n, $totals[$key]) -ErrorAction SilentlyContinue
  }
}
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} OVERNIGHT_PHASE2_DONE b2ar={1} d3={2} 5ht2b={3}" -f (Get-Date -Format 'HH:mm'), $totals['b2ar'], $totals['d3'], $totals['5ht2b'])