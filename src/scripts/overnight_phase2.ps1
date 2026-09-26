$ErrorActionPreference = 'Continue'
# src/scripts/overnight_phase2.ps1 -> 上两级 = 仓库根 (与 paths.py 的 ROOT 一致)
$Root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$TargetsDir = Join-Path $Root 'data\targets'
$ConfigsDir = Join-Path $Root 'configs'
$SrcDir = Join-Path $Root 'src'
# 解释器: 优先环境变量 EONMOL_PYTHON (与 paths.PYTHON 同源), 否则用 PATH 上的 python
$py = if ($env:EONMOL_PYTHON) { $env:EONMOL_PYTHON } else { 'python' }
$gen = Join-Path $SrcDir 'scripts\gen_sample_config.py'
$tmpl = Join-Path $ConfigsDir 'sample_for_pdb_guided_l3.yml'
$samplePy = Join-Path $SrcDir 'sample_for_pdb.py'
$out = Join-Path $Root 'outputs\overnight_phase2'
New-Item $out -ItemType Directory -Force | Out-Null
Set-Location $Root

# ===== 阶段1: 四靶点正式微调 =====
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE1_TRAIN_START" -f (Get-Date -Format 'HH:mm'))
& $py (Join-Path $Root 'train.py') --config ".\configs\train_multitarget_v2.yml" --logdir ".\logs" *>&1 | Out-File "$env:TEMP\phase2_train.log" -Encoding UTF8
$trainRc = $LASTEXITCODE
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE1_TRAIN_EXIT={1}" -f (Get-Date -Format 'HH:mm'), $trainRc)
if ($trainRc -ne 0) { Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value "PHASE1_FAIL -> 中止"; "PHASE1_FAIL"; exit 1 }
$ckptDir = Get-ChildItem (Join-Path $Root 'logs') -Directory -ErrorAction SilentlyContinue | Where-Object Name -like 'train_multitarget_v2*' | Sort-Object Name -Descending | Select-Object -First 1
$ckpt = ""
if ($ckptDir) {
  # 数值排序 (旧写法 Sort-Object Name 是字典序: 800.pt 会胜过 2000.pt)
  $cf = Get-ChildItem (Join-Path $ckptDir.FullName "checkpoints") -Filter "*.pt" -ErrorAction SilentlyContinue |
        Where-Object { $_.BaseName -match '^\d+$' } | Sort-Object { [int]$_.BaseName } -Descending | Select-Object -First 1
  if ($cf) { $ckpt = $cf.FullName }
}
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("CKPT={0}" -f $ckpt)
if (-not $ckpt) { Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value "PHASE1_NO_CKPT -> 中止"; "NO_CKPT"; exit 1 }

# ===== 阶段2: A/B 快评 =====
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE2_AB_START" -f (Get-Date -Format 'HH:mm'))
$pdb = (Get-ChildItem $TargetsDir -Filter "4EIY*.pdb" | Select-Object -First 1).FullName
& $py $samplePy --pdb_path $pdb --center " -0.4,8.5,17.1" --config ".\configs\sample_for_pdb_ab_baseline_305020.yml" --outdir "$out\abft_baseline" *>&1 | Out-File "$env:TEMP\abft_b.log" -Encoding UTF8
& $py $gen --template ".\configs\sample_for_pdb_ab_local_305020.yml" --out ".\configs\abft_finetuned.yml" --seed 2024 | Out-Null
$ft = [IO.File]::ReadAllText((Join-Path $ConfigsDir 'abft_finetuned.yml')).TrimStart([char]0xFEFF)
$ft = $ft -replace 'checkpoint: [^\r\n]+', ("checkpoint: " + ($ckpt -replace '\\','/'))
[IO.File]::WriteAllText((Join-Path $ConfigsDir 'abft_finetuned.yml'), $ft, (New-Object System.Text.UTF8Encoding($false)))
& $py $samplePy --pdb_path $pdb --center " -0.4,8.5,17.1" --config ".\configs\abft_finetuned.yml" --outdir "$out\abft_finetuned" *>&1 | Out-File "$env:TEMP\abft_f.log" -Encoding UTF8
& $py (Join-Path $SrcDir 'scripts\judge_promotion.py') "$out\abft_baseline" "$out\abft_finetuned" "$out\decision.txt"
$decision = (Get-Content "$out\decision.txt" -Encoding UTF8 | Select-String 'DECISION' | ForEach-Object { $_.Line })
if (-not $decision) { Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value "PHASE2_NO_DECISION -> 中止"; "NO_DECISION"; exit 1 }
# 注意: 'NOT_PROMOTED' 含子串 'PROMOTED', 必须用锚定正则, 否则闸门失效
$useFt = ($decision -match 'DECISION:\s*PROMOTED\s*$')
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE2_AB_DONE decision={1} useFt={2}" -f (Get-Date -Format 'HH:mm'), $decision, $useFt)

# ===== 阶段3: 三靶点轮转生成 (到次日 07:30) =====
$targets = @(
  @{ key='b2ar';  pdb=(Get-ChildItem $TargetsDir -Filter "2RH1*.pdb" | Select-Object -First 1).FullName; center=' -29.5,9.2,6.9';      seedBase=8100 },
  @{ key='d3';    pdb=(Get-ChildItem $TargetsDir -Filter "3PBL*.pdb" | Select-Object -First 1).FullName; center=' 0.085,-14.828,10.432'; seedBase=8200 },
  @{ key='5ht2b'; pdb=(Get-ChildItem $TargetsDir -Filter "4IB4*.pdb"  | Select-Object -First 1).FullName; center=' 22.448,18.284,11.726'; seedBase=8300 }
)
$end = Get-Date '07:30'
if ($end -lt (Get-Date)) { $end = $end.AddDays(1) }   # 次日 07:30 (旧写法会得到"今天07:30"或 +11h 误判)
$counts = @{ b2ar=0; d3=0; '5ht2b'=0 }
$totals = @{ b2ar=0; d3=0; '5ht2b'=0 }
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE3_GEN_START ckpt={1} end={2}" -f (Get-Date -Format 'HH:mm'), $(if ($useFt) { 'finetuned' } else { 'pretrained' }), $end.ToString('MM-dd HH:mm'))
while ((Get-Date) -lt $end) {
  $free = (Get-PSDrive D).Free / 1GB
  if ($free -lt 12) {
    Get-ChildItem (Join-Path $Root 'outputs') -Recurse -File -Filter "samples_*.pt" -ErrorAction SilentlyContinue | Where-Object { $_.Name -ne 'samples_all.pt' } | Remove-Item -Force -ErrorAction SilentlyContinue
    if (((Get-PSDrive D).Free / 1GB) -lt 6) { Start-Sleep -Seconds 900; continue }
  }
  foreach ($t in $targets) {
    if ((Get-Date) -ge $end) { break }
    $key = $t.key; $counts[$key]++
    $seed = $t.seedBase + $counts[$key]
    $cfg = Join-Path $ConfigsDir 'gen_batch.yml'
    $extra = @()
    if ($useFt -and $ckpt) { $extra = @("--template", ".\configs\abft_finetuned.yml") }
    & $py $gen --template $(if ($useFt -and $ckpt) { ".\configs\abft_finetuned.yml" } else { $tmpl }) --out $cfg --seed $seed --diversity-w 0.5 --guided 1 | Out-Null
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path $cfg)) { Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} {1} CONFIG_GEN_FAIL -> 中止" -f (Get-Date -Format 'HH:mm'), $key); break }
    $t0 = Get-Date
    & $py $samplePy --pdb_path $t.pdb --center $t.center --config $cfg --outdir "$out\gen_$key\batch_$($counts[$key])" *>&1 | Out-File "$env:TEMP\gen_${key}_$($counts[$key]).log" -Encoding UTF8
    $rc = $LASTEXITCODE
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
    Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} {1} batch_{2} seed={3} finished={4} total={5} exit={6}" -f (Get-Date -Format 'HH:mm'), $key, $counts[$key], $seed, $n, $totals[$key], $rc) -ErrorAction SilentlyContinue
  }
}
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} OVERNIGHT_PHASE2_DONE b2ar={1} d3={2} 5ht2b={3}" -f (Get-Date -Format 'HH:mm'), $totals['b2ar'], $totals['d3'], $totals['5ht2b'])