$ErrorActionPreference = 'Continue'
# src/scripts/phase4_retrain.ps1 -> 上两级 = 仓库根 (与 paths.py 的 ROOT 一致)
$Root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$TargetsDir = Join-Path $Root 'data\targets'
$ConfigsDir = Join-Path $Root 'configs'
$SrcDir = Join-Path $Root 'src'
# 解释器: 优先环境变量 EONMOL_PYTHON (与 paths.PYTHON 同源), 否则用 PATH 上的 python
$py = if ($env:EONMOL_PYTHON) { $env:EONMOL_PYTHON } else { 'python' }
$samplePy = Join-Path $SrcDir 'sample_for_pdb.py'
$out = Join-Path $Root 'outputs\overnight_phase2'
# 等待三靶点生成结束 (带超时上限: 旧写法 while($true) 在 phase2 崩溃/未启动时会永久挂起)
$deadline = (Get-Date).AddHours(14)
while ($true) {
  $done = Select-String -Path "$out\phase_log.txt" -Pattern 'OVERNIGHT_PHASE2_DONE' -ErrorAction SilentlyContinue
  if ($done) { break }
  if ((Get-Date) -gt $deadline) {
    Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value "PHASE4_WAIT_TIMEOUT -> 退出" -ErrorAction SilentlyContinue
    "PHASE4_WAIT_TIMEOUT"; exit 1
  }
  Start-Sleep -Seconds 300
}
Set-Location $Root
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE4_RETRAIN_START" -f (Get-Date -Format 'HH:mm'))
& $py (Join-Path $Root 'train.py') --config (Join-Path $ConfigsDir 'train_multitarget_v2.yml') --logdir (Join-Path $Root 'logs') *>&1 | Out-File "$env:TEMP\phase4_train.log" -Encoding UTF8
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE4_RETRAIN_EXIT={1}" -f (Get-Date -Format 'HH:mm'), $LASTEXITCODE)
# 数值排序 + 非数字名容错 (旧写法 Sort-Object {[int]BaseName} 遇 best.pt 会抛错)
$ckptDir = Get-ChildItem (Join-Path $Root 'logs') -Directory -ErrorAction SilentlyContinue | Where-Object Name -like 'train_multitarget_v2*' | Sort-Object Name -Descending | Select-Object -First 1
$ckpt = $null
if ($ckptDir) {
  $ckpt = Get-ChildItem (Join-Path $ckptDir.FullName "checkpoints") -Filter "*.pt" -ErrorAction SilentlyContinue |
          Where-Object { $_.BaseName -match '^\d+$' } |
          Sort-Object { [int]$_.BaseName } -Descending | Select-Object -First 1
}
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("CKPT={0}" -f ($(if ($ckpt) { $ckpt.FullName } else { "" })))
if ($ckpt) {
  $pdb = (Get-ChildItem $TargetsDir -Filter "4EIY*.pdb" | Select-Object -First 1).FullName
  & $py $samplePy --pdb_path $pdb --center " -0.4,8.5,17.1" --config (Join-Path $ConfigsDir 'sample_for_pdb_ab_baseline_305020.yml') --outdir "$out\re_abft_baseline" *>&1 | Out-File "$env:TEMP\re_abft_b.log" -Encoding UTF8
  $src = [IO.File]::ReadAllText((Join-Path $ConfigsDir 'sample_for_pdb_ab_local_305020.yml')).TrimStart([char]0xFEFF)
  $ft = ($src -replace 'checkpoint: [^\r\n]+', ("checkpoint: " + ($ckpt.FullName -replace '\\','/'))) -replace 'seed: \d+', 'seed: 2024'
  [IO.File]::WriteAllText((Join-Path $ConfigsDir 're_abft_finetuned.yml'), $ft, (New-Object System.Text.UTF8Encoding($false)))
  & $py $samplePy --pdb_path $pdb --center " -0.4,8.5,17.1" --config (Join-Path $ConfigsDir 're_abft_finetuned.yml') --outdir "$out\re_abft_finetuned" *>&1 | Out-File "$env:TEMP\re_abft_f.log" -Encoding UTF8
  & $py (Join-Path $SrcDir 'scripts\judge_promotion.py') "$out\re_abft_baseline" "$out\re_abft_finetuned" "$out\re_decision.txt"
  Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE4_AB_DONE" -f (Get-Date -Format 'HH:mm'))
}
