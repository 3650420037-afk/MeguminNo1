$ErrorActionPreference = 'Continue'
$repo = "D:\MMModel\Pocket2Mol"
$py = "D:\Miniconda3\envs\Pocket2Mol\python.exe"
$out = "$repo\outputs\overnight_phase2"
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
Set-Location $repo
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE4_RETRAIN_START" -f (Get-Date -Format 'HH:mm'))
& $py train.py --config ".\configs\train_multitarget_v2.yml" --logdir ".\logs" *>&1 | Out-File "$env:TEMP\phase4_train.log" -Encoding UTF8
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE4_RETRAIN_EXIT={1}" -f (Get-Date -Format 'HH:mm'), $LASTEXITCODE)
# 数值排序 + 非数字名容错 (旧写法 Sort-Object {[int]BaseName} 遇 best.pt 会抛错)
$ckptDir = Get-ChildItem "$repo\logs" -Directory -ErrorAction SilentlyContinue | Where-Object Name -like 'train_multitarget_v2*' | Sort-Object Name -Descending | Select-Object -First 1
$ckpt = $null
if ($ckptDir) {
  $ckpt = Get-ChildItem (Join-Path $ckptDir.FullName "checkpoints") -Filter "*.pt" -ErrorAction SilentlyContinue |
          Where-Object { $_.BaseName -match '^\d+$' } |
          Sort-Object { [int]$_.BaseName } -Descending | Select-Object -First 1
}
Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("CKPT={0}" -f ($(if ($ckpt) { $ckpt.FullName } else { "" })))
if ($ckpt) {
  $pdb = (Get-ChildItem "D:\MMModel\靶点结构" -Filter "4EIY*.pdb" | Select-Object -First 1).FullName
  & $py sample_for_pdb.py --pdb_path $pdb --center " -0.4,8.5,17.1" --config ".\configs\sample_for_pdb_ab_baseline_305020.yml" --outdir "$out\re_abft_baseline" *>&1 | Out-File "$env:TEMP\re_abft_b.log" -Encoding UTF8
  $src = [IO.File]::ReadAllText("$repo\configs\sample_for_pdb_ab_local_305020.yml").TrimStart([char]0xFEFF)
  $ft = ($src -replace 'checkpoint: [^\r\n]+', ("checkpoint: " + ($ckpt.FullName -replace '\\','/'))) -replace 'seed: \d+', 'seed: 2024'
  [IO.File]::WriteAllText("$repo\configs\re_abft_finetuned.yml", $ft, (New-Object System.Text.UTF8Encoding($false)))
  & $py sample_for_pdb.py --pdb_path $pdb --center " -0.4,8.5,17.1" --config ".\configs\re_abft_finetuned.yml" --outdir "$out\re_abft_finetuned" *>&1 | Out-File "$env:TEMP\re_abft_f.log" -Encoding UTF8
  & $py scripts\judge_promotion.py "$out\re_abft_baseline" "$out\re_abft_finetuned" "$out\re_decision.txt"
  Add-Content "$out\phase_log.txt" -Encoding UTF8 -Value ("{0} PHASE4_AB_DONE" -f (Get-Date -Format 'HH:mm'))
}