param([string]$outDir = "D:\MMModel\outputs\ppt_png\v3")
Get-Process POWERPNT -ErrorAction SilentlyContinue | Stop-Process -Force
Start-Sleep -Seconds 2
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
$src = "C:\Users\31908\Desktop\gpcr (1).pptx"
$ppt = New-Object -ComObject PowerPoint.Application
$pres = $ppt.Presentations.Open($src, $true, $false, $false)
$count = $pres.Slides.Count
for ($n = 1; $n -le $count; $n++) {
  $f = Join-Path $outDir ("p{0:d2}.png" -f $n)
  if (Test-Path $f) { Remove-Item $f -Force }
  $pres.Slides.Item($n).Export($f, "PNG", 1600, 900)
}
$pres.Close()
$ppt.Quit()
[System.Runtime.InteropServices.Marshal]::ReleaseComObject($ppt) | Out-Null
[GC]::Collect()
Write-Output ("exported " + $count + " slides to " + $outDir)
