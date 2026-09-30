param([string]$outDir = "D:\MMModel\outputs\ppt_png\v4", [string]$slides = "10,15,16,22")
Get-Process POWERPNT -ErrorAction SilentlyContinue | Stop-Process -Force
Start-Sleep -Seconds 2
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
$src = "C:\Users\31908\Desktop\gpcr (1).pptx"
$ppt = New-Object -ComObject PowerPoint.Application
$pres = $ppt.Presentations.Open($src, $true, $false, $false)
foreach ($n in $slides.Split(",")) {
  $k = [int]$n
  $f = Join-Path $outDir ("p{0:d2}.png" -f $k)
  if (Test-Path $f) { Remove-Item $f -Force }
  $pres.Slides.Item($k).Export($f, "PNG", 1600, 900)
}
$pres.Close(); $ppt.Quit()
Write-Output ("exported: " + $slides)
