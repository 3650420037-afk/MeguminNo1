# 构建 GPCR 药物生成器 (两种模式)
#   1) 单文件版  dist\7-eonmol.exe             —— 便于分发, 启动需解压 (较慢)
#   2) 快速启动版 dist\fast\7-eonmol\...exe    —— 免解压目录版, 启动快数倍
# 用法: powershell -ExecutionPolicy Bypass -File src\gui\build_exe.ps1
$ErrorActionPreference = "Stop"

# 仓库根: src\gui\build_exe.ps1 -> src\gui -> src -> 仓库根
$Root   = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$GUI    = Join-Path $Root "src\gui"
# 解释器: 优先环境变量 EONMOL_PYTHON (与 src/paths.py 的 PYTHON 同一覆盖方式),
# 其次当前激活的 conda 环境, 最后 PATH 上的 python
$PY = $env:EONMOL_PYTHON
if (-not $PY -and $env:CONDA_PREFIX) { $PY = Join-Path $env:CONDA_PREFIX "python.exe" }
if (-not $PY) { $PY = (Get-Command python.exe).Source }
# Tcl/Tk 运行库位于解释器所在环境的 Library\bin
$TCLBIN = Join-Path (Split-Path -Parent $PY) "Library\bin"
if (-not (Test-Path $TCLBIN)) {
    throw "解释器 $PY 不带 Tcl/Tk (找不到 $TCLBIN)。请先激活 Pocket2Mol 环境, 或设置环境变量 EONMOL_PYTHON 指向该环境的 python.exe 后重试。"
}
$NAME   = "7-eonmol"

# 明确不用的模块: 排除后体积更小、启动更快
$EXCLUDES = @(
    "--exclude-module", "numpy", "--exclude-module", "torch", "--exclude-module", "rdkit",
    "--exclude-module", "PIL", "--exclude-module", "yaml", "--exclude-module", "pptx",
    "--exclude-module", "lmdb", "--exclude-module", "scipy", "--exclude-module", "pandas",
    "--exclude-module", "matplotlib", "--exclude-module", "IPython",
    "--exclude-module", "setuptools", "--exclude-module", "pip", "--exclude-module", "wheel",
    "--exclude-module", "unittest", "--exclude-module", "pydoc", "--exclude-module", "doctest",
    "--exclude-module", "lib2to3", "--exclude-module", "distutils",
    "--exclude-module", "ssl", "--exclude-module", "_ssl",
    "--exclude-module", "hashlib", "--exclude-module", "_hashlib",
    "--exclude-module", "sqlite3", "--exclude-module", "curses"
)

function Trim-TclData([string]$Root) {
    # Tcl/Tk 默认带完整时区/多语言/全编码数据 (500+ 文件), 启动时要逐个扫描 —— 裁掉
    $tcl = Join-Path $Root "_internal\_tcl_data"
    if (-not (Test-Path $tcl)) { $tcl = Join-Path $Root "_tcl_data" }
    if (-not (Test-Path $tcl)) { return }
    foreach ($d in @("tzdata", "msgs")) {
        $p = Join-Path $tcl $d
        if (Test-Path $p) { Remove-Item $p -Recurse -Force -ErrorAction SilentlyContinue }
    }
    $enc = Join-Path $tcl "encoding"
    if (Test-Path $enc) {
        $keep = @("ascii.enc", "cp1252.enc", "cp936.enc", "gb2312.enc", "gb12345.enc",
                  "iso8859-1.enc", "utf-8.enc", "unicode.enc", "cp437.enc", "macRoman.enc")
        Get-ChildItem $enc -File | Where-Object { $keep -notcontains $_.Name } |
            Remove-Item -Force -ErrorAction SilentlyContinue
    }
}

Write-Host "== 1/5 生成图标 ==" -ForegroundColor Cyan
& $PY "$GUI\make_icon.py"
if ($LASTEXITCODE -ne 0) { throw "图标生成失败" }

Write-Host "== 2/5 备份旧版 ==" -ForegroundColor Cyan
$Legacy = Join-Path $GUI "dist\legacy"
New-Item -ItemType Directory -Force -Path $Legacy | Out-Null
$OldExe = Join-Path $GUI "dist\$NAME.exe"
if (Test-Path $OldExe) {
    $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
    Move-Item -Force $OldExe (Join-Path $Legacy "${NAME}_$stamp.exe")
    Write-Host "  旧版已归档: legacy\${NAME}_$stamp.exe"
}
Remove-Item (Join-Path $GUI "dist\fast") -Recurse -Force -ErrorAction SilentlyContinue

$Common = @("--noconfirm", "--clean", "--windowed", "--icon", "$GUI\icon.ico",
            "--paths", (Join-Path $Root "src"),
            "--hidden-import", "ui_kit", "--hidden-import", "icon_b64",
            "--hidden-import", "paths",
            "--add-binary", "$TCLBIN\tcl86t.dll;.", "--add-binary", "$TCLBIN\tk86t.dll;.") + $EXCLUDES

Write-Host "== 3/5 打包单文件版 ==" -ForegroundColor Cyan
& $PY -m PyInstaller @Common --onefile --name $NAME `
    --distpath "$GUI\dist" --workpath "$GUI\build\onefile" --specpath "$GUI\build" `
    "$GUI\app_pro.py"
if ($LASTEXITCODE -ne 0) { throw "单文件打包失败" }

Write-Host "== 4/5 打包快速启动版 ==" -ForegroundColor Cyan
& $PY -m PyInstaller @Common --onedir --name $NAME `
    --distpath "$GUI\dist\fast" --workpath "$GUI\build\onedir" --specpath "$GUI\build" `
    "$GUI\app_pro.py"
if ($LASTEXITCODE -ne 0) { throw "快速启动版打包失败" }
$FastDir = Join-Path $GUI "dist\fast\$NAME"
Trim-TclData $FastDir

Write-Host "== 5/5 创建桌面快捷方式 ==" -ForegroundColor Cyan
$FastExe = Join-Path $FastDir "$NAME.exe"
if (Test-Path $FastExe) {
    $ws = New-Object -ComObject WScript.Shell
    $lnk = $ws.CreateShortcut((Join-Path ([Environment]::GetFolderPath("Desktop")) "$NAME.lnk"))
    $lnk.TargetPath = $FastExe
    $lnk.WorkingDirectory = $FastDir
    $lnk.IconLocation = "$GUI\icon.ico"
    $lnk.Description = "7-eonmol (免解压目录版, 启动更快)"
    $lnk.Save()
    Write-Host "  快捷方式: 桌面\$NAME.lnk"
}

$one = Join-Path $GUI "dist\$NAME.exe"
$oneMB = [math]::Round((Get-Item $one).Length / 1MB, 1)
$fastFiles = (Get-ChildItem $FastDir -Recurse -File | Measure-Object).Count
$fastMB = [math]::Round((Get-ChildItem $FastDir -Recurse -File | Measure-Object Length -Sum).Sum / 1MB, 1)
Write-Host ("完成:`n  单文件版   {0}  ({1} MB)`n  快速启动版 {2}  ({3} 个文件, {4} MB)" -f `
    $one, $oneMB, $FastExe, $fastFiles, $fastMB) -ForegroundColor Green



