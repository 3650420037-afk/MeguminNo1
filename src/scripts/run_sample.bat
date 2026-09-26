@echo off
chcp 65001 >nul
setlocal
REM src\scripts\run_sample.bat -> 上两级 = 仓库根 (与 paths.py 的 ROOT 一致)
pushd "%~dp0..\.."
set "Root=%CD%"
REM 解释器: 优先环境变量 EONMOL_PYTHON (与 paths.PYTHON 同源), 否则用 PATH 上的 python
set "py=%EONMOL_PYTHON%"
if not defined py set "py=python"
echo ========================================
echo   Pocket2Mol 一键分子生成脚本
echo ========================================
echo.
echo  [用法] 用记事本编辑本文件,修改下面两行:
echo    1. PDB 路径  : 换成您的蛋白结构文件
echo    2. --center  : 换成您的口袋中心坐标
echo       ★ 注意:第一个数字前必须有一个空格!
echo.
echo  当前配置(示例蛋白 4yhj.pdb):
echo    PDB   : data\example\4yhj.pdb
echo    CENTER: " 32.0,28.0,36.0"
echo.
echo  开始运行...
echo ========================================
echo.

"%py%" "src\sample_for_pdb.py" ^
    --pdb_path "data\example\4yhj.pdb" ^
    --center " 32.0,28.0,36.0" ^
    --config "configs\sample_for_pdb.yml" ^
    --outdir "outputs"

echo.
echo ========================================
echo   运行结束!结果在 outputs\ 目录下
echo ========================================
popd
pause
