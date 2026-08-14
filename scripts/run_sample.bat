@echo off
chcp 65001 >nul
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
echo    PDB   : ./example/4yhj.pdb
echo    CENTER: " 32.0,28.0,36.0"
echo.
echo  开始运行...
echo ========================================
echo.

cd /d D:\MMModel\Pocket2Mol

D:\Miniconda3\envs\Pocket2Mol\python.exe sample_for_pdb.py ^
    --pdb_path ./example/4yhj.pdb ^
    --center " 32.0,28.0,36.0" ^
    --config ./configs/sample_for_pdb.yml ^
    --outdir ./outputs

echo.
echo ========================================
echo   运行结束!结果在 outputs\ 目录下
echo ========================================
pause

