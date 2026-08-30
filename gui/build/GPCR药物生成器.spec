# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['../app_easy.py'],
    pathex=[],
    binaries=[('D:/Miniconda3/envs/Pocket2Mol/Library/bin/tcl86t.dll', '.'), ('D:/Miniconda3/envs/Pocket2Mol/Library/bin/tk86t.dll', '.')],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='GPCR药物生成器',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
