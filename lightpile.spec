# -*- mode: python -*-

# Analysis figures out what to include (script, dependency modules,
# binaries).
a = Analysis(['lightpile\\app.py'],
             hiddenimports=[],
             hookspath=['.\\hooks\\'],
             runtime_hooks=None)
pyz = PYZ(a.pure)

# EXE contains what will be in the windows executable
#   'exclude_binaries': should binaries be excluded of the executalbe (better
#                       put them into the directory)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='lightpile.exe',
          debug=False,
          strip=None,
          upx=True,
          console=True )

# COLLECT creates the output folder with the EXE and all other content.
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas + [('README',        'README',    'DATA'),
                          ('LICENSE',       'LICENSE',   'DATA'),
                          ('AGPLv3',        'AGPLv3',    'DATA'),
                          ('CHANGES',       'CHANGES',   'DATA'),
                          ('INSTALL',       'INSTALL',   'DATA')],
               strip=None,
               upx=True,
               name='lightpile')
