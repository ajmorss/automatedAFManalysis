# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 19:14:31 2014

@author: Andrew
"""

import sys
from cx_Freeze import setup, Executable

sys.path.append('C:\\Users\\heinrichlab\\Desktop\\Andrew\\Python')
base = None
if sys.platform == 'win32':
    base = 'Win32GUI'
includes = ['ForceDatAFM']
options = {
    'build_exe': {'packages': ['scipy', "matplotlib.backends.backend_tkagg", 'sqlalchemy', 'MySQLdb', 'sklearn'],
        'excludes': [] ,
        'includes': includes,
        "include_files": [('C:\\Program Files\\WinPython-64bit-2.7.6.3\\python-2.7.6.amd64\\Lib\\site-packages\\scipy\\special\\_ufuncs.pyd','_ufuncs.pyd')] # Sometimes a little finetuning is needed
    }
}

executables = [
    Executable('Single_Bond_Automatic_Analysis_App.py', base=base)
]

setup(name='Data Viewer',
      version='0.1',
      description='SBAAA',
      executables=executables,
      options=options
      )
