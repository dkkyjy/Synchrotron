import numpy as np
from distutils.core import setup, Extension

Synchrotron_module = Extension('_Synchrotron',
                             sources=['Synchrotron_wrap.c', 'Synchrotron.c'],
                             include_dirs=['/usr/local/include', np.get_include()],
                             library_dirs=['/usr/local/lib'],
                             libraries=['gsl', 'gslcblas', 'numpy'],
                            )

setup(name = 'Synchrotron',
      version = '0.1',
      author = 'DuanKaikai',
      description = 'Synchrotron emission',
      ext_modules = [Synchrotron_module],
      py_modules = ['Synchrotron'],
     )
