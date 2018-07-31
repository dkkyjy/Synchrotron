import numpy as np
from distutils.core import setup, Extension

Synchrotron_module = Extension('Synchrotron',
                             sources=['Synchrotron.c', 'ElectronDistribution.c'],
                             include_dirs=['/usr/local/include', np.get_include()],
                             library_dirs=['/usr/local/lib'],
                             libraries=['gsl', 'gslcblas'],
                            )

setup(name = 'Synchrotron',
      version = '0.1',
      author = 'Duan Kaikai',
      description = 'Synchrotron emission',
      ext_modules = [Synchrotron_module],
      py_modules = ['Synchrotron'],
     )
