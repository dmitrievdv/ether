from numpy.distutils.core import Extension

ext = Extension(name = 'ethercalc', sources = ['fsrc/frill.f90'],
	 f2py_options = ['only:','compute', ':'])

if __name__ == '__main__':
	from numpy.distutils.core import setup
	setup(name = 'f2py-ethercalc', ext_modules = [ext])
					 # ('coords', dict(sources=['fsrc/coords.f90'])),
					 # ('surfmod', dict(sources=['fsrc/surfmod.f90'])),
					 # ('constants', dict(sources=['fsrc/constants.f90'])),
					 # ('inter2d', dict(sources=['fsrc/2dinter.f90']))] 