from numpy.distutils.core import Extension

ext = Extension(name = 'ethercalc', sources = ['fsrc/compute.f90'],
	 libraries = ['fbauhube'],
	 f2py_options = ['only:','frill', 'fundamental', ':'])

if __name__ == '__main__':
	from numpy.distutils.core import setup
	setup(name = 'f2py-ethercalc', ext_modules = [ext],
		  libraries = [('fbauhube', dict(sources=['fsrc/fbauhube.f90']))])
					 # ('coords', dict(sources=['fsrc/coords.f90'])),
					 # ('surfmod', dict(sources=['fsrc/surfmod.f90'])),
					 # ('constants', dict(sources=['fsrc/constants.f90'])),
					 # ('inter2d', dict(sources=['fsrc/2dinter.f90']))] 