#!/usr/bin/env python
## from distutils.core import setup
from numpy.distutils.core import setup, Extension
with open('lagrange/VERSION') as f: VERSION = f.read().strip()

description = '''\
Lagrange is a Python package implementing likelihood models for
geographic range evolution on phylogenetic trees, with methods for
inferring rates of dispersal and local extinction and ancestral
ranges.'''

## def configure_expokit(parent_package='', top_path=None):
##     from numpy.distutils.misc_util import Configuration, get_info
##     config = Configuration('lagrange', parent_package, top_path)
##     src = ['lagrange/expokit.f90']
##     config.add_extension(
##         'expokit', sources=src, extra_link_args=['-llapack', '-lblas']## ,
##         ## f2py_options=['--link-lapack', '--link-blas']
##         )
##     return config

expokit = Extension(name='expokit', sources=['lagrange/expokit.f90'],
                    extra_link_args=['-llapack', '-lblas'])
    
setup(
    name='Lagrange',
    description='Likelihood analysis of geographic range evolution',
    license='GPL',
    platforms='ALL',
    long_description=description,
    version=VERSION,
    author='Richard Ree',
    author_email='rree@fieldmuseum.org',
    url='https://github.com/rhr/lagrange-python',
    packages=['lagrange'],
    package_data=dict(lagrange=['VERSION']),
    ext_modules=[expokit]
    )

## if __name__ == '__main__':
##     from numpy.distutils.core import setup
##     setup(configuration=configure_expokit)
