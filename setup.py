#!/usr/bin/env python

from distutils.core import setup

with open('lagrange/VERSION') as f: VERSION = f.read().strip()

description = '''\
Lagrange is a Python package implementing likelihood models for
geographic range evolution on phylogenetic trees, with methods for
inferring rates of dispersal and local extinction and ancestral
ranges.'''

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
    package_data={'lagrange': ['VERSION']}
    )
