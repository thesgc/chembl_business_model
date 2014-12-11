#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'mnowotka'

import sys

try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

if sys.version_info < (2, 7, 3) or sys.version_info >= (2, 7, 9):
    raise Exception('ChEMBL software stack requires python 2.7.3 - 2.7.9')

setup(
    name='chembl_business_model',
    version='0.5.18',
    author='Michal Nowotka',
    platforms=['Linux'],
    author_email='mnowotka@ebi.ac.uk',
    description='Python package providing chembl webservices API.',
    url='https://www.ebi.ac.uk/chembl/',
    license='CC BY-SA 3.0',
    packages=['chembl_business_model',
              'chembl_business_model.models'],
    long_description=open('README.rst').read(),
    install_requires=[
        'requests',
        'chembl_core_model>=0.5.8',
        'django-celery',
        'BeautifulSoup==3.2.1',
        'Pillow>=2.2.1',
        'cairocffi>=0.5.1',
 #       'gdb>=0.2.2',
        'numpy>=1.8.0',
    ],
    include_package_data=False,
    classifiers=['Development Status :: 2 - Pre-Alpha',
                 'Environment :: Web Environment',
                 'Framework :: Django',
                 'Intended Audience :: Developers',
                 'License :: Creative Commons :: Attribution-ShareAlike 3.0 Unported',
                 'Operating System :: POSIX :: Linux',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Chemistry'],
    zip_safe=False,
)
