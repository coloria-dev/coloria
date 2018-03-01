# -*- coding: utf-8 -*-
#
import codecs
import os

from setuptools import setup, find_packages

# https://packaging.python.org/single_source_version/
base_dir = os.path.abspath(os.path.dirname(__file__))
about = {}
with open(os.path.join(base_dir, 'colorio', '__about__.py'), 'rb') as f:
    # pylint: disable=exec-used
    exec(f.read(), about)


def read(fname):
    try:
        content = codecs.open(
            os.path.join(base_dir, fname), encoding='utf-8'
            ).read()
    except FileNotFoundError:
        content = ''
    return content


setup(
    name='colorio',
    version=about['__version__'],
    packages=find_packages(),
    package_data={'colorio': [
        'data/ebner_fairchild.yaml',
        'data/hung-berns/*.yaml',
        'data/munsell/real.yaml',
        'data/macadam1942/*.yaml',
        'data/observers/*.yaml',
        'data/illuminants/*.yaml',
        'data/luo-rigg/*.yaml',
        ]},
    url='https://github.com/nschloe/colorio',
    download_url='https://pypi.python.org/pypi/colorio',
    author=about['__author__'],
    author_email=about['__email__'],
    install_requires=[
        'matplotlib',
        'meshio',
        'meshzoo',
        'numpy',
        ],
    description='tools for color models',
    long_description=read('README.rst'),
    license=about['__license__'],
    classifiers=[
        about['__license__'],
        about['__status__'],
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics'
        ]
    )
