
from setuptools import setup
import os
import codecs
import sys
import site
site.ENABLE_USER_SITE = "--user" in sys.argv[1:]


# https://packaging.python.org/single_source_version/
base_dir = os.path.abspath(os.path.dirname(__file__))

setup(
    name='tikzplotlib',
    version="0.9.8",
    packages=['tikzplotlib'],
    url='https://github.com/nschloe/tikzplotlib',
    download_url='https://pypi.python.org/pypi/tikzplotlib',
    install_requires=[
        'matplotlib >=1.4.0',
        'numpy',
        'Pillow >= 3.0.0',
        'six',
        ],
    description='convert matplotlib figures into TikZ/PGFPlots',
    classifiers=[
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Multimedia :: Graphics :: Graphics Conversion',
        'Topic :: Scientific/Engineering :: Visualization'
        ]
    )
