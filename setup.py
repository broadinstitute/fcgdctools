#TBD
# long_description from README
# version - either from version file or from python function 

from setuptools import setup, find_packages



setup(
    name = 'fcgdctools',
    packages = ['fcgdctools'],
    version = '0.1.3',
    description = 'Utilities for integrating FireCloud and GDC',
    author = 'Chet Birger',
    author_email = 'birger@broadinstitute.org',
    url = 'https://github.com/broadinstitute/fcgdctools',
    classifiers = [
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Development Status :: 2 - Pre-Alpha"
        ],
    entry_points={
        'console_scripts': [
            'genFcWsLoadFiles=fcgdctools.fc_loadfiles:main',
        ],
    },
    install_requires=['requests']
)    
