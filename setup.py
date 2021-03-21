#TBD
# long_description from README
# version - either from version file or from python function 

from setuptools import setup, find_packages



setup(
    name = 'fcgdctools',
    packages = find_packages(),
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
            'gdcLoadFiles=fcgdctools.fc_loadfiles:main',
            'gdcWorkspace=fcgdctools.ws_builder:main'
        ],
    },
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=[
        'requests[security]',
        'firecloud>=0.16.31']
)    
