from setuptools import setup

setup(
    name = 'constants',
    version = '0.1.0',    
    description = 'Handy package for easy administration of physical constants and unit conversions',
    url = 'https://github.com/JonKDahl/physical_constants_and_units',
    author = 'Jon Kristian Dahl',
    author_email = 'jonkd@uio.no',
    license = 'The Unlicense',
    packages = ['constants'],
    install_requires = ['numpy'],

    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)