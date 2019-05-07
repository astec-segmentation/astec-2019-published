from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='pkl2tlp',
    version='1.0',
    description='converter from pkl ASTEC output to tlp file',
    long_description=long_description,
    url='https://github.com/leoguignard/pkl2tlp',
    author='Leo Guignard',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
    ],
    install_requires=['numpy'],
    scripts=['pkl2tlp.py', 'pkl_and_correction2tlp.py']
)
