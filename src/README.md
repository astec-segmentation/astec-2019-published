# astec-package
Adaptive Segmentation and Tracking of Embryonic Cells 

## Authors:

- Emmanuel Faure
- Leo Guignard
- Gregoire Malandain
- Gael Michelin

# 1. CONTENTS

The folder contains the following elements:
- definitions.py: the definitions of the folder process and the RAWDATA parameters (time / angle etc.)
- ASTEC: the folder  containing all the libraries 
- 1-fuse.py: script to fuse raw data from light-sheet microscope
- 2-mars.py: marker-based watershed segmentation
- 3-manualcorrection.py:
- 4-astec.py: segmentation by propagation
- 5-postcorrection.py
- 6-named.py
- 7-virtualembryo.py 
- checklineage.py : tools for check lineage safety 
- rename.py : tools to rename file 
- README : the licence terms you accept by using the workflow
- license.txt   : the licence terms you accept by using the workflow

# 2. INSTALLATION

This code has been tested and used both on MacOs and Linux.

## 2.1 C code

### Authors:

- Gregoire Malandain
- Gael Michelin

### 2.1.1 Requirements

* a c compiler: Should be installed, refer with command line (on linux)

        dpkg --list | grep compiler

    to get list of avaiblable C compiler version or install gcc (on
    linux) with

        sudo apt-get install gcc

* cmake, a software to build source code (http://www.cmake.org)

        sudo apt-get install cmake

* zlib, a compression library	(http://www.zlib.net) :

        sudo apt-get install zlib1g-dev

### 2.1.2 Compilation / installation

* Compilation of the main code

        cd <astec-package>/ASTEC/CommunFunctions/cpp/vt/
        mkdir build
        cd build
        cmake ../
        make

    Executables will be built in `<astec-package>/ASTEC/CommunFunctions/cpp/vt/build/bin`

## 2.2 Python code

### Authors:

- Emmanuel Faure
- Leo Guignard
- Gregoire Malandain
- Gael Michelin

### 2.2.1 Requirements

* python 2.7 (python 3.x is not supported)

* pip, an installer for python (https://pypi.python.org/pypi/pip)

        sudo apt-get install python-pip python-dev build-essential


* zlib, a compression library (http://www.zlib.net):

        sudo apt-get install zlib1g-dev

* the libraries numpy, scipy (http://www.numpy.org, http://www.scipy.org)

        sudo apt-get install python-numpy python-scipy

* libhdf5-dev, cython, h5py  a library to read/write hdf5 images (https://www.hdfgroup.org/HDF5/):

        sudo apt-get install libhdf5-dev
        sudo pip install cython 
	    sudo pip install h5py

### 2.2.2 Compilation / installation

* pylibtiff

        cd <astec-package>/ASTEC/CommunFunctions/libtiff-0.4.0/
        sudo python setup.py install



















