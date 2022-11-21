# Muscle Simulation 

We seek to implement the algorithm proposed by V. Modi et al. in their paper [EMU - Efficient Muscle Simulation](https://www.dgp.toronto.edu/projects/efficient-muscles/emu.pdf).

# Installing Dolfinx

In order to accomplish this project, we must install fenicsx.


### Create virtual environment

To start off, let's create the virtual environment in the root of our `muscle-simulation` repo:

`$ python3 -m venv env`

### For **Debian**:
- `$ sudo apt install fenicsx`

### For **Ubuntu**:
- `$ add-apt-repository ppa:fenics-packages/fenics`
- `$ apt update`
- `$ apt install fenicsx`

# C++ dependenicies

You will need:
- C++ compiler
- Boost
- CMake
- pkg-config
- FFCx
- Basix
- Pugixml
- PETSc

## Boost:

To install Boost:

`$ sudo apt install libboost-all-dev`

## CMake:

To install CMake:

`$ sudo apt install cmake`

## FFCx

To install FFCx:

- `$ git clone https://github.com/FEniCS/ffcx ~/ffcx`
- `$ cd ~/ffcx`
- `$ cmake -B build-dir -S cmake/` 
- `$ cmake --build build-dir` 
- `$ cmake --install build-dir`

## pkg-config:

To install pkg-config:

`$ sudo apt install pkg-config`

## Basix:

To install Basix:

- `$ git clone git@github.com:FEniCS/basix.git ~/basix`
- `$ cd ~/basix/cpp`
- `$ cmake -DCMAKE_BUILD_TYPE=Release -B build-dir -S .`
- `$ cmake --build build-dir`
- `$ sudo cmake --install build-dir`

## Pugixml:

Download and extract the archive from the [pugixml official website](https://pugixml.org/).

To install Pugixml:

- `$ tar -xf pugixml-1.13.tar.gz`
- `$ cd pugixml-1.13/`
- `$ mkdir build/`
- `$ cd build/`
- `$ cmake ..`
- `$ make install`

## PETSc:

*Note*: Make sure the `mpicc` and `mpic++` compilers are installed on your system! They should be already installed by default.

To install PETSc:

- `$ git clone -b release https://gitlab.com/petsc/petsc.git ~/petsc`
- `$ cd ~/petsc`
- `$ ./configure --with-cc=mpicc --with-cxx=mpic++ --with-fc=0`
- `$ make all check`
- `$ export PETSC_DIR=/home/$USER/petsc`
- `$ export PETSC_ARCH=arch-linux-cxx-debug`

**Do** ***NOT*** **delete the petsc folder after you have installed it! The Python version will still need it!!!**

# Build dolfinx

Now that we have installed the C++ dependencies for Dolfinx, we can start building the C++ core components.

Finally install the Python version onto our virtual environment, using pip.

## Build C++ core components

We assume that we clone dolfinx in the home directory `~`

- `$ git clone git@github.com:FEniCS/dolfinx.git`
- `$ cd ~/dolfinx/`
- `$ mkdir cpp/build`
- `$ cd cpp/build`
- `$ cmake ..`
- `$ sudo make install`

## Python installation

Once done, we can start installing Dolfinx for Python, as well as the requirements.

We change the directory to the root of our `muscle-simulation` git repo.

Assuming it is located in the home directory `~`

- `$ cd ~/muscle-simulation/`
- `$ source env/bin/activate`
- `$ pip install -r requirements.txt`
- `$ cd ~/dolfinx/python/`
- `$ pip install .`

With this, Dolfinx should now be installed on your virtual environment!

# Roadmap

### Data processing
* [x] Tetrahedralize 2D mesh
* [ ] Labelize meshes

### Algorithm implementation
* [ ] Implement neo-Hookean function

`... more to come ...`

* [ ] Implement EMU algorithm
