# Muscle Simulation - NumPy EMU implementation 

An implementation of the EMU algorithm proposed by V. Modi et al. in their paper [EMU - Efficient Muscle Simulation](https://www.dgp.toronto.edu/projects/efficient-muscles/emu.pdf). We have implemented the various functions and concepts defined in the paper.

Here is a summary of what the different files implement:

- `physics_E.py`:
    - Discretized energy formula (Eq.5)
    - Derivative of discretized energy (Eq.9)
    - Computationally optimized hessian using low-rank approximation (Eq.14)
    - Inverse hessian using Woodbury matrix identity (Eq.17)
- `physics_Ec.py`:
    - As-continuous-as-possible (ACAP) energy (Eq.3)
    - ACAP energy minimization solve (Eq.4)
    - Algorithm to generate the $G$ matrix
    - Derivative of the continuous energy (Eq.10)
- `physics_psi.py`:
    - Neo-Hookean elasticity functions for muscle, bone and tendon (Eq.8)
    - Derivatives and hessians of Neo-Hookean functions
- `physics_EMU.py`:
    - Implementation of EMU algorithm

There is a wrapper (`physics_lib.py`) which essentially includes all the physics functions, except the EMU algorithm, into one file. 

## Setup

### Prerequisites

`python3-venv` will be necessary in order to create the virtual environment. To install on Ubuntu/Debian:
```console
$ sudo apt-get update -y
$ sudo apt-get install -y python3-venv 
```

### Initializing the virtual environment

To setup the virtual environment:

```console
$ python3 -m venv env
$ source env/bin/activate
$ pip install -r requirements.txt
```

To enter the virtual environment:
```console
$ source env/bin/activate
```

To exit the virtual environment:
```
$ deactivate
```

## Running the algorithm

The program can be executed in two different modes: convergence mode, and iterative mode.

Convergence mode:
- By default, it is run in convergence mode, which means that the algorithm will automatically halt once it has converged.
- To execute the program on a mesh in convergence mode:
```console
$ ./src/main.py simple_tet.obj
```

Iterative mode:
- The algorithm can also be run in iterative mode, which sets an upper bound for the number of iterations that the algorithm will perform.
- *Note*: The algorithm can converge before it reaches the upper bound.
- To execute the program on a mesh in iterative mode:
```console
$ ./src/main.py simple_tet.obj -i n #with n an integer
```

The program will generate the meshes in the `data/` folder, formatted as `EMU_mesh_*_*.obj`. To remove the generated meshes, simply run:

```console
$ make clean
```

We provide a mesh of a single tetrahedron in `data/simple_tet.obj` which can be used for testing and debugging purposes.

## Viewing the results

To view the results, we suggest using blender to import the .obj files, and view manually. Unfortunately, we cannot animate the movement of the mesh using blender, so you will have to view each individual mesh by hand to study its deformation.

## Authors
- [Philippe Bouchet](https://github.com/sudomane): philippe.bouchet@epita.fr
- Jean-Baptiste Deloges: jean-baptiste.deloges@epita.fr
- Sebastien Barbier: sebastien1.barbier@epita.fr