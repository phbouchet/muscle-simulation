# Muscle Simulation 

We seek to implement the algorithm proposed by V. Modi et al. in their paper [EMU - Efficient Muscle Simulation](https://www.dgp.toronto.edu/projects/efficient-muscles/emu.pdf)

## Setting up the virtual environment

To setup the virtual environment:

* `python3 -m venv env`
* `source env/bin/activate`
* `pip install -r requirements.txt`

Once setup, the environment does not need to be modified, and can be automatically detected with certain IDEs, such as Visual Studio Code.

To exit the virtual environment, execute:

* `deactivate`

## Project roadmap


### Data processing
* [x] Tetrahedralize 2D mesh
* [ ] Labelize meshes

### Algorithm implementation
* [ ] Implement neo-Hookean function
`... more to come ...`
* [ ] Implement EMU algorithm
