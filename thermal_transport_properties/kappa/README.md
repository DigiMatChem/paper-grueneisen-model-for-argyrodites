This repository provides scripts to compute two-channel–based lattice thermal conductivity (κ) using different approaches.

# Analytical Model 

- `Fitting.py` fits experimental data and evaluates κ using the two-channel model.

# Derived from elasticity calc
- Used the provided `elasticity_env.yml` to create an python environment to install necessary dependencies
- `run_elasticity_workflow.py` example script to create workflows to compute elasticity for a structure.
- `get_v_and_kappa_min.ipynb` reads the locally saved output Taskdoc of atomate2 elasticity workflow and extracts sound velocity and computes minimum lattice thermal conductivy using Cahill and Agne model 

# Grüneisen Model 

-`conductivity.py` and `get_kappa.py` compute the lattice thermal conductivity components starting from the Xia (Minikappa) model, but with phonon lifetimes calculated from the Grüneisen parameter. 
-`collect_kappa.py` collects kappa components

# Minimum lattice thermal conductivity using harmonic phonon data (Xia model)

- Code to compute this is available [here](https://github.com/yimavxia/Minikappa)