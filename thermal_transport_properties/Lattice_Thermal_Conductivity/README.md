# This repository provides scripts to compute two-channel–based lattice thermal conductivity (κ) using different approaches.
---
## Analytical Model 

- `AM_comparison.ipynb` fits experimental data and evaluates κ using the two-channel model. Additionlly, it provides a comparison in the phonon channel subtracting point-defect and boundary scattering.
---

## Derived from elasticity calc

- Used the provided `elasticity_env.yml` to create an python environment to install necessary dependencies
- `run_elasticity_workflow.py` example script to create workflows to compute elasticity for a structure.
- `get_v_and_kappa_min.ipynb` reads the locally saved output Taskdoc of atomate2 elasticity workflow and extracts sound velocity and computes minimum lattice thermal conductivy using Cahill and Agne model 

---

## Grüneisen Model 

- `conductivity.py` and `get_kappa.py` compute the lattice thermal conductivity components starting from the Xia (Minikappa) model, but with phonon lifetimes calculated from the Grüneisen parameter. 
- `collect_kappa.py` collects kappa components

---

## MLIP_MACE-MP-03b

- Used the provided `phono3py_ML.yml` to create an python environment to install necessary dependencies
- `get_AgGeSb_results.py` example script to computa lattice thermal conductivity based on the foundationa ML potential MACE-MP-03b 
- `Get_kappa_values.ipynb` collects the kappa values from MACE-MP-03b computations 
---

## Minimum lattice thermal conductivity using harmonic phonon data (Xia model)

- Code to compute this is available [here](https://github.com/yimavxia/Minikappa)