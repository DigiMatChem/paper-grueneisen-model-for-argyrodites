This repository provides scripts to compute two-channel–based lattice thermal conductivity (κ) using different approaches.

# Analytical Model 

- `Fitting.py` fits experimental data and evaluates κ using the two-channel model.

# Grüneisen Model 

-`conductivity.py` and `get_kappa.py` compute the lattice thermal conductivity components starting from the Xia (Minikappa) model, but with phonon lifetimes calculated from the Grüneisen parameter. 
-`collect_kappa.py` collects kappa components

# Minikappa

- `conductivity.py` and `get_minikappa.py` calculate the minimal lattice thermal conductivity according to the Xia model. 