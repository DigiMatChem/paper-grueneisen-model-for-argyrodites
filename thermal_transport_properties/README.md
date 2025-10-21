# Thermal Transport Properties of Ag<sub>8</sub>TS<sub>6</sub> (T = Si, Ge, Sn) Argyrodites

This repository contains the outputs and analysis scripts related to the thermal transport properties of Ag<sub>8</sub>TS<sub>6</sub> (T = Si, Ge, Sn) argyrodites.

Raw data necessary for some of this scripts is available in the linked Zenodo repository.

---

## Group Velocity

- `Ag8TS6/velocity.py`: extracts group velocities from the `mesh.yaml` files.  
- `Ag8TS6/plot_velocity_freq.py`: plots the group velocity as a function of frequency.

---

## Grüneisen Parameter
Grüneisen parameter raw data are provided in [ZENODO](https://zenodo.org/uploads/17399976). This directory includes all necessary files to generate the gruneisen_mesh.yaml file using the provided script `run_gru_mesh.sh` (located in each respective argyrodite folder). The resulting files are subsequently used to compute the average Grüneisen parameter and produce related plots.

- Use the provided `phonon.yml` file to create a Python environment and install all necessary dependencies.  
- `Grüneisen_parameter.ipynb`: example notebook for calculating the Grüneisen parameters using **pymatgen**.  
- `Grüneisen_plots.ipynb`: notebook for generating customized Grüneisen parameter plots as a function of frequency.

---

## Lattice Thermal Conductivity (κ)

- This directory contains all the models used to compute the **lattice thermal conductivity** of the Ag<sub>8</sub>TS<sub>6</sub> argyrodites.

---

## Phonon Lifetime

- `Tau_comparison.ipynb`: notebook for comparing phonon lifetimes obtained from the Grüneisen model, MLIP (MACE-MP-03b), and analytical models.

---
