# Thermal Transport Properties of Ag<sub>8</sub>TS<sub>6</sub> (T = Si, Ge, Sn) Argyrodites

This repository includes the outputs and analysis scripts related to the thermal transport properties of Ag<sub>8</sub>TS<sub>6</sub> (T = Si, Ge, Sn) argyrodites.

The corresponding raw data required to reproduce these results are available in our [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17399976.svg)](https://doi.org/10.5281/zenodo.17399976) and [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17397457.svg)](https://doi.org/10.5281/zenodo.17397457) repositories. For convenience, selected output files are also included here.

---
## Experimental Data

- It includes all experimental datasets for the Ag<sub>8</sub>TS<sub>6</sub> (T = Si, Ge, Sn) argyrodites discussed in the manuscript and Supporting Information.

## Group Velocity


- `Ag8TS6/velocity.py`: extracts group velocities from the `mesh.yaml` files.  
- `Ag8TS6/plot_velocity_freq.py`: plots the group velocity as a function of frequency.

---

## Grüneisen Parameter

- Use the provided `phonon.yml` file to create a Python environment and install all necessary dependencies.  
- `Grüneisen_parameter.ipynb`: example notebook for calculating the Grüneisen parameters using **pymatgen**.  
- `Grüneisen_plots.ipynb`: notebook for generating customized Grüneisen parameter plots as a function of frequency.

---

## Lattice Thermal Conductivity (κ)

- This directory includes the models used to compute the **lattice thermal conductivity** of the Ag<sub>8</sub>TS<sub>6</sub> argyrodites. Please refer to the `README` for further details on the models.

---

## Phonon Lifetime

- `Tau_comparison.ipynb`: notebook for comparing phonon lifetimes obtained from the Grüneisen model, MLIP-MACE-MP-03b, and analytical models.

---
