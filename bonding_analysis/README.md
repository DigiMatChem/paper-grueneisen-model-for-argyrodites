
# Pre requisites

- use the provided `requirements.txt` to install the necessary dependencies in a conda env with python v3.10

# Descriptions of scripts

Raw data necessary for this scripts is available in the linked [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17397457.svg)](https://doi.org/10.5281/zenodo.17397457) repository.

- `two_center/bonding_analysis.ipynb` script will reproduce the Fig 2 of main manuscript. One just needs set the appropriate `calc_dir` path.
- To visualize and save the Summary of Local environments detected using bond strengths as criterion use the `two_center/app.py` script. Please use `python two_center/app.py --help` for details on its usage. Hovering over the bonds in the app will show corresponding neighbor *ICOHPs* or *ICOBIs* depending on the file being read. Can be used to reproduce the Figures and Tables S2-S5 of SI.
- `multi_center/ncicobi_bonding_analysis.ipynb` script will reproduce the Fig 3 of main manuscript. One just needs set the appropriate path with the LOBSTER data 

### Note on Spin Channels in the Multi-Center LOBSTER Output

The output files resulting from multicenter bonding calculations produced by LOBSTER v4.1.0 occasionally report SPIN down-channel values inaccurately. Consequently, the SPIN up-channel was duplicated during analysis to get the correct multicenter ICOBI values.
