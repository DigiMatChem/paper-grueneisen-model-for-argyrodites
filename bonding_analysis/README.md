# Pre requisites

- use the provided `requirements.txt` to install the necessary dependencies in a conda env with python v3.10

# Descriptions of scripts

Raw data necessary for this scripts is available in the linked Zenodo repository.

- `two_center/bonding_analysis.ipynb` script will reproduce the Fig 2 of main manuscript. One just needs set the appropriate `calc_dir` path.
- To visualize and save the Summary of Local environments detected using bond strengths as criterion use the `two_center/app.py` script. Please use `python two_center/app.py --help` for details on its usage. Hovering over the bonds in the app will show corresponding neighbor *ICOHPs* or *ICOBIs* depending on the file being read. Can be used to reproduce the Figures and Tables S2-S5 of SI.
- `multi_center/ncicobi_bonding_analysis.ipynb` script will reproduce the Fig 3 of main manuscript. One just needs set the appropriate path with the LOBSTER data 
