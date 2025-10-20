from mace.calculators import MACECalculator, mace_mp
from phono3py import Phono3py
from pymatgen.io.phonopy import get_pmg_structure, get_phonopy_structure
from pymatgen.core.structure import Structure
import numpy as np
from ase.filters import FrechetCellFilter
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import Trajectory as AseTrajectory
from ase.optimize import BFGS, FIRE, LBFGS, BFGSLineSearch, LBFGSLineSearch, MDMin
from mace.calculators import MACECalculator

# set a device
device = "cuda"
strct = Structure.from_file("path/to/the/structure")

print(strct)
unitcell = strct.to_ase_atoms()
unitcell.calc= MACECalculator(model_path='mace-mp-0b3-medium.model', device='cuda') # include the MLIP model 

unitcell_filter = FrechetCellFilter(unitcell)
opt=LBFGS(unitcell_filter)
traj= AseTrajectory("save_AgSnS.traj","w",unitcell )
opt.attach(traj)

opt.run(fmax=1e-5, steps=30000)


unitcell=AseAtomsAdaptor.get_structure(unitcell)
print(unitcell)
unitcell=get_phonopy_structure(unitcell)
ph3 = Phono3py(unitcell, supercell_matrix=[1, 2, 2], primitive_matrix='auto')

ph3.generate_displacements()
set_of_forces = []
dummy_force = np.zeros((len(ph3.supercells_with_displacements[0].get_scaled_positions()), 3))
print(len(ph3.supercells_with_displacements))
for displacement in ph3.supercells_with_displacements:
    if displacement is not None:
        struct = get_pmg_structure(displacement)
        cell = struct.to_ase_atoms()

        
        # Copyright by Atsushi Togo
        cell.set_calculator(MACECalculator(model_path='mace-mp-0b3-medium.model', device='cuda')) # include the MLIP model 
        forces = cell.get_forces()
        drift_force = forces.sum(axis=0)
        print(("[Phonopy] Drift force:" + "%11.5f" * 3) % tuple(drift_force))
        for force in forces:
            force -= drift_force / forces.shape[0]
        set_of_forces.append(forces)
    else:
        set_of_forces.append(dummy_force)
ph3.save("phono3py.yaml")



from phono3py.interface.phono3py_yaml import Phono3pyYaml

ph3yml = Phono3pyYaml()

ph3yml.read("phono3py.yaml")

disp_dataset = ph3yml.dataset
ph3.dataset = disp_dataset


set_of_forces = np.array(set_of_forces)
ph3.forces = set_of_forces.reshape(-1, len(ph3.supercell), 3)

ph3.produce_fc3()

ph3.save("results.yaml")

ph3.mesh_numbers = [6, 14, 10]
ph3.init_phph_interaction()
temperatures=range(0,600,10)

cmd="phono3py-load results.yaml --lbte --wigner --tmin 0 --tmax 600 --tstep 10  --mesh 6 14 10 --wgp"
import os
os.system(cmd)

import yaml

with open('ir_grid_points.yaml', 'r') as file:
    data = yaml.safe_load(file)

    # Extract the list of grid_point values
    grid_points = [entry['grid_point'] for entry in data.get('ir_grid_points', [])]

print(grid_points)

for grid in grid_points:
    ph3.run_thermal_conductivity(temperatures=temperatures,is_isotope=False, boundary_mfp=100000, conductivity_type="wigner",is_LBTE=False, write_pp=True,grid_points=[grid])

ph3.run_thermal_conductivity(temperatures=temperatures,is_isotope=False, boundary_mfp=100000, conductivity_type="wigner",is_LBTE=False, read_pp=True,write_kappa=True)


