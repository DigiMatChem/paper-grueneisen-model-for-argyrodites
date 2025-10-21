import numpy as np
from conductivity import class_kappa
import multiprocessing as mp
from tqdm.autonotebook import tqdm
"""
author: @yixia
An example script to perform predictions of kL_min for the diamond example at 600 K
"""

#kappa = class_kappa()
temp_range = list(range(300, 620, 20))
#kappa.get_minikappa_phonopy(mesh_in = [8, 24, 16],                        
#                            sc_mat = np.eye(3)*[1, 3, 2],
#                            pm_mat = [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
#                            list_temp = temp_range,
#                            name_pcell = "POSCAR-unitcell",
#                            name_ifc2nd = "FORCE_CONSTANTS",
#                            list_taufactor = [2.0])

def condt_temp(temp):
    kappa = class_kappa()
    kappa.get_minikappa_phonopy(mesh_in = [10, 19, 14],
                            sc_mat = np.eye(3)*[1, 3, 2],
                            pm_mat = [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                            list_temp = [temp],
                            name_pcell = "POSCAR-unitcell",
                            name_ifc2nd = "FORCE_CONSTANTS",
                            list_taufactor = [1.0]
                            )

row = []
with mp.Pool(processes=64, maxtasksperchild=1) as pool, tqdm(total=len(temp_range), desc="Saving thermal conductivity data" ) as pbar:
    for _, result in enumerate(pool.imap_unordered(condt_temp, temp_range, chunksize=1)):
        pbar.update()
        row.append(result)
