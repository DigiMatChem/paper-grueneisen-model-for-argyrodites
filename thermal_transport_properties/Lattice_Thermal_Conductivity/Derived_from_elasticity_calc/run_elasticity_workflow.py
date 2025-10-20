from termios import PENDIN

from pymatgen.core import Structure
from atomate2.vasp.powerups import update_user_incar_settings
from atomate2.vasp.flows.elastic import ElasticMaker
from jobflow_remote import submit_flow

struct = Structure.from_file("pth/to/the/structure") #load you structure


elastic_flow = ElasticMaker(bulk_relax_maker=None).make(structure=struct) #initialize your workflow. Modify settings if needed

new_flow = update_user_incar_settings(elastic_flow, {'ALGO': 'Normal', 'MAGMOM':None,'ISPIN':1,'ENCUT':520, 'GGA' : 'PE',
                                                    'ENAUG':None, 'SIGMA' : 0.05, 'NPAR' :12}) 

resources = {"nodes": 4, "ntasks": 192, "time": "08:00:00",
             "qverbatim": "#SBATCH --get-user-env"}   

# change worker and project name to what you setup in the jobflow remote config file.
print(submit_flow(new_flow, worker="phonon_worker", resources=resources, project="phonon"))