import numpy as np
import scipy.constants as const
from pymatgen.io.phonopy import get_gruneisenparameter
from pymatgen.core.units import amu_to_kg
import math
import h5py
import pandas as pd
import matplotlib.pyplot as plt

from Phono3pyPowerTools.Phono3pyIO import Phono3pyKappaHDF5

#with Phono3pyKappaHDF5('kappa-m61410.hdf5') as kappa_hdf5:
#        freqs = kappa_hdf5.GetModeFreqs()
#        tau = kappa_hdf5.GetModeTau(temp=300)


#frequency= freqs.flatten()
#tau = tau.flatten()

pi = np.pi
gru = get_gruneisenparameter(gruneisen_path="gruneisen.yaml")
volpc = (gru.structure.volume*(1E-30))
debye_freq = gru.debye_temp_phonopy()* const.value("Boltzmann constant in Hz/K") 
debye_freq = debye_freq*2*pi
freq=gru.frequencies*const.tera*2*pi

#print('Debye_freq', debye_freq)

debye_omega = (gru.debye_temp_phonopy() * const.k)/const.hbar
#print('Debye_omega', debye_omega)

freq = gru.frequencies*2*np.pi*const.tera

def life_time(theta_d=None, t=300):
	if gru.structure is None:
		raise ValueError("Structure is not defined.")

	average_mass = np.mean([s.specie.atomic_mass for s in gru.structure]) * amu_to_kg
	#print(average_mass)
	if theta_d is None:
		theta_d = gru.acustic_debye_temp
	mean_g = gru.average_gruneisen(t=t, squared=True, limit_frequencies='acoustic')
	v = (debye_omega)/ ((((6*np.pi ** 2) * (gru.structure.num_sites/volpc))) **(1/3)) #in m/s
	print('Sound velocity', v)
	p1 = (((1 - 0.514 * mean_g ** -1) + (0.228 * mean_g ** -2)) / 0.0948)
	p2 = ((const.hbar ** 2) * (mean_g ** 2)) / (const.k * theta_d * average_mass * v * ((volpc**(1/3))) )
	p= p1*p2
	tau = ((p * (freq)**2 * (t / theta_d)) *math.exp (-(theta_d) / (3 * t))) #Lifetime in s or
	tau_ps = (tau/const.tera)  #picosecondas
	#print('tau', tau)
	return tau_ps

# get lifetimes at different temps
lifetimes = life_time(theta_d=gru.acoustic_debye_temp, t=300)
# lifetimes with gru multiplicity
lifetime_mat = []
freq_mat = []
for inx, qgru in enumerate(lifetimes):
	life_arr = []
	fre_arr = []
	for i in range(len(qgru)):
		life_arr.extend([qgru[i]] * gru.multiplicities[i])
		fre_arr.extend([gru.frequencies[inx][i]]*gru.multiplicities[i])

	lifetime_mat.append(life_arr)
	freq_mat.append(fre_arr)

xs = 1/np.array(lifetime_mat).flatten()
ys = np.array(freq_mat).flatten()

MACE= pd.read_csv('MACE-3b0.csv')
AM_tau=pd.read_csv('AM_tau.csv')
AM_tau_subtracted_pd=pd.read_csv('AM_tau_sub_pd.csv')


plt.scatter(ys, xs, label='GM')
plt.scatter(MACE['Frequency'], MACE['Lifetime'], label='mace-mp03b')
plt.scatter(AM_tau['omega'], AM_tau['tau'], label='AM')
plt.scatter(AM_tau_subtracted_pd['omega'], AM_tau_subtracted_pd['tau'], label='AM_sub_pd')
plt.legend()
plt.ylim(0,8)
plt.ylabel(rf"$\mathrm{{(Lifetimes)}}$")
plt.xlabel(r"$\mathrm{Frequencies}$")
plt.show()
plt.savefig('Tau_comp.tif')

