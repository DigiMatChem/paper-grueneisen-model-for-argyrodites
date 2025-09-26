import numpy as np
import math
import cmath
import os
import sys
import scipy.constants as const
from pymatgen.io.phonopy import get_gruneisenparameter
from pymatgen.core.units import amu_to_kg

#from gru_lifetimes import life_time

"""
author: @yixia
A python class with a method to compute kL_min
The calcualtion of kL_min relies on modified "group_velocity.py" for Phonopy version of 2.7.1
"""

gru = get_gruneisenparameter(gruneisen_path="gruneisen.yaml")
debye_freq = gru.debye_temp_phonopy()* const.value("Boltzmann constant in Hz/K") 
debye_freq = debye_freq*2*np.pi
print('Debye_freq', debye_freq)

class class_kappa:
    """
    input: obj_poscar
    """
    def __init__(self):
        pass

    def get_minikappa_phonopy(self,
                          mesh_in = [8,8,8],
                          sc_mat = np.eye(3)*2,
                          pm_mat = np.eye(3),
                          list_temp = [300.0],
                          name_pcell = "POSCAR-unitcell",
                          name_ifc2nd = "FORCE_CONSTANTS",
                          list_taufactor = [2.0]):
        try:
            import phonopy
            from phonopy import Phonopy
            from phonopy.structure.atoms import PhonopyAtoms
            from phonopy.interface.calculator import read_crystal_structure
        except:
            print ("Phonopy API version of 2.7.1 is required!")
            
        # load phonopy object
        phonon = phonopy.load(supercell_matrix = sc_mat,
                              primitive_matrix = pm_mat,
                              unitcell_filename = name_pcell,
                              is_symmetry = False,
                              force_constants_filename = name_ifc2nd)
        primcell = (phonon.get_primitive()).cell.T
        volpc = np.abs(np.dot(np.cross(primcell[1],primcell[2]),primcell[0]))
        print (f"Number of atoms in primitive cell: {len(phonon.primitive.masses)}")

        # phonon
        phonon.run_mesh(mesh_in,
                        with_eigenvectors = False,
                        is_gamma_center = True,
                        with_group_velocities = True,
                        is_time_reversal = False,
                        is_mesh_symmetry = False)
        mesh_dict = phonon.get_mesh_dict()
        qpoints = mesh_dict['qpoints']
        weights = mesh_dict['weights']
        freqs = mesh_dict['frequencies']
        eigs = mesh_dict['eigenvectors']
        gvfull = mesh_dict['group_velocities']

        # unit
        hbar = const.hbar #possible J.s
        kB = const.k #J.K-1
        pi = np.pi
        volpc = volpc*(1E-30) #A³ to m³
        gvfull = gvfull*100 # m/s
        freqs = (freqs*const.tera)*pi*2
        nband = len(freqs[0])
        freqcf = 0.1

        # kappa
        for temp in list_temp:
            def life_time(theta_d=None, t=None):
                average_mass = np.mean([s.specie.atomic_mass for s in gru.structure]) * amu_to_kg
                #print(average_mass)
                if theta_d is None:
                    theta_d = gru.acustic_debye_temp
                mean_g = gru.average_gruneisen(t=temp, squared=True)
                #v=3320 #Experimental sound velocity
                v = (debye_freq)/ ((((6*np.pi ** 2) * (gru.structure.num_sites/volpc))) **(1/3)) #in m/s
                print('Sound velocity', v)
                p1 = (((1 - 0.514 * mean_g ** -1) + (0.228 * mean_g ** -2)) / 0.0948)
                p2 = ((const.hbar ** 2) * (mean_g ** 2)) / (const.k * theta_d * average_mass * v * ((volpc**(1/3))) )
                p= p1*p2
                tau_s = (p * (gru.frequencies*const.tera*2*pi)**2 * (t / theta_d) * math.exp(-theta_d / (3 * t)))  #Lifetime in s or
                tau = (tau_s*const.tera)  #picoseconds
                #print('tau', tau)
                return tau_s
            # get lifetimes at different temps
            lifetimes = life_time(theta_d=gru.acoustic_debye_temp, t=temp)
            # lifetimes with gru multiplicity
            lifetime_mat = []
            for inx, qgru in enumerate(lifetimes):
                life_arr = []
                for i in range(len(qgru)):
                    life_arr.extend([qgru[i]] * gru.multiplicities[i])
                lifetime_mat.append(life_arr)

            life_times = np.array(lifetime_mat)  # use this lifetimes

            for factor in list_taufactor:
                kappaband=np.zeros((nband,nband,3,3), dtype=np.complex128, order='C')
                nqpt=len(qpoints)
                nband=len(freqs[0])
                for i in range(nband):
                    for j in range(nband):
                        for iq in range(nqpt):
                            for k in range(3):
                                for kp in range(3):
                                    omega1=freqs[iq,i]
                                    omega2=freqs[iq,j]
                                    if omega1>freqcf and omega2>freqcf:
                                        if (freqs[iq, i] / 2 / pi) > 0:
                                            Gamma1_w=freqs[iq,i]/2/pi*factor
                                            Gamma1 = (life_times[i][iq]*factor)*2*pi
                                            #print(Gamma1, Gamma1_w)
                                        else:
                                            Gamma1=1E10
                                        if (freqs[iq, j] / 2 / pi) > 0:
                                            Gamma2_w=freqs[iq,j]/2/pi*factor
                                            Gamma2 = (life_times[j][iq]*factor)*2*pi
                                            #print(Gamma2, Gamma2_w)
                                        else:
                                            Gamma2=1E10
                                        fBE1=1.0/(np.exp(hbar*omega1/kB/temp)-1.0)
                                        fBE2=1.0/(np.exp(hbar*omega2/kB/temp)-1.0)
                                        tmpv=(gvfull[iq,i,j,k]*gvfull[iq,j,i,kp]).real
                                        #print('vel',tmpv)
                                        kappaband[i,j,k,kp]=kappaband[i,j,k,kp]+(omega1+omega2)/2* \
                                            (fBE1*(fBE1+1)*omega1+fBE2*(fBE2+1)*omega2)*tmpv \
                                            /(4*(omega1-omega2)**2+(Gamma1+Gamma2)**2)*\
                                            (Gamma1+Gamma2)
                                        #print(Gamma1)
                # conversion
                print (f"Factor: {factor}")
                kappaband=kappaband*(hbar**2)/(kB*temp*temp*volpc*nqpt)

                kappaD  = np.zeros((3,3), dtype=np.complex128, order='C')
                kappaOD = np.zeros((3,3), dtype=np.complex128, order='C')
                kappaF  = np.zeros((3,3), dtype=np.complex128, order='C')
                for i in range(nband):
                    for j in range(nband):
                        kappaF=kappaF+kappaband[i,j]
                        if i==j:
                            kappaD=kappaD+kappaband[i,j]
                        else:
                            kappaOD=kappaOD+kappaband[i,j]
                f = open("kappa"+"-"+str(temp)+"-"+str(factor)+".dat", "w")
                f.write("   ".join(map(str, np.round(kappaD.real.reshape((9,1)).flatten(),decimals=8) ))+"\n")
                f.write("   ".join(map(str, np.round(kappaOD.real.reshape((9,1)).flatten(),decimals=8) ))+"\n")
                f.write("   ".join(map(str, np.round(kappaF.real.reshape((9,1)).flatten(),decimals=8) ))+"\n")
                f.close()
                if True:
                    print ("Diagonal part of thermal conductivity: ")
                    print (kappaD.real[0,0])
                    print ("Off-diagonal part of thermal conductivity: ")
                    print (kappaOD.real[0,0])
                    print ("Total thermal conductivity: ")
                    print (kappaF.real[0,0])
        #return [kappaD.real[0,0], kappaOD.real[0,0]]

