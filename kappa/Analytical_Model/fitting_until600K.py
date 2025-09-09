import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import hbar, k
from scipy.integrate import trapezoid as trapz
from scipy.optimize import curve_fit

# --- Functions ---

def HarmonicVs(vTrans, vLong):
    return (3 * (vLong**-3 + 2 * vTrans**-3) ** -1) ** (1/3)

def C_omega(angular_freq, Temp):
    x = (hbar * angular_freq) / (k * (Temp + 1e-10))
    return x**2 * np.exp(x) / (np.exp(x) - 1)**2

def omega_IR(n, vs, N):
    return (6 * n * np.pi**2)**(1/3) * vs * N**(-1/3)

def kappa_2ch(Temp, C1, A, C3, P):
    w_IR = omega_IR(n, vs, N)
    Prop_vals = np.where(omega < w_IR)
    Diff_vals = np.where(omega >= w_IR)
    kappa_tot = []

    for T in Temp:
        inverse_tau = (C1 * 1e-16 * omega[Prop_vals]**2 * T) + \
                      (A * omega[Prop_vals]) + \
                      (C3 * 1e-40 * omega[Prop_vals]**4)
        kappa_pr_omega = DOS[Prop_vals] * C_omega(omega[Prop_vals], T) * vs**2 / inverse_tau
        kappa_diff_omega = DOS[Diff_vals] * C_omega(omega[Diff_vals], T) * P * omega[Diff_vals]
        kappa_pr = (1 / 3) * (3 * n * k) * trapz(kappa_pr_omega, omega[Prop_vals])
        kappa_diff = (1 / np.pi) * n**(1/3) * k * trapz(kappa_diff_omega, omega[Diff_vals])
        kappa_tot.append(kappa_pr + kappa_diff)

    return kappa_tot

def both_kappas(Temp, C1, A, C3, P):
    omega_IR_val = omega_IR(n, vs, N)
    Prop_vals = np.where(omega < omega_IR_val)
    Diff_vals = np.where(omega >= omega_IR_val)

    kappa_tot, prop, diff = [], [], []
    for T in Temp:
        inverse_tau = (C1 * 1e-16 * omega[Prop_vals]**2 * T) + \
                      (A * omega[Prop_vals]) + \
                      (C3 * 1e-40 * omega[Prop_vals]**4)
        kappa_pr_omega = DOS[Prop_vals] * C_omega(omega[Prop_vals], T) * vs**2 / inverse_tau
        kappa_diff_omega = DOS[Diff_vals] * C_omega(omega[Diff_vals], T) * P * omega[Diff_vals]
        kappa_pr = (1 / 3) * (3 * n * k) * trapz(kappa_pr_omega, omega[Prop_vals])
        kappa_diff = (1 / np.pi) * n**(1/3) * k * trapz(kappa_diff_omega, omega[Diff_vals])
        prop.append(kappa_pr)
        diff.append(kappa_diff)
        kappa_tot.append(kappa_pr + kappa_diff)

    df = pd.DataFrame({'kappa_pr': prop, 'kappa_diff': diff, 'kappa_total': kappa_tot})
    df.to_csv('kappa_pd.csv', index=False)
    return prop, diff

#Fit Data
# dos model
dat = np.loadtxt('total_dos.dat')#,delimiter='\t',skiprows=1,dtype=float)
freq= dat[:,0] #THz unit
DOS =dat[:,1] # not normalized
omega= 2*np.pi*freq*1E12 #Hz
print(omega.shape)
intDOS = trapz(DOS,omega)
DOS = DOS/intDOS # normalized DOS
# material props
vl = 2871
vt = 1427
vs = HarmonicVs(vt,vl)
print("Average SoS = ", vs) # m/s harmonic mean speed of sound
N = 60 # atoms per primitive unit cell
n = 4.64E28 #atoms/m^3
ThetaD = (hbar/k) * (6*np.pi**2*n)**(1/3) * vs
print("Theta= ", (hbar/k) * (6*np.pi**2*n)**(1/3) * vs)

# load data
data = pd.read_csv('sample_temp_conductivity.csv')
#data = pd.read_csv('Our_data.csv')  # Replace with your file path
#T=data.index
T = data['Tem'].values  # Temperature data (in K)
kappa_exp = data['Conductivity'].values  # Experimental thermal conductivity data

# --- Fit only up to 300 K ---
mask_fit = T <= 300
T_fit = T[mask_fit]
kappa_fit = kappa_exp[mask_fit]

val, var = curve_fit(kappa_2ch, T_fit, kappa_fit, p0=[1.60, 0.02, 0.5, 1],
                     bounds=((0, 0.015, 0, 0), (np.inf, 0.1, np.inf, 1)))

stdev = np.sqrt(np.diag(var))

# --- Print fit results ---
print('\nFitted Parameters:')
print(f'Ioffe-Regel Frequency = {omega_IR(n, vs, N) * 1e-12 / (2 * np.pi):.2f} THz')
print(f'C1 = {val[0]:.3f} ± {stdev[0]:.3f}  × 10⁻¹⁶')
print(f'A  = {val[1]:.3f} ± {stdev[1]:.3f}')
print(f'C3 = {val[2]:.3f} ± {stdev[2]:.3f}  × 10⁻⁴⁰')
print(f'P  = {val[3]:.3f} ± {stdev[3]:.3f}')

# --- Predict up to 600 K ---
T_predict = np.linspace(00, 600, 1000)
kappa_model = kappa_2ch(T_predict, *val)
kappa_pr, kappa_diff = both_kappas(T_predict, *val)

# --- Plot ---
plt.figure()
plt.plot(T, kappa_exp, 'ks', markerfacecolor='w', label='Experimental')
plt.plot(T_predict, kappa_model, 'r-', label='2ch Model (fit + extrapolated)')
plt.plot(T_predict, kappa_pr, '-', color='orange', label='Phonon')
plt.plot(T_predict, kappa_diff, '-', color='green', label='Diff')
plt.axvline(300, color='gray', linestyle='--', label='Fit limit (300 K)')
plt.xlabel('Temperature / K')
plt.ylabel('Thermal Conductivity / W m$^{-1}$ K$^{-1}$')
plt.legend()
plt.xlim([-10, 620])
plt.tight_layout()
plt.savefig('2ch_model_extrapolated_to_600K.tif', dpi=300)
plt.show()

# --- Optional: Save prediction ---
df_pred = pd.DataFrame({'T (K)': T_predict, 'kappa_total': kappa_model,
                        'kappa_pr': kappa_pr, 'kappa_diff': kappa_diff})
df_pred.to_csv('predicted_kappa_up_to_600K.csv', index=False)

# --- Plot DOS ---
plt.figure()
plt.plot(omega * 1e-12 / (2 * np.pi), DOS, 'k-')  # omega in THz
plt.axvline(omega_IR(n, vs, N) * 1e-12 / (2 * np.pi), linestyle=':', color='k', label=r'$\omega_{IR}$')
plt.xlabel(r'$\omega$ (THz)')
plt.ylabel('Normalized DOS')
plt.legend()
plt.tight_layout()
plt.savefig('DOS_plot_with_IR_cutoff.tif', dpi=300)
plt.show()

