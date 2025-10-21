import os
import pandas as pd
import matplotlib.pyplot as plt

cwd = os.getcwd()

path_parts = cwd.strip(os.sep).split(os.sep)
last_two = os.path.join(path_parts[-2], path_parts[-1]) 

temps = list(range(0, 620, 20))
columns = ['phxx', 'phyy', 'phzz', 'diffx', 'diffy', 'diffz', 'Totalxx', 'Totalyy', 'Totalzz']
Ge = pd.DataFrame(index=temps, columns=columns)

file_names = [f for f in os.listdir() if f.startswith('kappa-')]

for file in file_names:
    temp = int(file.split('-')[1])
    with open(file, 'r') as f:
        lines = f.readlines()
        ph = lines[0].strip().split()
        diff = lines[1].strip().split()
        total = lines[-1].strip().split()
        
        Ge.loc[temp, 'phxx'] = float(ph[0])
        Ge.loc[temp, 'phyy'] = float(ph[4])
        Ge.loc[temp, 'phzz'] = float(ph[-1])
        
        Ge.loc[temp, 'diffx'] = float(diff[0])
        Ge.loc[temp, 'diffy'] = float(diff[4])
        Ge.loc[temp, 'diffz'] = float(diff[-1])
        
        Ge.loc[temp, 'Totalxx'] = float(total[0])
        Ge.loc[temp, 'Totalyy'] = float(total[4])
        Ge.loc[temp, 'Totalzz'] = float(total[-1])

# Compute averages
Ge['ph'] = Ge[['phxx', 'phyy', 'phzz']].mean(axis=1)
Ge['diff'] = Ge[['diffx', 'diffy', 'diffz']].mean(axis=1)
Ge['Total'] = Ge[['Totalxx', 'Totalyy', 'Totalzz']].mean(axis=1)

# Save to CSV
#Ge.to_csv('kappa_Ge.csv')
csv_filename = f'kappa_{last_two.replace("/", "_")}.csv'
Ge.to_csv(csv_filename)

# Plotting
plt.plot(Ge.index, Ge['Total'], ':', marker='o', markerfacecolor="white", color='blue', label='Total')
plt.plot(Ge.index, Ge['diff'], ':', marker='o', markerfacecolor="white", color='red', label='Diffuson')
plt.plot(Ge.index, Ge['ph'], ':', marker='o', markerfacecolor="white", color='black', label='Phonon')
plt.xlabel('Temperature (K)')
plt.ylabel('Thermal Conductivity (W/mK)')
plt.legend()
plt.tight_layout()
plt.savefig('kappa.pdf')
plt.close()
