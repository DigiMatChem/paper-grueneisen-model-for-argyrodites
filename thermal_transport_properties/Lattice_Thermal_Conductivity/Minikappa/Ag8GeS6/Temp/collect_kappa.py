import os
import pandas as pd
import matplotlib.pyplot as plt

# Get the current working directory and project name for CSV output
cwd = os.getcwd()
path_parts = cwd.strip(os.sep).split(os.sep)
last_two = os.path.join(path_parts[-2], path_parts[-1]) 

# Get all kappa files and extract temperatures
file_names = [f for f in os.listdir() if f.startswith('minikappa-') and f.endswith('.dat')]
temperatures = sorted({int(f.split('-')[1]) for f in file_names})

# Initialize DataFrame
columns = ['phxx', 'phyy', 'phzz', 'diffx', 'diffy', 'diffz', 'Totalxx', 'Totalyy', 'Totalzz']
kappa = pd.DataFrame(index=temperatures, columns=columns)

# Fill the DataFrame
for file in file_names:
    try:
        temp = int(file.split('-')[1])
        with open(file, 'r') as f:
            lines = f.readlines()
            ph = lines[0].strip().split()
            diff = lines[1].strip().split()
            total = lines[-1].strip().split()
            
            kappa.loc[temp, 'phxx'] = float(ph[0])
            kappa.loc[temp, 'phyy'] = float(ph[4])
            kappa.loc[temp, 'phzz'] = float(ph[-1])
            
            kappa.loc[temp, 'diffx'] = float(diff[0])
            kappa.loc[temp, 'diffy'] = float(diff[4])
            kappa.loc[temp, 'diffz'] = float(diff[-1])
            
            kappa.loc[temp, 'Totalxx'] = float(total[0])
            kappa.loc[temp, 'Totalyy'] = float(total[4])
            kappa.loc[temp, 'Totalzz'] = float(total[-1])
    except Exception as e:
        print(f"Skipping file {file} due to error: {e}")

# Compute averages
kappa['ph'] = kappa[['phxx', 'phyy', 'phzz']].mean(axis=1)
kappa['diff'] = kappa[['diffx', 'diffy', 'diffz']].mean(axis=1)
kappa['Total'] = kappa[['Totalxx', 'Totalyy', 'Totalzz']].mean(axis=1)

# Save to CSV
filename_csv = f'kappa_{last_two.replace("/", "_")}.csv'
kappa.to_csv(filename_csv)

# Plotting
plt.plot(kappa.index, kappa['Total'], ':', markerfacecolor="white", color='blue', label='Total')
plt.plot(kappa.index, kappa['diff'], ':', markerfacecolor="white", color='red', label='Diffuson')
plt.plot(kappa.index, kappa['ph'], ':', markerfacecolor="white", color='black', label='Phonon')
plt.xlabel('Temperature (K)')
plt.ylabel('Thermal Conductivity (W/mK)')
plt.legend()
plt.tight_layout()
filename_pdf = f'kappa_{last_two.replace("/", "_")}.pdf'
plt.savefig(filename_pdf)
plt.close()
