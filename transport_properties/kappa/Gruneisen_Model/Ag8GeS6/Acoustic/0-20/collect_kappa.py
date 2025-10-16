import os
import pandas as pd
import matplotlib.pyplot as plt

Ge=pd.DataFrame(index= list(range(0, 620, 20)), columns= ['ph','diff','Total'])
file_names = [i for i in os.listdir() if i.startswith('kappa')]

for file in file_names:
    temp = int(file.split('-')[1:])
    with open(file, 'r') as f:
        ph = f.readlines()[0].strip('\n').split('  ')
    Ge.loc[temp,'phxx_pd'] = float(ph[0])
    Ge.loc[temp,'phyy_pd'] = float(ph[4])
    Ge.loc[temp,'phzz_pd'] = float(ph[-1])
    Ge['ph_pd'] = Ge[['phxx_pd', 'phyy_pd','phzz_pd']].mean(axis=1)
    
for file in file_names:
    temp = int(file.split('-')[1])
    with open(file, 'r') as f:
        diff = f.readlines()[1].strip('\n').split('  ')
    Ge.loc[temp,'diffx'] = float(diff[0])
    Ge.loc[temp,'diffy'] = float(diff[4])
    Ge.loc[temp,'diffz'] = float(diff[-1])
    Ge['diff'] = Ge[['diffx', 'diffy','diffz']].mean(axis=1)
    
for file in file_names:
    temp = int(file.split('-')[1])
    with open(file, 'r') as f:
        Total = f.readlines()[-1].strip('\n').split('  ')
    Ge.loc[temp,'Totalxx'] = float(Total[0])
    Ge.loc[temp,'Totalyy'] = float(Total[4])
    Ge.loc[temp,'Totalzz'] = float(Total[-1])
    Ge['Total'] = Ge[['Totalxx', 'Totalyy','Totalzz']].mean(axis=1)
    
    
    Ge.to_csv('boundary.csv')
    
plt.plot(Ge.index, Ge['Total'],':', markerfacecolor="white", color='blue', label='Total' )
plt.plot(Ge.index, Ge['diff'],':', markerfacecolor="white", color='red', label='diff' )
plt.plot(Ge.index, Ge['ph'],':', markerfacecolor="white", color='k', label='ph' )
