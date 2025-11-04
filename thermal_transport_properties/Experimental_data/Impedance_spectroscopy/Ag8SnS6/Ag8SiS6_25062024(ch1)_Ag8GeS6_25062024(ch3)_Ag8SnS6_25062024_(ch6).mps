EC-LAB SETTING FILE

Number of linked techniques : 4
Grouped channel(s) : 1, 3, 6

EC-LAB for windows v11.52 (software)
Internet server v11.52 (firmware)
Command interpretor v11.52 (firmware)

Filename : C:\Users\AKZMESS001\Documents\EC-Lab\Data\Anupama\Ag8SiS6_25062024(ch1)_Ag8GeS6_25062024(ch3)_Ag8SnS6_25062024_(ch6)\Ag8SiS6_25062024(ch1)_Ag8GeS6_25062024(ch3)_Ag8SnS6_25062024_(ch6).mps

Device : VSP-300
Electrode connection : standard
Potential control : Ewe
Ewe ctrl range : min = -10,00 V, max = 10,00 V
Ewe,I filtering : <None>
Safety Limits :
	Do not start on E overload
Channel : Grounded
Electrode material : 
Initial state : 
Electrolyte : 
Comments : Na1_25 10_C01
Comments : Na1_75 10_C02
Comments : Na1_5 11_C01
Cable : standard
Electrode surface area : 0,001 cm²
Characteristic mass : 0,001 g
Equivalent Weight : 0,000 g/eq.
Density : 0,000 g/cm3
Volume (V) : 0,001 cm³
Record Ece
Record Power
Record EIS quality indicators
Cycle Definition : Charge/Discharge alternance
TCU : 
   Chamber name : REGULATOR BINDER-MB2
   Delta time : 60
   Delta temperature : 0,100000001490116
   Has end temperature : True
   End temperature : 25
Turn to OCV between techniques

Technique : 1
Temperature Control Unit
Temp List (°C)      25/-40/-30/-20/-10/0.0/10/20/25/30/40/50/60
Stab. deltaT (°C)   0,5
Stab. dT/dt (°C/mn) 0,2
Stab. tmax (s)      3600

Technique : 2
Wait
select              1                   
td (h:m:s)          1:00:0,0000         
from                0                   
tech. num.          1                   
date (m/d/y)        07/03/09            
date (h:m:s)        11:40:55            
record              0                   
dE (mV)             0,00                
dI                  0,000               
unit dI             µA                  
dt (s)              0,0000              

Technique : 3
Potentio Electrochemical Impedance Spectroscopy
Mode                Single sine         
E (V)               0,0000              
vs.                 Ref                 
tE (h:m:s)          0:00:0,0000         
record              0                   
dI                  0,000               
unit dI             mA                  
dt (s)              0,000               
fi                  5,000               
unit fi             MHz                 
ff                  1,000               
unit ff             Hz                  
Nd                  15                  
Points              per decade          
spacing             Logarithmic         
Va (mV)             10,0                
pw                  0,40                
Na                  3                   
corr                0                   
E range min (V)     -10,000             
E range max (V)     10,000              
I Range             Auto                
Bandwidth           8                   
nc cycles           0                   
goto Ns'            0                   
nr cycles           0                   
inc. cycle          0                   

Technique : 4
Loop
goto Ne             1                   
nt times            12                  
