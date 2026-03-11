# -*- coding: utf-8 -*-
"""
Created on Wed May 21 14:38:36 2025

@author: SAINATH
"""
import pandas as pd
import numpy as np
from numpy import array
import math

#import the output file from OTSun simulation with ASTMG173
df=pd.read_csv("/..../Receiver_Diameter/1740_70_V_G_CO_CSR.csv") #write the filepath
astmg=pd.read_csv("/..../ASTMG173.csv")
Er = df.iloc[:,:-10].values
Xr = df.iloc[:,1:-9].values
yr = df.iloc[:,2:-8].values
Zr = df.iloc[:,3:-7].values
W = df.iloc[:,10:].values
C1 = astmg.iloc[:,:-1].values
C2 = astmg.iloc[:,1:].values

N= astmg['C1'].count()
nl = df['Column1'].count()

# Compute Eo

Eo = np.zeros((nl, 1))
for j in range(N):
    for i in range(nl):
        if W[i] == C1[j]:
            Eo[i] = C2[j].item()

# Compute Eray and Wray

Eray = (59.9332 * Eo) / 3000 # Here put area of solar plane and no. of rays per lambda. Both the figures are in OTSun results file. 
Wray = Eray * Er

# Compute theta
Zrmin = (df['Column4'].min())-0.001
theta = np.arctan(Xr / (Zr - Zrmin)) * (180 / np.pi)

# Compute th
th = np.where(theta >= 0, (90 - theta) * 2, -(90 + theta) * 2)

# Sorting and table creation
ddf = pd.DataFrame({'th': th.ravel(), 'Wray': Wray.ravel()}).sort_values(by='th')
S1 = ddf['th'].values
S2 = ddf['Wray'].values

# Compute Flux
R =  (df['Column2'].max())/1000  #Receiver radius
stp = 5
Ang = np.arange(-180, 181, stp)
nf = len(Ang)
As = (math.pi / 180) * stp * R * 10

Flux = np.zeros(nf)
FX = np.sum(S2)
for j in range(nf):
    for i in range(nl):
        if S1[i] <= Ang[j]:
            Flux[j] = np.sum(S2[:i+1])

# Compute F1
F1 = np.zeros(nf-1)
for j in range(0, nf-1):
    F1[j] = (Flux[j+1] - Flux[j]) / As

# Plot results
import matplotlib.pyplot as plt
Ang1 = np.arange(-180+(stp/2), 180+(stp/2), stp)
plt.plot(Ang1, F1)
plt.xlabel('Receiver Angle (°)')
plt.ylabel('Flux (W/m²)')
plt.show()

# Write results to file
with open('Data_PTC_1740_70_V_G_CO_CSR', 'w') as file:


    file.write("\n\nAngle Flux(W/m^2)\n")
    for i in range(len(Ang)):
        file.write(f"{Ang1[i]} {F1[i]:.4f}\n")

    file.write(f"\nTotal Energy Flux (W/m^2):{FX:.4f}\n")

