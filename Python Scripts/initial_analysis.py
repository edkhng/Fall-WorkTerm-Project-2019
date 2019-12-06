"""
Creates a 3D plot showing detector and the hits detected.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/50 TeV/tau_50TeV23_nt_Ntuple.csv'

data = np.loadtxt(fname, delimiter=',', comments='#')
layerID = data[:,0]
columnID = data[:,1]
cellID = data[:,2]
time = data[:,3]/500  # [ns]
x = data[:,4]/1000  # [m]
y = data[:,5]/1000  # [m]
z = data[:,6]/1000  # [m]
energy = data[:,7]

colors = 1/(time)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y, z, c=colors, cmap='plasma')
plt.show()
