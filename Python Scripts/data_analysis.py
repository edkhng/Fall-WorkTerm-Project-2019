import numpy as np
import matplotlib.pyplot as plt
import analysis_functions as af

simID = 1
energy = 100
size = 150

#fname = '{} TeV/tau_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)
#fname = '{} TeV/e_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)
fname = '{} TeV/had_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)
#fname = '{}TeV_{}x{}/muon{}_nt_Ntuple.csv'.format(energy, size, size, simID)

#fname = '{} TeV/tau_had_{}.csv'.format(energy, simID)
#fname = '{} TeV/e_had_{}.csv'.format(energy, simID)

layerID, columnID, cellID, time, x, y, z, energy = af.get_data(fname)
print(len(x))

'''
print(len(data))

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(x,y,z, c=time**0.2, cmap='hsv')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()


select1 = columnID == 0
select2 = columnID == 1
select3 = columnID == 2
select4 = columnID == 3
select5 = columnID == 4
select6 = columnID == 5

time1 = time[select1]
x1 = x[select1]
y1 = y[select1]
z1 = z[select1]

time2 = time[select2]
x2 = x[select2]
y2 = y[select2]
z2 = z[select2]

time3 = time[select3]
x3 = x[select3]
y3 = y[select3]
z3 = z[select3]

time4 = time[select4]
x4 = x[select4]
y4 = y[select4]
z4 = z[select4]

time5 = time[select5]
x5 = x[select5]
y5 = y[select5]
z5 = z[select5]

time6 = time[select6]
x6 = x[select6]
y6 = y[select6]
z6 = z[select6]

plt.figure()
plt.plot(time1, x1, '.')
plt.plot(time2, x2, '.')
plt.plot(time3, x3, '.')
plt.plot(time4, x4, '.')
plt.plot(time5, x5, '.')
plt.plot(time6, x6, '.')
plt.ylabel('X Label')

plt.figure()
plt.plot(time1, y1, '.')
plt.plot(time2, y2, '.')
plt.plot(time3, y3, '.')
plt.plot(time4, y4, '.')
plt.plot(time5, y5, '.')
plt.plot(time6, y6, '.')
plt.ylabel('Y Label')

plt.figure()
plt.plot(time1, z1, '.')
plt.plot(time2, z2, '.')
plt.plot(time3, z3, '.')
plt.plot(time4, z4, '.')
plt.plot(time5, z5, '.')
plt.plot(time6, z6, '.')
plt.ylabel('Z Label')
#plt.show()


from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
plt.plot(x1, y1, z1, '.')
plt.plot(x2, y2, z2, '.')
plt.plot(x3, y3, z3, '.')
plt.plot(x4, y4, z4, '.')
plt.plot(x5, y5, z5, '.')
plt.plot(x6, y6, z6, '.')
#ax.scatter(x,y,z, c=time**0.2, cmap='hsv')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()


PMT_ID = [1, 1, 10]

layerID, columnID, cellID, time, x, y, z, energy = af.get_data(fname, PMT_ID)

print(PMT_ID)

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(x,y,z, c=time**0.2, cmap='hsv')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()


plt.figure()
plt.hist(time, bins=70, range=(75,200))
plt.show()
'''
