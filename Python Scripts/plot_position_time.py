import numpy as np
import matplotlib.pyplot as plt
import analysis_functions as af

simID = 1
energy = 60
size = 150

#fname = '{}TeV_{}x{}/tau{}_nt_Ntuple.csv'.format(energy, size, size, simID)
#fname = '{}TeV_{}x{}/e{}_nt_Ntuple.csv'.format(energy, size, size, simID)
#fname = '{}TeV_{}x{}/had{}_nt_Ntuple.csv'.format(energy, size, size, simID)
#fname = '{}TeV_{}x{}/muon{}_nt_Ntuple.csv'.format(energy, size, size, simID)

fname = '{}TeV_{}x{}/tau_had_{}_merge.csv'.format(energy, size, size, simID)
#fname = '{}TeV_{}x{}/e_had_{}.csv'.format(energy, size, size, simID)

layerID, columnID, cellID, time, x, y, z, energy = af.get_data(fname)

#print(len(x))

# select1 = columnID == 0
# select2 = columnID == 1
# select3 = columnID == 2
select1 = (columnID == 0) & (layerID == 2)
select2 = (columnID == 1) & (layerID == 2)
select3 = (columnID == 2) & (layerID == 2)
# select4 = columnID == 3
# select5 = columnID == 4
# select6 = columnID == 5

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

# time4 = time[select4]
# x4 = x[select4]
# y4 = y[select4]
# z4 = z[select4]
#
# time5 = time[select5]
# x5 = x[select5]
# y5 = y[select5]
# z5 = z[select5]
#
# time6 = time[select6]
# x6 = x[select6]
# y6 = y[select6]
# z6 = z[select6]

plt.figure()
plt.plot(time1, x1, '.', label='y= -50 m')
plt.plot(time2, x2, '.', label='y= 0 m')
plt.plot(time3, x3, '.', label='y= 50 m')
# plt.plot(time4, x4, '.')
# plt.plot(time5, x5, '.')
# plt.plot(time6, x6, '.')
plt.xlim(0, 800)
plt.legend()
plt.ylabel('x (m)')
plt.xlabel('Time (ns)')
#
# #
# # plt.figure()
# # plt.plot(time1, y1, '.')
# # plt.plot(time2, y2, '.')
# # plt.plot(time3, y3, '.')
# # # plt.plot(time4, y4, '.')
# # # plt.plot(time5, y5, '.')
# # # plt.plot(time6, y6, '.')
# # plt.ylabel('Y Label')
#
plt.figure()
plt.plot(time1, z1, '.', label='y= -50 m')
plt.plot(time2, z2, '.', label='y= 0 m')
plt.plot(time3, z3, '.', label='y= 50 m')
# plt.plot(time4, z4, '.')
# plt.plot(time5, z5, '.')
# plt.plot(time6, z6, '.')
plt.xlim(0, 800)
plt.ylabel('z (m)')
plt.xlabel('Time (ns)')
plt.legend()
#plt.show()

# x = np.array([-50, 0, 50])
# y = np.array([-50, 0, 50])
# z = np.linspace(-68.75, 68.75, 12)
# xx, yy, zz = np.meshgrid(x, y, z)


from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection = '3d')
plt.plot(time1, x1, z1, '.', label='y= -50 m')
plt.plot(time2, x2, z2, '.', label='y= 0 m')
plt.plot(time3, x3, z3, '.', label='y= 50 m')
# plt.plot(x1, y1, z1, '.')
# plt.plot(x2, y2, z2, '.')
# plt.plot(x3, y3, z3, '.')
# plt.plot([-20], [-10], [-62.5], 'x', c='r', ms=10, label='Event Vertex')
# plt.plot([-20, -20], [-10, -10], [-62.5, 62.5], 'r--', lw=1, label='Particle Trajectory')
# ax.scatter(xx, yy, zz, color='C0', marker='o', label='PMT')
#ax.scatter(x,y,z, c=time**0.2, cmap='hsv')
ax.set_xlim(0, 800)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('x (m)')
ax.set_zlabel('z (m)')
# ax.set_xlabel('x (m)')
# ax.set_ylabel('y (m)')
# ax.set_zlabel('z (m)')
plt.legend()
plt.show()
