"""
This script creates a single nice plot of the time residuals. The user
will need to choose whether it plots the tau or e- event. Comment out
everything numbered 3 and 4 to plot tau events and comment out 1 and 4
to plot e- events.
"""

import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import os

energy = 100

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_{}TeV377_nt_Ntuple.csv'.format(energy, energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/had_{}TeV411_nt_Ntuple.csv'.format(energy, energy)
# fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_{}TeV371_nt_Ntuple.csv'.format(energy, energy)
fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_had_merge_{}_TeV.csv'.format(energy, energy)
# fname5 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_had_merge_{}_TeV.csv'.format(energy, energy)

data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')
# data3 = np.loadtxt(fname3, delimiter=',', comments='#')
data4 = np.loadtxt(fname4, delimiter=',', comments='#')
# data5 = np.loadtxt(fname5, delimiter=',', comments='#')

plt.figure(figsize=(10,7))

PMT_ID = [0, 1, 5]
layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(data1, PMT_ID)
layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(data2, PMT_ID)
# layerID3, columnID3, cellID3, time3, x3, y3, z3, energy3 = get_data(data3, PMT_ID)
layerID4, columnID4, cellID4, time4, x4, y4, z4, energy4 = get_data(data4, PMT_ID)
# layerID5, columnID5, cellID5, time5, x5, y5, z5, energy5 = get_data(data5, PMT_ID)

# print(len(time1), len(time2))

PMT_pos = PMT_ID_to_pos(PMT_ID)
dx, dy, dz = seperation_vector(PMT_pos)
d = distance_to_vertex(dx, dy, dz)
theta = angle_to_vertex(dx, dy, dz)

N = min(len(time1), len(time2))
# print(N)
# N = min(len(time3), len(time2))

# time_range = get_range(time4)
# print(time_range)

bin_size = 0.05
# print(bin_size)
bins = int((max(time4) - min(time4))/bin_size)
range = (min(time4), max(time4))

n1, bins1, patches1 = plt.hist(time1, bins=bins, range=range, histtype='step', color='C0', label='tau')
n2, bins2, patches2 = plt.hist(time2, bins=bins, range=range, histtype='step', color='C1', label='had')
n4, bins4, patches4 = plt.hist(time4, bins=bins, range=range, histtype='step', color='C2', label='Sum')
plt.legend(fontsize=12)
tmin, tmax = get_range_plot(n4, bins4)
plt.xlim(tmin, tmax)
plt.title('Sensor Distance = {:.2f} m'.format(d), fontsize=14)
plt.xlabel('Time [ns]')
plt.ylabel('Photon Count')
# plt.yscale('log')

# bins = int((max(time5) - min(time5))/bin_size)
# range = (min(time5), max(time5))
#
# n3, bins3, patches3 = plt.hist(time3, bins=bins, range=range, histtype='step', color='C0', label='e-')
# n2, bins2, patches2 = plt.hist(time2, bins=bins, range=range, histtype='step', color='C1', label='had')
# n5, bins5, patches5 = plt.hist(time5, bins=bins, range=range, histtype='step', color='C2', label='Sum')
# plt.legend(fontsize=12)
# tmin, tmax = get_range_plot(n5, bins5)
# plt.xlim(tmin, tmax)
# plt.title('Sensor Distance = {:.2f} m'.format(d), fontsize=14)
# plt.xlabel('Time [ns]')
# plt.ylabel('Photon Count')
# plt.yscale('log')

plt.show()
