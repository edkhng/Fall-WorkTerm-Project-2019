import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import os

neutrino_energy = 100

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_{}TeV377_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/had_{}TeV411_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
# fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_{}TeV371_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)

fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)
# fname5 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)


plt.figure(figsize=(10,7))

PMT_ID = [0, 1, 5]
layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(fname1, PMT_ID)
layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(fname2, PMT_ID)
# layerID3, columnID3, cellID3, time3, x3, y3, z3, energy3 = get_data(fname3, PMT_ID)
layerID4, columnID4, cellID4, time4, x4, y4, z4, energy4 = get_data(fname4, PMT_ID)
# layerID5, columnID5, cellID5, time5, x5, y5, z5, energy5 = get_data(fname5, PMT_ID)

# print(len(time1), len(time2))

PMT_pos = PMT_ID_to_pos(PMT_ID)
dx, dy, dz = seperation_vector('A', PMT_pos)
d = distance_to_vertex(dx, dy, dz)
theta = angle_to_vertex('A', dx, dy, dz)

# time1 = clean_data(time1, d)
# time2 = clean_data(time2, d)
# time3 = clean_data(time3, d)
# time4 = clean_data(time4, d)
# time5 = clean_data(time5, d)

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
