import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import os

neutrino_energy = 100

# fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_{}TeV377_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/had_{}TeV411_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_{}TeV371_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)

# fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)
fname5 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)


plt.figure(figsize=(10,7))

PMT_ID = [1, 1, 6]
# layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(fname1, PMT_ID)
layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(fname2, PMT_ID)
layerID3, columnID3, cellID3, time3, x3, y3, z3, energy3 = get_data(fname3, PMT_ID)
# layerID4, columnID4, cellID4, time4, x4, y4, z4, energy4 = get_data(fname4, PMT_ID)
layerID5, columnID5, cellID5, time5, x5, y5, z5, energy5 = get_data(fname5, PMT_ID)

PMT_pos = PMT_ID_to_pos(PMT_ID)
dx, dy, dz = seperation_vector('A', PMT_pos)
d = distance_to_vertex(dx, dy, dz)
theta = angle_to_vertex('A', dx, dy, dz)

# time1 = clean_data(time1, d)
time2 = clean_data(time2, d)
time3 = clean_data(time3, d)
# time4 = clean_data(time4, d)
time5 = clean_data(time5, d)

bin_size = 1

bins = []
# time = [time1, time2, time4]
time = [time3, time2, time5]
for j in range(3):
    bin = int((max(time[j]) - min(time[j])) / bin_size)
    bins.append(bin)

# plt.hist(time1, bins=bins[0], range=(min(time1), max(time1)), histtype='step', color='C0', label='tau')
# plt.hist(time2, bins=bins[1], range=(min(time2), max(time2)), histtype='step', color='C1', label='had')
# plt.hist(time4, bins=bins[2], range=(min(time4), max(time4)), histtype='step', color='C2', label='Sum')
# plt.legend(fontsize=12)
# plt.xlim(min(time4), min(time4)+60)
# plt.title('Sensor Distance = {:.2f} m'.format(d), fontsize=14)
# plt.xlabel('Time [ns]')
# plt.ylabel('Photon Count')


plt.hist(time3, bins=bins[0], range=(min(time3), max(time3)), histtype='step', color='C0', label='e-')
plt.hist(time2, bins=bins[1], range=(min(time2), max(time2)), histtype='step', color='C1', label='had')
plt.hist(time5, bins=bins[2], range=(min(time5), max(time5)), histtype='step', color='C2', label='Sum')
plt.legend(fontsize=12)
plt.xlim(min(time5), min(time5)+60)
plt.title('Sensor Distance = {:.2f} m'.format(d), fontsize=14)
plt.xlabel('Time [ns]')
plt.ylabel('Photon Count')

plt.show()
