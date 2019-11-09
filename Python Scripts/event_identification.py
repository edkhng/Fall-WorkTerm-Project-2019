import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import timeit

neutrino_energy = 100

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_gaussian_tau_all.csv'.format(neutrino_energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_gaussian_e_all.csv'.format(neutrino_energy)
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_bi_gaussian_tau_all.csv'.format(neutrino_energy)
fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_bi_gaussian_e_all.csv'.format(neutrino_energy)

data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')
data3 = np.loadtxt(fname3, delimiter=',', comments='#')
data4 = np.loadtxt(fname4, delimiter=',', comments='#')

pos1a, pos1b, wid1a, wid1b, amp1a, amp1b = data1[:,0], data1[:,1], data1[:,2], data1[:,3], data1[:,4], data1[:,5]
pos2a, pos2b, wid2a, wid2b, amp2a, amp2b = data2[:,0], data2[:,1], data2[:,2], data2[:,3], data2[:,4], data2[:,5]
pos3a, pos3b, wid3a, wid3b, amp3a, amp3b, r3a, r3b = data3[:,0], data3[:,1], data3[:,2], data3[:,3], data3[:,4], data3[:,5], data3[:,6], data3[:,7]
pos4a, pos4b, wid4a, wid4b, amp4a, amp4b, r4a, r4b = data4[:,0], data4[:,1], data4[:,2], data4[:,3], data4[:,4], data4[:,5], data4[:,6], data4[:,7]

time_diff1 = pos1a - pos1b
time_diff2 = pos2a - pos2b
time_diff3 = pos3a - pos3b
time_diff4 = pos4a - pos4b

amp_diff1 = amp1a - amp1b
amp_diff2 = amp2a - amp2b
amp_diff3 = amp3a - amp3b
amp_diff4 = amp4a - amp4b

wid_diff1 = wid1a - wid1b
wid_diff2 = wid2a - wid2b
wid_diff3 = wid3a - wid3b
wid_diff4 = wid4a - wid4b

wid_ratio1 = wid1a/wid1b
wid_ratio2 = wid2a/wid2b
wid_ratio3 = wid3a/wid3b
wid_ratio4 = wid4a/wid4b
