import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import timeit

neutrino_energy = 100

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/fitting_params_gaussian_tau_good.csv'
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/fitting_params_gaussian_e_good.csv'
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/fitting_params_bi_gaussian_tau_good.csv'
fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/fitting_params_bi_gaussian_e_good.csv'

data1 = np.loadtxt(fname1, delimiter=',', comments='#', usecols=(0,2,3,4,5,6,7,8,9,10))
data2 = np.loadtxt(fname2, delimiter=',', comments='#', usecols=(0,2,3,4,5,6,7,8,9,10))
data3 = np.loadtxt(fname3, delimiter=',', comments='#', usecols=(0,2,3,4,5,6,7,8,9,10,11,12))
data4 = np.loadtxt(fname4, delimiter=',', comments='#', usecols=(0,2,3,4,5,6,7,8,9,10,11,12))

energy1, pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, chisq1, bins1, params1 = data1[:,0], data1[:,1], data1[:,2], data1[:,3], data1[:,4], data1[:,5], data1[:,6], data1[:,7], data1[:,8], data1[:,9]
energy2, pos2a, pos2b, wid2a, wid2b, amp2a, amp2b, chisq2, bins2, params2 = data1[:,0], data2[:,1], data2[:,2], data2[:,3], data2[:,4], data2[:,5], data2[:,6], data2[:,7], data2[:,8], data2[:,9]
energy3, pos3a, pos3b, wid3a, wid3b, amp3a, amp3b, r3a, r3b, chisq3, bins3, params3 = data3[:,0], data3[:,1], data3[:,2], data3[:,3], data3[:,4], data3[:,5], data3[:,6], data3[:,7], data3[:,8], data3[:,9], data3[:,10], data3[:,11]
energy4, pos4a, pos4b, wid4a, wid4b, amp4a, amp4b, r4a, r4b, chisq4, bins4, params4 = data4[:,0], data4[:,1], data4[:,2], data4[:,3], data4[:,4], data4[:,5], data4[:,6], data4[:,7], data4[:,8], data4[:,9], data4[:,10], data4[:,11]

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

# print(energy1[:5])
# print(pos1a[:5])
# print(pos1b[:5])

plt.figure(figsize=(10,7))
plt.scatter(time_diff3, wid_diff3, label='tau events')
plt.scatter(time_diff4, wid_diff4, label='e- events')
plt.xlabel('Time difference [ns]')
plt.ylabel('FWHM difference [ns]')
plt.legend()

plt.figure(figsize=(10,7))
plt.scatter(time_diff3, amp_diff3, label='tau events')
plt.scatter(time_diff4, amp_diff4, label='e- events')
plt.xlabel('Time difference [ns]')
plt.ylabel('Amp difference')
plt.legend()

plt.figure(figsize=(10,7))
plt.scatter(amp_diff3, wid_diff3, label='tau events')
plt.scatter(amp_diff4, wid_diff4, label='e- events')
plt.xlabel('Amp difference')
plt.ylabel('FWHM difference [ns]')
plt.legend()

plt.show()
