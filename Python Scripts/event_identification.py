"""
This script analyzes the csv of the fitting parameters, chisq and p-value
produced by fit_data_e.py and fit_data_tau.py. Bad fits are filtered out
by removing fits with a time difference > 100, a amplitude ratio > 4 and
a FWHM ratio > 4. The analysis is seperated by the hit range and decay
time range. The distribution for each parameter is plotted to see if
there is seperation between tau and e- event for that variable.
"""

import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare

energy = 100

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_double_bi-gaussian_tau.csv'.format(energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_double_bi-gaussian_e.csv'.format(energy)
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_single_bi-gaussian_tau.csv'.format(energy)
fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_single_bi-gaussian_e.csv'.format(energy)

data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')
data3 = np.loadtxt(fname3, delimiter=',', comments='#')
data4 = np.loadtxt(fname4, delimiter=',', comments='#')

hits_range = (200000, 360000)
decay_time_range = (30, 60)

print("Total number of fits for each function")
print(np.shape(data1)[1])
print(np.shape(data2)[1])
print(np.shape(data3)[1])
print(np.shape(data4)[1])

clean_data1 = find_good_fits(data1)
clean_data2 = find_good_fits(data2)

event_ID1, decay_time1, hits1, layerID1, columnID1, cellID1, pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b, upos1a, upos1b, uwid1a, uwid1b, uamp1a, uamp1b, ur1a, ur1b, chisq1, bins1, params1, p1 = get_fit_data(clean_data1, hits_range, decay_time_range)
event_ID2, decay_time2, hits2, layerID2, columnID2, cellID2, pos2a, pos2b, wid2a, wid2b, amp2a, amp2b, r2a, r2b, upos2a, upos2b, uwid2a, uwid2b, uamp2a, uamp2b, ur2a, ur2b, chisq2, bins2, params2, p2 = get_fit_data(clean_data2)
event_ID3, decay_time3, hits3, layerID3, columnID3, cellID3, pos3, wid3, amp3, r3, upos3, uwid3, uamp3, ur3, chisq3, bins3, params3, p3 = get_fit_data(data3, hits_range, decay_time_range)
event_ID4, decay_time4, hits4, layerID4, columnID4, cellID4, pos4, wid4, amp4, r4, upos4, uwid4, uamp4, ur4, chisq4, bins4, params4, p4 = get_fit_data(data4)

print("\nNumber of fits to analyze after filtering")
print(len(event_ID1))
print(len(event_ID2))
print(len(event_ID3))
print(len(event_ID4))

print("\nNumber of events in these fits")
print(len(np.unique(event_ID1)))
print(len(np.unique(event_ID2)))
print(len(np.unique(event_ID3)))
print(len(np.unique(event_ID4)))

chisqdof1 = chisq1/(bins1 - params1)
chisqdof2 = chisq2/(bins2 - params2)
chisqdof3 = chisq3/(bins3 - params3)
chisqdof4 = chisq4/(bins4 - params4)

time_diff1 = abs(pos1a - pos1b)
time_diff2 = abs(pos2a - pos2b)

amp_diff1 = abs(amp1a - amp1b)
amp_diff2 = abs(amp2a - amp2b)

amp_ratio1 = []
amp_ratio2 = []

for i in range(len(amp1a)):
    amp_ratio1.append(max(amp1a[i]/amp1b[i], amp1b[i]/amp1a[i]))

for i in range(len(amp2a)):
    amp_ratio2.append(max(amp2a[i]/amp2b[i], amp2b[i]/amp2a[i]))

wid_diff1 = abs(wid1a - wid1b)
wid_diff2 = abs(wid2a - wid2b)

wid_ratio1 = []
wid_ratio2 = []

for i in range(len(wid1a)):
    wid_ratio1.append(max(wid1a[i]/wid1b[i], wid1b[i]/wid1a[i]))

for i in range(len(wid2a)):
    wid_ratio2.append(max(wid2a[i]/wid2b[i], wid2b[i]/wid2a[i]))


plt.figure(figsize=(10,7))
plt.hist(time_diff1, bins=10, range=(0,36), histtype='step', label='tau events')
plt.hist(time_diff2, bins=10, range=(0,36), histtype='step', label='e- events')
plt.xlabel('Time difference [ns]')
plt.ylabel('Frequency')
plt.xlim(0)
plt.title('Time difference Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/time_diff.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(amp_ratio1, bins=5, range=(1,4), histtype='step', label='tau events')
plt.hist(amp_ratio2, bins=5, range=(1,4), histtype='step', label='e- events')
plt.xlabel('Amplitude Ratio')
plt.ylabel('Frequency')
plt.xlim(1)
plt.title('Amplitude Ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/amp_ratio.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(wid_ratio1, bins=5, range=(1,3), histtype='step', label='tau events')
plt.hist(wid_ratio2, bins=5, range=(1,3), histtype='step', label='e- events')
plt.xlabel('FWHM Ratio')
plt.ylabel('Frequency')
plt.xlim(1)
plt.title('FWHM Ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/wid_diff.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(chisq1, bins=18, range=(0,100), histtype='step', label='tau events')
plt.hist(chisq2, bins=18, range=(0,100), histtype='step', label='e- events')
plt.xlabel('Chi^2')
plt.ylabel('Frequency')
plt.xlim(0)
plt.title('Double Bi_gaussian Chi^2 Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/chisq_double_bi_gaussian.png'.format(energy))
#
plt.figure(figsize=(10,7))
plt.hist(chisq3, bins=20, range=(0,300), histtype='step', label='tau events')
plt.hist(chisq4, bins=20, range=(0,300), histtype='step', label='e- events')
plt.xlabel('Chi^2')
plt.ylabel('Frequency')
plt.xlim(0)
plt.title('Single Bi_gaussian Chi^2 Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/chisq_single_bi_gaussian.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(chisqdof1, bins=15, range=(0,30), histtype='step', label='tau events')
plt.hist(chisqdof2, bins=15, range=(0,30), histtype='step', label='e- events')
plt.xlabel('Chi^2/dof')
plt.ylabel('Frequency')
plt.xlim(0)
plt.title('Double Bi_gaussian Chi^2/dof Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/chisqdof_double_bi_gaussian.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(chisqdof3, bins=20, range=(0,500), histtype='step', label='tau events')
plt.hist(chisqdof4, bins=20, range=(0,500), histtype='step', label='e- events')
plt.xlabel('Chi^2/dof')
plt.ylabel('Frequency')
plt.xlim(0)
plt.title('Single Bi_gaussian Chi^2/dof Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/chisqdof_single_bi_gaussian.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(p1, bins=15, range=(0,0.03), histtype='step', label='tau events')
plt.hist(p2, bins=15, range=(0,0.03), histtype='step', label='e- events')
plt.xlabel('p-value')
plt.ylabel('Frequency')
plt.xlim(0)
plt.title('Double Bi_gaussian p-value Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/p-value_double_bi_gaussian.png'.format(energy))

# trying chi^2 single/ chi^2 double, also p double/ p single

event_ID_tau = []
event_ID_e = []
keep_tau = []
keep_e = []

for i in range(len(event_ID1)):
    id = [event_ID1[i], layerID1[i], columnID1[i], cellID1[i]]
    event_ID_tau.append(id)

for i in range(len(event_ID2)):
    id = [event_ID2[i], layerID2[i], columnID2[i], cellID2[i]]
    event_ID_e.append(id)

for i in range(len(event_ID3)):
    id = [event_ID3[i], layerID3[i], columnID3[i], cellID3[i]]
    if id in event_ID_tau:
        keep_tau.append(True)
    else:
        keep_tau.append(False)

for i in range(len(event_ID4)):
    id = [event_ID4[i], layerID4[i], columnID4[i], cellID4[i]]
    if id in event_ID_e:
        keep_e.append(True)
    else:
        keep_e.append(False)

event_ID3, decay_time3, hits3, layerID3, columnID3, cellID3, pos3, wid3, amp3, r3, upos3, uwid3, uamp3, ur3, chisq3, bins3, params3, p3 = event_ID3[keep_tau], decay_time3[keep_tau], hits3[keep_tau], layerID3[keep_tau], columnID3[keep_tau], cellID3[keep_tau], pos3[keep_tau], wid3[keep_tau], amp3[keep_tau], r3[keep_tau], upos3[keep_tau], uwid3[keep_tau], uamp3[keep_tau], ur3[keep_tau], chisq3[keep_tau], bins3[keep_tau], params3[keep_tau], p3[keep_tau]
event_ID4, decay_time4, hits4, layerID4, columnID4, cellID4, pos4, wid4, amp4, r4, upos4, uwid4, uamp4, ur4, chisq4, bins4, params4, p4 = event_ID4[keep_e], decay_time4[keep_e], hits4[keep_e], layerID4[keep_e], columnID4[keep_e], cellID4[keep_e], pos4[keep_e], wid4[keep_e], amp4[keep_e], r4[keep_e], upos4[keep_e], uwid4[keep_e], uamp4[keep_e], ur4[keep_e], chisq4[keep_e], bins4[keep_e], params4[keep_e], p4[keep_e]

chisq1 = chisq1[p3 > 0]
chisq2 = chisq2[p4 > 0]
chisq3 = chisq3[p3 > 0]
chisq4 = chisq4[p4 > 0]
p1 = p1[p3 > 0]
p2 = p2[p4 > 0]
p3 = p3[p3 > 0]
p4 = p4[p4 > 0]

# print(p4)
# print(event_ID2[12])
# print(layerID2[12])
# print(columnID2[12])
# print(cellID2[12])

plt.figure(figsize=(10,7))
plt.hist(chisq1*chisq3, histtype='step', label='tau events')
plt.hist(chisq2*chisq4, histtype='step', label='e- events')
plt.xlabel('Chi^2 Single Bi_gaussian / Chi^2 Double Bi_gaussian')
plt.ylabel('Frequency')
# plt.xlim(0)
plt.title('Chi^2 ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/chi^2_ratio.png'.format(energy))

bins = np.linspace(1, np.log10(max(p1/p3)), 15)
plt.figure(figsize=(10,7))
plt.hist(p1/p3, bins=10**bins, histtype='step', label='tau events')
plt.hist(p2/p4, bins=10**bins, histtype='step', label='e- events')
plt.xlabel('p-value Single Bi_gaussian/ p-value Double Bi_gaussian')
plt.ylabel('Frequency')
# plt.xlim(0)
plt.xscale('log')
plt.title('p-value ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/bi_gaussian/p-value_ratio.png'.format(energy))

n1, bins1 = np.histogram(time_diff1, bins=20, range=(0,100))
n2, bins2 = np.histogram(time_diff2, bins=20, range=(0,100))

fpr1 = []
tpr1 = []
tau = e = 0
N_tau = len(time_diff1)
N_e = len(time_diff2)
N_tot = N_e + N_tau
for i in range(len(bins1)-1):
    tau += n1[-i - 1]
    e += n2[-i - 1]
    tpr1.append(tau/N_tau)
    fpr1.append(e/N_e)

fpr2 = []
tpr2 = []
bins = np.linspace(1, np.log10(max(p1/p3)*1.2), 20)
n1, bins1 = np.histogram(p1/p3, bins=10**bins)
n2, bins2 = np.histogram(p2/p4, bins=10**bins)

tau = e = 0
N_tau = len(p1)
N_e = len(p2)
N_tot = N_e + N_tau
for i in range(len(bins1)-1):
    tau += n1[-i - 1]
    e += n2[-i - 1]
    tpr2.append(tau/N_tau)
    fpr2.append(e/N_e)

fpr3 = []
tpr3 = []
bins = np.linspace(1, np.log10(max(p1/p3)), 20)
n1, bins1 = np.hist(p1/p3, bins=10**bins)
n1, bins1 = np.hist(p2/p4, bins=10**bins)

tau = e = 0
N_tau = len(time_diff1)
N_e = len(time_diff2)
N_tot = N_e + N_tau
for i in range(len(bins1)-1):
    tau += n1[-i - 1]
    e += n2[-i - 1]
    tpr3.append(tau/N_tau)
    fpr3.append(e/N_e)

fpr2 = []
tpr2 = []
bins = np.linspace(1, np.log10(max(p1/p3)), 20)
n1, bins1 = np.hist(p1/p3, bins=10**bins)
n1, bins1 = np.hist(p2/p4, bins=10**bins)

tau = e = 0
N_tau = len(time_diff1)
N_e = len(time_diff2)
N_tot = N_e + N_tau
for i in range(len(bins1)-1):
    tau += n1[-i - 1]
    e += n2[-i - 1]
    tpr2.append(tau/N_tau)
    fpr2.append(e/N_e)


plt.plot(fpr1, tpr1, label='ROC time_difference')
plt.plot(fpr2, tpr2, label='ROC p-value ratio')
plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend()
plt.show()
