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
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/50 TeV/fitting_params_double_bi-gaussian_e.csv'#.format(energy)
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params_single_bi-gaussian_tau.csv'.format(energy)
fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/50 TeV/fitting_params_single_bi-gaussian_e.csv'#.format(energy)

data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')
data3 = np.loadtxt(fname3, delimiter=',', comments='#')
data4 = np.loadtxt(fname4, delimiter=',', comments='#')

hits_range = (40000, 360000)
decay_time_range = (0, 60)

print("Number of inf in chi^2")
print(len(data1[:,22][np.isinf(data1[:,22])]))
print(len(data2[:,22][np.isinf(data2[:,22])]))
print(len(data3[:,14][np.isinf(data3[:,14])]))
print(len(data4[:,14][np.isinf(data4[:,14])]))

print("\nTotal number of fits for each function")
print(np.shape(data1)[0])
print(np.shape(data2)[0])
print(np.shape(data3)[0])
print(np.shape(data4)[0])

clean_data1, clean_data3 = find_good_fits(data1, data3)
clean_data2, clean_data4 = find_good_fits(data2, data4)

print("\nNumber of fits after filtering")
print(len(clean_data1[0][:]))
print(len(clean_data2[0][:]))
print(len(clean_data3[0][:]))
print(len(clean_data4[0][:]))

event_ID1, decay_time1, hits1, layerID1, columnID1, cellID1, pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b, upos1a, upos1b, uwid1a, uwid1b, uamp1a, uamp1b, ur1a, ur1b, chisq1, bins1, params1, p1 = get_fit_data(clean_data1, hits_range, decay_time_range)
event_ID2, decay_time2, hits2, layerID2, columnID2, cellID2, pos2a, pos2b, wid2a, wid2b, amp2a, amp2b, r2a, r2b, upos2a, upos2b, uwid2a, uwid2b, uamp2a, uamp2b, ur2a, ur2b, chisq2, bins2, params2, p2 = get_fit_data(clean_data2)
event_ID3, decay_time3, hits3, layerID3, columnID3, cellID3, pos3, wid3, amp3, r3, upos3, uwid3, uamp3, ur3, chisq3, bins3, params3, p3 = get_fit_data(clean_data3, hits_range, decay_time_range)
event_ID4, decay_time4, hits4, layerID4, columnID4, cellID4, pos4, wid4, amp4, r4, upos4, uwid4, uamp4, ur4, chisq4, bins4, params4, p4 = get_fit_data(clean_data4)

# event_ID1, decay_time1, hits1, layerID1, columnID1, cellID1, pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b, upos1a, upos1b, uwid1a, uwid1b, uamp1a, uamp1b, ur1a, ur1b, chisq1, bins1, params1, p1 = event_ID1[:-1], decay_time1[:-1], hits1[:-1], layerID1[:-1], columnID1[:-1], cellID1[:-1], pos1a[:-1], pos1b[:-1], wid1a[:-1], wid1b[:-1], amp1a[:-1], amp1b[:-1], r1a[:-1], r1b[:-1], upos1a[:-1], upos1b[:-1], uwid1a[:-1], uwid1b[:-1], uamp1a[:-1], uamp1b[:-1], ur1a[:-1], ur1b[:-1], chisq1[:-1], bins1[:-1], params1[:-1], p1[:-1]

# event_ID3, decay_time3, hits3, layerID3, columnID3, cellID3, pos3, wid3, amp3, r3, upos3, uwid3, uamp3, ur3, chisq3, bins3, params3, p3 = event_ID3[:-3], decay_time3[:-3], hits3[:-3], layerID3[:-3], columnID3[:-3], cellID3[:-3], pos3[:-3], wid3[:-3], amp3[:-3], r3[:-3], upos3[:-3], uwid3[:-3], uamp3[:-3], ur3[:-3], chisq3[:-3], bins3[:-3], params3[:-3], p3[:-3]

print("\nNumber of fits to analyze in hit range and decay time range")
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
plt.hist(time_diff1, bins=20, range=(0,100), density=True, histtype='step', label='tau events')
plt.hist(time_diff2, bins=20, range=(0,100), density=True, histtype='step', label='e- events')
plt.xlabel('Time difference [ns]')
plt.ylabel('Normalized bincount')
plt.xlim(0)
plt.title('Time difference Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/time_diff.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(amp_ratio1, bins=5, range=(0.25,4), density=True, histtype='step', label='tau events')
plt.hist(amp_ratio2, bins=5, range=(0.25,4), density=True, histtype='step', label='e- events')
plt.xlabel('Amplitude Ratio')
plt.ylabel('Normalized bincount')
plt.xlim(0.25,4)
plt.title('Amplitude Ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/amp_ratio.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(wid_ratio1, bins=5, range=(0.25,4), density=True, histtype='step', label='tau events')
plt.hist(wid_ratio2, bins=5, range=(0.25,4), density=True, histtype='step', label='e- events')
plt.xlabel('FWHM Ratio')
plt.ylabel('Normalized bincount')
plt.xlim(0.25,4)
plt.title('FWHM Ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/wid_ratio.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(chisq1, bins=20, range=(0,200), density=True, histtype='step', label='tau events')
plt.hist(chisq2, bins=20, range=(0,200), density=True, histtype='step', label='e- events')
plt.xlabel('Chi^2')
plt.ylabel('Normalized bincount')
plt.xlim(0)
plt.title('Double Bi_gaussian Chi^2 Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chisq_double_bi_gaussian.png'.format(energy))
#
plt.figure(figsize=(10,7))
plt.hist(chisq3, bins=20, range=(0,1000), density=True, histtype='step', label='tau events')
plt.hist(chisq4, bins=20, range=(0,1000), density=True, histtype='step', label='e- events')
plt.xlabel('Chi^2')
plt.ylabel('Normalized bincount')
plt.xlim(0)
plt.title('Single Bi_gaussian Chi^2 Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chisq_single_bi_gaussian.png'.format(energy))

# plt.figure(figsize=(10,7))
# plt.hist(chisqdof1, bins=60, range=(0,30), density=True, histtype='step', label='tau events')
# plt.hist(chisqdof2, bins=60, range=(0,30), density=True, histtype='step', label='e- events')
# plt.xlabel('Chi^2/dof')
# plt.ylabel('Normalized bincount')
# plt.xlim(0)
# plt.title('Double Bi_gaussian Chi^2/dof Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
# plt.legend()
# plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chisqdof_double_bi_gaussian.png'.format(energy))
#
# plt.figure(figsize=(10,7))
# plt.hist(chisqdof3, bins=50, range=(0,100), density=True, histtype='step', label='tau events')
# plt.hist(chisqdof4, bins=50, range=(0,100), density=True, histtype='step', label='e- events')
# plt.xlabel('Chi^2/dof')
# plt.ylabel('Normalized bincount')
# plt.xlim(0)
# plt.title('Single Bi_gaussian Chi^2/dof Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
# plt.legend()
# plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chisqdof_single_bi_gaussian.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(p1, bins=20, range=(0,1), density=True, histtype='step', label='tau events')
plt.hist(p2, bins=20, range=(0,1), density=True, histtype='step', label='e- events')
plt.xlabel('p-value')
plt.ylabel('Normalized bincount')
plt.xlim(0)
plt.title('Double Bi_gaussian p-value Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/p-value_double_bi_gaussian.png'.format(energy))

plt.figure(figsize=(10,7))
plt.hist(p3, bins=50, range=(0,1), density=True, histtype='step', label='tau events')
plt.hist(p4, bins=50, range=(0,1), density=True, histtype='step', label='e- events')
plt.xlabel('p-value')
plt.ylabel('Normalized bincount')
plt.xlim(0)
plt.title('Single Bi_gaussian p-value Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/p-value_single_bi_gaussian.png'.format(energy))

# removing where p3, p4 = 0 so no divide by zero error
chisq1 = chisq1[(p3 > 0) & (p1 > 0)]
chisq2 = chisq2[(p4 > 0) & (p2 > 0)]
chisq3 = chisq3[(p3 > 0) & (p1 > 0)]
chisq4 = chisq4[(p4 > 0) & (p2 > 0)]

time_diff1 = time_diff1[(p3 > 0) & (p1 > 0)]
time_diff2 = time_diff2[(p4 > 0) & (p2 > 0)]

params1 = params1[(p3 > 0) & (p1 > 0)]
params2 = params2[(p4 > 0) & (p2 > 0)]
params3 = params3[(p3 > 0) & (p1 > 0)]
params4 = params4[(p4 > 0) & (p2 > 0)]

bins1 = bins1[(p3 > 0) & (p1 > 0)]
bins2 = bins2[(p4 > 0) & (p2 > 0)]
bins3 = bins3[(p3 > 0) & (p1 > 0)]
bins4 = bins4[(p4 > 0) & (p2 > 0)]

select1 = (p3 > 0) & (p1 > 0)
select2 = (p4 > 0) & (p2 > 0)

p1 = p1[select1]
p2 = p2[select2]
p3 = p3[select1]
p4 = p4[select2]

chisqdof1 = chisq1/(bins1 - params1)
chisqdof2 = chisq2/(bins2 - params2)
chisqdof3 = chisq3/(bins3 - params3)
chisqdof4 = chisq4/(bins4 - params4)

print("\nNumber of fits to analyze after removing ones with p=0")
print(len(p1))
print(len(p2))
print(len(p1))
print(len(p2))

plt.figure(figsize=(10,7))
plt.hist(chisq1/chisq3, bins=40, range=(0,2), density=True, histtype='step', label='tau events')
plt.hist(chisq2/chisq4, bins=40, range=(0,2), density=True, histtype='step', label='e- events')
plt.xlabel('Chi^2 Double Bi_gaussian / Chi^2 Single Bi_gaussian')
plt.ylabel('Normalized bincount')
plt.xlim(0)
plt.title('Chi^2 ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chi^2_ratio.png'.format(energy))

# plt.figure(figsize=(10,7))
# plt.hist(chisqdof1/chisqdof3, bins=50, range=(0,5), density=True, histtype='step', label='tau events')
# plt.hist(chisqdof2/chisqdof4, bins=50, range=(0,5), density=True, histtype='step', label='e- events')
# plt.xlabel('Chi^2/dof Double Bi_gaussian / Chi^2/dof Single Bi_gaussian')
# plt.ylabel('Normalized bincount')
# plt.xlim(0)
# plt.title('Chi^2/dof ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
# plt.legend()
# plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chi^2_dof_ratio.png'.format(energy))


bins = np.linspace(0.5, np.log10(1.1*max(p1/p3)), 20)
plt.figure(figsize=(10,7))
plt.hist(np.log10(p1/p3), bins=bins, density=True, histtype='step', label='tau events')
plt.hist(np.log10(p2/p4), bins=bins, density=True, histtype='step', label='e- events')
plt.xlabel('p-value Double Bi_gaussian/ p-value Single Bi_gaussian log scale')
plt.ylabel('Normalized bincount')
plt.xlim(0)
plt.title('p-value ratio Distribution; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/p-value_ratio.png'.format(energy))


#### Scatterplots ###
plt.figure(figsize=(10,7))
plt.scatter(p1, p3, label='tau events')
plt.scatter(p2, p4, label='e- events')
plt.xlabel('p-value double')
plt.ylabel('p-value single')
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend()
plt.title('p-value Scatterplot')
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/p-value_scatterplot.png'.format(energy))

plt.figure(figsize=(10,7))
plt.scatter(time_diff1, np.log10(p1/p3), label='tau events')
plt.scatter(time_diff2, np.log10(p2/p4), label='e- events')
plt.xlabel('time difference [ns]')
plt.ylabel('p-value ratio log')
plt.xlim(0,100)
plt.ylim(0,np.log10(1.1*max(p1/p3)))
plt.legend()
plt.title('p-value ratio vs time_diff Scatterplot')
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/p-value_time_diff_scaterplot.png'.format(energy))

plt.figure(figsize=(10,7))
plt.scatter(chisq1, chisq3, label='tau events')
plt.scatter(chisq2, chisq4, label='e- events')
plt.xlabel('chisq double')
plt.ylabel('chisq single')
plt.xlim(0,200)
plt.ylim(0,1000)
plt.legend()
plt.title('chisq Scatterplot')
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chisq_scaterplot.png'.format(energy))

plt.figure(figsize=(10,7))
plt.scatter(chisq1, time_diff1, label='tau events')
plt.scatter(chisq2, time_diff2, label='e- events')
plt.xlabel('chisq double')
plt.ylabel('time_diff')
plt.xlim(0,200)
plt.ylim(0,100)
plt.legend()
plt.title('chisq double vs time diff Scatterplot')
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chisq_time_diff_scaterplot.png'.format(energy))

plt.figure(figsize=(10,7))
plt.scatter(chisq1/chisq3, time_diff1, label='tau events')
plt.scatter(chisq2/chisq4, time_diff2, label='e- events')
plt.xlabel('chisq ratio')
plt.ylabel('time_diff')
plt.xlim(0,2)
plt.ylim(0,100)
plt.legend()
plt.title('chisq ratio vs time diff Scatterplot')
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chisq_ratio_time_diff_scaterplot.png'.format(energy))

plt.figure(figsize=(10,7))
plt.scatter(chisq1/chisq3, np.log10(p1/p3), label='tau events')
plt.scatter(chisq2/chisq4, np.log10(p2/p4), label='e- events')
plt.xlabel('chisq ratio')
plt.ylabel('p-value ratio log')
plt.xlim(0,2)
plt.ylim(0, np.log10(1.1*max(p1/p3)))
plt.legend()
plt.title('chisq ratio vs p-value ratio Scatterplot')
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/chisq_p_ratio_scaterplot.png'.format(energy))

# plt.show()


N_tau = len(p1)
N_e = len(p2)
N_tot = N_e + N_tau

bins = np.linspace(0.5, np.log10(1.1*max(p1/p3)), 40)
n1, bins1 = np.histogram(time_diff1, bins=25, range=(0,100))
n2, bins2 = np.histogram(time_diff2, bins=25, range=(0,100))
n3, bins3 = np.histogram(chisq1/chisq3, bins=40, range=(0,2))
n4, bins4 = np.histogram(chisq2/chisq4, bins=40, range=(0,2))
n5, bins5 = np.histogram(np.log10(p1/p3), bins=bins)
n6, bins6 = np.histogram(np.log10(p2/p4), bins=bins)

fpr1 = []
tpr1 = []
tau = e = 0
for i in range(len(bins1)-1):
    tau += n1[-i - 1]
    e += n2[-i - 1]
    tpr1.append(tau/N_tau)
    fpr1.append(e/N_e)

fpr2 = [0]
tpr2 = [0]
tau = e = 0
for i in range(len(bins3)-1):
    tau += n3[i]
    e += n4[i]
    tpr2.append(tau/N_tau)
    fpr2.append(e/N_e)

fpr3 = []
tpr3 = []
tau = e = 0
for i in range(len(bins5)-1):
    tau += n5[-i - 1]
    e += n6[-i - 1]
    tpr3.append(tau/N_tau)
    fpr3.append(e/N_e)

tpr1.append(1)
fpr1.append(1)
tpr2.append(1)
fpr2.append(1)
tpr3.append(1)
fpr3.append(1)

auc1 = np.trapz(tpr1, fpr1)
auc2 = np.trapz(tpr2, fpr2)
auc3 = np.trapz(tpr3, fpr3)

plt.figure(figsize=(10,7))
plt.plot(fpr1, tpr1, label='ROC time_difference AUC={:.2f}'.format(auc1))
plt.plot(fpr2, tpr2, label='ROC chi^2 ratio AUC={:.2f}'.format(auc2))
plt.plot(fpr3, tpr3, label='ROC p-value ratio AUC={:.2f}'.format(auc3))
plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve; decay_range:{}-{} ns, hit_range: {}-{}'.format(*decay_time_range, *hits_range))
plt.legend()
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/roc_curve.png'.format(energy))
# plt.show()
