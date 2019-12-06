"""
This script fits the time residuals of the tau and e- events for a single
sensor. The type of fit and the event needs to be selected by the user
by commenting or uncommenting chunks of the code. The script takes two
parameters the energy and event ID. Mostly just used if for a
larger/nicer plot of a specific fit.
"""

import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare

energy = int(sys.argv[1])
event_ID = int(sys.argv[2])

fname1 = '{} TeV Tau/tau_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, event_ID)
fname2 = '{} TeV Tau/{}_TeV_tau_decay_times.csv'.format(energy, energy)

data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')

decay_times = data2[:,2]
decay_time = decay_times[event_ID]

# specify sensor here
a = 1
b = 1
PMT_ID = [a, b, 1]

fname3 = '{} TeV Final Results/base_e_had_event/string{}{}/{}TeV_e_had_merge_bins_PMT_ID{}{}{}.csv'.format(energy, a, b, energy, *PMT_ID)

data3 = np.loadtxt(fname3, delimiter=',', comments='#')

e_n, had_n, bins = data3[:,0], data3[:,1], data3[:,2]
layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(data1, PMT_ID)

bin_size = bins[1] - bins[0]
last_bin = bins[-1] + bin_size
bins = np.append(bins, last_bin)

PMT_pos = PMT_ID_to_pos(PMT_ID)
dx, dy, dz = seperation_vector(PMT_pos)
d = distance_to_vertex(dx, dy, dz)
theta = angle_to_vertex(dx, dy, dz)

np.random.seed(event_ID)
# e- and had are smooth averages of 100 events, alter by drawing from a poisson
# distribution with lambda equal to n
e_n = np.random.poisson(e_n)
had_n = np.random.poisson(had_n)

n1, bins1 = np.histogram(time1, bins=bins)
n2 = had_n
n3 = e_n
n4 = n1 + had_n
n5 = e_n + had_n

peak_height1, peak_pos1 = get_peak(n1, bins)
peak_height2, peak_pos2 = get_peak(n2, bins)
peak_height3, peak_pos3 = get_peak(n3, bins)
peak_height4, peak_pos4 = get_peak(n4, bins)
peak_height5, peak_pos5 = get_peak(n5, bins)

FWHM1 = get_FWHM(n1, bins)
FWHM2 = get_FWHM(n2, bins)
FWHM3 = get_FWHM(n3, bins)
FWHM4 = get_FWHM(n4, bins)
FWHM5 = get_FWHM(n5, bins)

tmin4, tmax4 = get_range_fit(n4, bins, peak_height4)
tmin5, tmax5 = get_range_fit(n5, bins, peak_height5)

bins = bins[:-1]

range4 = (bins >= tmin4) & (bins <= tmax4)
range5 = (bins >= tmin5) & (bins <= tmax5)

bins4 = bins[range4]
bins5 = bins[range5]

n4 = n4[range4]
n5 = n5[range5]

if min(bins4) < 6:

sigma4 = np.sqrt(n4)
sigma5 = np.sqrt(n5)

guesses1 = [peak_pos1, peak_pos2, FWHM1, FWHM2, peak_height1, peak_height2]
guesses2 = [peak_pos4, FWHM4, peak_height4]
guesses3 = [peak_pos3, peak_pos2, FWHM3, FWHM2, peak_height3, peak_height2]
guesses4 = [peak_pos5, FWHM5, peak_height5]

# guesses1 = [peak_pos1, peak_pos2, FWHM1, FWHM2, peak_height1, peak_height2, 0.3, 0.3]
# guesses2 = [peak_pos4, FWHM4, peak_height4, 0.3]
# guesses3 = [peak_pos3, peak_pos2, FWHM3, FWHM2, peak_height3, peak_height2, 0.3, 0.3]
# guesses4 = [peak_pos5, FWHM5, peak_height5, 0.3]

fitparams1, fitcov1 = curve_fit(fit_function_gaussian, bins4, n4, sigma=sigma4, bounds=(0., np.inf), p0=guesses1, maxfev=100000)
fitparams2, fitcov2 = curve_fit(gaussian, bins4, n4, sigma=sigma4, bounds=(0., np.inf), p0=guesses2, maxfev=100000)
fitparams3, fitcov3 = curve_fit(fit_function_gaussian, bins5, n5, sigma=sigma5, bounds=(0., np.inf), p0=guesses3, maxfev=100000)
fitparams4, fitcov4 = curve_fit(gaussian, bins5, n5, sigma=sigma5, bounds=(0., np.inf), p0=guesses4, maxfev=100000)

# fitparams1, fitcov1 = curve_fit(fit_function, bins4, n4, sigma=sigma4, bounds=(0., np.inf), p0=guesses1, maxfev=100000)
# fitparams2, fitcov2 = curve_fit(bi_gaussian, bins4, n4, sigma=sigma4, bounds=(0., np.inf), p0=guesses2, maxfev=100000)
# fitparams3, fitcov3 = curve_fit(fit_function, bins5, n5, sigma=sigma5, bounds=(0., np.inf) ,p0=guesses3, maxfev=100000)
# fitparams4, fitcov4 = curve_fit(bi_gaussian, bins5, n5, sigma=sigma5, bounds=(0., np.inf), p0=guesses4, maxfev=100000)

fitparams1_error = np.sqrt(np.diag(fitcov1))
fitparams2_error = np.sqrt(np.diag(fitcov2))
fitparams3_error = np.sqrt(np.diag(fitcov3))
fitparams4_error = np.sqrt(np.diag(fitcov4))

dof1 = len(bins4) - len(guesses1)
dof2 = len(bins4) - len(guesses2)
dof3 = len(bins5) - len(guesses3)
dof4 = len(bins5) - len(guesses4)

chisq1, p1 = chisquare(n4, fit_function_gaussian(bins4, *fitparams1), ddof=dof1)
chisq2, p2 = chisquare(n4, gaussian(bins4, *fitparams2), ddof=dof2)
chisq3, p3 = chisquare(n5, fit_function_gaussian(bins5, *fitparams3), ddof=dof3)
chisq4, p4 = chisquare(n5, gaussian(bins5, *fitparams4), ddof=dof4)

# chisq1, p1 = chisquare(n4, fit_function(bins4, *fitparams1), ddof=dof1)
# chisq2, p2 = chisquare(n4, bi_gaussian(bins4, *fitparams2), ddof=dof2)
# chisq3, p3 = chisquare(n5, fit_function(bins5, *fitparams3), ddof=dof3)
# chisq4, p4 = chisquare(n5, bi_gaussian(bins5, *fitparams4), ddof=dof4)

chisqdof1 = chisq1/dof1
chisqdof2 = chisq2/dof2
chisqdof3 = chisq3/dof3
chisqdof4 = chisq4/dof4

pos1a, pos1b, wid1a, wid1b, amp1a, amp1b = fitparams1
pos2a, pos2b, wid2a, wid2b, amp2a, amp2b = fitparams3

# pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b = fitparams1
# pos2a, pos2b, wid2a, wid2b, amp2a, amp2b, r2a, r2b = fitparams3

# print("1")
# print(fitparams1)
# print(fitparams1_error)
# print("\n2")
# print(fitparams2)
# print(fitparams2_error)
# print("\n3")
# print(fitparams3)
# print(fitparams3_error)
# print("\n4")
# print(fitparams4)
# print(fitparams4_error)
# print(chisq1, chisq2, chisq3, chisq4)
# print(p1, p2, p3, p4)


plt.figure(figsize=(12,8))
plt.suptitle('{} TeV Time Residuals with Gaussian Fits; tau_decay_time = {:.2f}\nPMT_ID = [{}, {}, {}], distance = {:.4f} m, angle={:.4f} degrees'.format(energy, decay_time, *PMT_ID, d, theta), fontsize=14)

plt.subplot(221)
plt.plot(bins4, gaussian(bins4, pos1a, wid1a, amp1a), label='Gaussian 1')
plt.plot(bins4, gaussian(bins4, pos1b, wid1b, amp1b), label='Gaussian 2')
plt.step(bins4, n4, where='mid', color='C2', label='Sum')
plt.plot(bins4, fit_function_gaussian(bins4, *fitparams1), color='C3', label='Double Gaussian Fit')
plt.xlim(tmin4, tmax4)
plt.ylim(0)
plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}, p = {:.2f}'.format(chisq1, chisqdof1, len(bins4), len(guesses1), p1), fontsize=11)
plt.legend(fontsize=10)

plt.subplot(222)
plt.plot(bins5, gaussian(bins5, pos2a, wid2a, amp2a), label='Gaussian 1')
plt.plot(bins5, gaussian(bins5, pos2b, wid2b, amp2b), label='Gaussian 2')
plt.step(bins5, n5, where='mid', color='C2', label='Sum')
plt.plot(bins5, fit_function_gaussian(bins5, *fitparams3), color='C3', label='Double Gaussian Fit')
plt.xlim(tmin5, tmax5)
plt.ylim(0)
plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}, p = {:.2f}'.format(chisq3, chisqdof3, len(bins5), len(guesses3), p3), fontsize=11)
plt.legend(fontsize=10)

plt.subplot(223)
plt.step(bins4, n4, where='mid', color='C2', label='Sum')
plt.plot(bins4, gaussian(bins4, *fitparams2), color='C3', label='Single Gaussian Fit')
plt.xlim(tmin4, tmax4)
plt.ylim(0)
plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}, p = {:.2f}'.format(chisq2, chisqdof2, len(bins4), len(guesses2), p2), fontsize=11)
plt.legend(fontsize=10)

plt.subplot(224)
plt.step(bins5, n5, where='mid', color='C2', label='Sum')
plt.plot(bins5, gaussian(bins5, *fitparams4), color='C3', label='Single Gaussian Fit')
plt.xlim(tmin5, tmax5)
plt.ylim(0)
plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}, p = {:.2f}'.format(chisq4, chisqdof4, len(bins5), len(guesses4), p4), fontsize=11)
plt.legend(fontsize=10)
plt.savefig('fit_test{}.png'.format(event_ID))

##########################################################################
# plt.suptitle('{} TeV Time Residuals with Bi_gaussian Fits; tau_decay_time = {:.2f}\nPMT_ID = [{}, {}, {}], distance = {:.4f} m, angle={:.4f} degrees'.format(energy, decay_time, *PMT_ID, d, theta), fontsize=14)
#
# plt.subplot(221)
# plt.plot(bins4, bi_gaussian(bins4, pos1a, wid1a, amp1a, r1a), label='Gaussian 1')
# plt.plot(bins4, bi_gaussian(bins4, pos1b, wid1b, amp1b, r1b), label='Gaussian 2')
# plt.step(bins4, n4, where='mid', color='C2', label='Sum')
# plt.plot(bins4, fit_function(bins4, *fitparams1), color='C3', label='Double Bi_gaussian Fit')
# plt.xlim(tmin4, tmax4)
# plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}, p = {:.2f}'.format(chisq1, chisqdof1, len(bins4), len(guesses1), p1), fontsize=11)
# plt.legend(fontsize=10)
#
# plt.subplot(222)
# plt.plot(bins5, bi_gaussian(bins5, pos2a, wid2a, amp2a, r2a), label='Gaussian 1')
# plt.plot(bins5, bi_gaussian(bins5, pos2b, wid2b, amp2b, r2b), label='Gaussian 2')
# plt.step(bins5, n5, where='mid', color='C2', label='Sum')
# plt.plot(bins5, fit_function(bins5, *fitparams3), color='C3', label='Double Bi_gaussian Fit')
# plt.xlim(tmin5, tmax5)
# plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}, p = {:.2f}'.format(chisq3, chisqdof3, len(bins5), len(guesses3), p3), fontsize=11)
# plt.legend(fontsize=10)
#
# plt.subplot(223)
# plt.step(bins4, n4, where='mid', color='C2', label='Sum')
# plt.plot(bins4, bi_gaussian(bins4, *fitparams2), color='C3', label='Single Bi_gaussian Fit')
# plt.xlim(tmin4, tmax4)
# plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}, p = {:.2f}'.format(chisq2, chisqdof2, len(bins4), len(guesses2), p2), fontsize=11)
# plt.legend(fontsize=10)
#
# plt.subplot(224)
# plt.step(bins5, n5, where='mid', color='C2', label='Sum')
# plt.plot(bins5, bi_gaussian(bins5, *fitparams4), color='C3', label='Single Bi_gaussian Fit')
# plt.xlim(tmin5, tmax5)
# plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}, p = {:.2f}'.format(chisq4, chisqdof4, len(bins5), len(guesses4), p4), fontsize=11)
# plt.legend(fontsize=10)
#
# plt.subplots_adjust(top=0.85, bottom=0.05, left=0.08, right=0.93, hspace=0.35, wspace=0.15)
# plt.savefig('fit_test{}_bi.png'.format(event_ID))
