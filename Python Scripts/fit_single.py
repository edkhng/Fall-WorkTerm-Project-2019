import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare

neutrino_energy = 100

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_{}TeV377_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/had_{}TeV411_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
# fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_{}TeV371_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)

fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)
# fname5 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)

PMT_ID = [1, 1, 3]
layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(fname1, PMT_ID)
layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(fname2, PMT_ID)
# layerID3, columnID3, cellID3, time3, x3, y3, z3, energy3 = get_data(fname3, PMT_ID)
layerID4, columnID4, cellID4, time4, x4, y4, z4, energy4 = get_data(fname4, PMT_ID)
# layerID5, columnID5, cellID5, time5, x5, y5, z5, energy5 = get_data(fname5, PMT_ID)

PMT_pos = PMT_ID_to_pos(PMT_ID)
dx, dy, dz = seperation_vector('A', PMT_pos)
d = distance_to_vertex(dx, dy, dz)
theta = angle_to_vertex('A', dx, dy, dz)

time1 = clean_data(time1, d)
time2 = clean_data(time2, d)
# time3 = clean_data(time3, d)
time4 = clean_data(time4, d)
# time5 = clean_data(time5, d)

bin_size = 1

bins = []
time = [time1, time2, time4]
time = [time2, time3, time5]
for j in range(3):
    bin = int((max(time[j]) - min(time[j])) / bin_size)
    bins.append(bin)


n1, bins1 = np.histogram(time1, bins=bins[0], range=(min(time1), max(time1)))
n2, bins2 = np.histogram(time2, bins=bins[1], range=(min(time2), max(time2)))
# n3, bins3 = np.histogram(time3, bins=bins[2], range=(min(time3), max(time3)))
n4, bins4 = np.histogram(time4, bins=bins[3], range=(min(time4), max(time4)))
# n5, bins5 = np.histogram(time5, bins=bins[4], range=(min(time5), max(time5)))

peak_height1, peak_pos1 = get_peak(n1, bins1)
peak_height2, peak_pos2 = get_peak(n2, bins2)
# peak_height3, peak_pos3 = get_peak(n3, bins3)
peak_height4, peak_pos4 = get_peak(n4, bins4)
# peak_height5, peak_pos5 = get_peak(n5, bins5)

FWHM1 = get_FWHM(n1, bins1)
FWHM2 = get_FWHM(n2, bins2)
# FWHM3 = get_FWHM(n3, bins3)
FWHM4 = get_FWHM(n4, bins4)
# FWHM5 = get_FWHM(n5, bins5)

tmin1, tmax1 = select_time_range(n1, bins1, peak_height1)
tmin2, tmax2 = select_time_range(n2, bins2, peak_height2)
# tmin3, tmax3 = select_time_range(n3, bins3, peak_height3)
tmin4, tmax4 = select_time_range(n4, bins4, peak_height4)
# tmin5, tmax5 = select_time_range(n5, bins5, peak_height5)

N1 = range_size(time1, tmin1, tmax1)
N2 = range_size(time2, tmin2, tmax2)
# N3 = range_size(time3, tmin3, tmax3)
N4 = range_size(time4, tmin4, tmax4)
# N5 = range_size(time5, tmin5, tmax5)

bin1 = int((tmax1-tmin1)/bin_size)
bin2 = int((tmax2-tmin2)/bin_size)
# bin3 = int((tmax3-tmin3)/bin_size)
bin4 = int((tmax4-tmin4)/bin_size)
# bin5 = int((tmax5-tmin5)/bin_size)

# print("\nThe number of bins: [{}, {}, {}, {}, {}]".format(bin1, bin2, bin3, bin4, bin5))


n1, bins1 = np.histogram(time1, bins=bin1, range=(tmin1, tmax1))
n2, bins2 = np.histogram(time2, bins=bin2, range=(tmin2, tmax2))
# n3, bins3 = np.histogram(time3, bins=bin3, range=(tmin3, tmax3))
n4, bins4 = np.histogram(time4, bins=bin4, range=(tmin4, tmax4))
# n5, bins5 = np.histogram(time5, bins=bin5, range=(tmin5, tmax5))

guesses1 = [peak_pos1, peak_pos2, FWHM1, FWHM2, peak_height1, peak_height2]
guesses2 = [peak_pos4, FWHM4, peak_height4]
guesses3 = [peak_pos1, peak_pos2, FWHM1, FWHM2, peak_height1, peak_height2, 0.3, 0.3]
guesses4 = [peak_pos4, FWHM4, peak_height4, 0.3]

# guesses1 = [peak_pos3, peak_pos2, FWHM3, FWHM2, peak_height3, peak_height2]
# guesses2 = [peak_pos5, FWHM5, peak_height5]
# guesses3 = [peak_pos3, peak_pos2, FWHM3, FWHM2, peak_height3, peak_height2, 0.3, 0.3]
# guesses4 = [peak_pos5, FWHM5, peak_height5, 0.3]

fitparams1, fitcov1 = curve_fit(fit_function_gaussian, bins4[:-1], n4, p0=guesses1, maxfev=100000)
fitparams2, fitcov2 = curve_fit(gaussian, bins4[:-1], n4, p0=guesses2, maxfev=100000)
fitparams3, fitcov3 = curve_fit(fit_function, bins4[:-1], n4, p0=guesses3, maxfev=100000)
fitparams4, fitcov4 = curve_fit(bi_gaussian, bins4[:-1], n4, p0=guesses4, maxfev=100000)

# fitparams1, fitcov1 = curve_fit(fit_function_gaussian, bins5[:-1], n5, p0=guesses1, maxfev=100000)
# fitparams2, fitcov2 = curve_fit(gaussian, bins5[:-1], n5, p0=guesses2, maxfev=100000)
# fitparams3, fitcov3 = curve_fit(fit_function, bins5[:-1], n5, p0=guesses3, maxfev=100000)
# fitparams4, fitcov4 = curve_fit(bi_gaussian, bins5[:-1], n5, p0=guesses4, maxfev=100000)

dof1 = bin4 - len(guesses1)
dof2 = bin5 - len(guesses2)
dof3 = bin4 - len(guesses3)
dof4 = bin5 - len(guesses4)

chisq1, p1 = chisquare(n4, fit_function_gaussian(bins4[:-1], *fitparams1), ddof=dof1)
chisq3, p3 = chisquare(n4, gaussian(bins4[:-1], *fitparams2), ddof=dof2)
chisq2, p2 = chisquare(n4, fit_function(bins4[:-1], *fitparams3), ddof=dof3)
chisq4, p4 = chisquare(n4, bi_gaussian(bins4[:-1], *fitparams4), ddof=dof4)

# chisq1, p1 = chisquare(n5, fit_function_gaussian(bins5[:-1], *fitparams1), ddof=dof1)
# chisq3, p3 = chisquare(n5, gaussian(bins5[:-1], *fitparams2), ddof=dof2)
# chisq2, p2 = chisquare(n5, fit_function(bins5[:-1], *fitparams3), ddof=dof3)
# chisq4, p4 = chisquare(n5, bi_gaussian(bins5[:-1], *fitparams4), ddof=dof4)

chisqdof1 = chisq1/dof3
chisqdof2 = chisq2/dof4
chisqdof3 = chisq3/dof3
chisqdof4 = chisq4/dof4

# print(fitparams)
# print(p)


plt.figure(figsize=(10,7))
# plt.suptitle('{} TeV Time Residuals with Bi_gaussian Fits;\nPMT_ID = [{}, {}, {}],distance = {:.4f} m, angle={:.4f} degrees'.format(neutrino_energy, *PMT_ID, d, theta), fontsize=14)

plt.hist(time1, bins=bins[0], range=(min(time1), max(time1)), histtype='step', color='C0', label='tau')
plt.hist(time2, bins=bins[1], range=(min(time2), max(time2)), histtype='step', color='C1', label='had')
plt.hist(time4, bins=bins[2], range=(min(time4), max(time4)), histtype='step', color='C2', label='Sum')
plt.plot(bins4, fit_function_gaussian(bins4, *fitparams1), color='C3', label='Double Gaussian Fit')
# plt.plot(bins4, gaussian(bins4, *fitparams2), color='C3', label='Single Gaussian Fit')
# plt.plot(bins4, fit_function(bins4, *fitparams3), color='C3', label='Double Bi_gaussian Fit')
# plt.plot(bins4, bi_gaussian(bins4, *fitparams4), color='C3', label='Single Bi_gaussian Fit')
plt.xlim(tmin4, tmax4)
plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq1, chisqdof1, bin4, len(guesses1)), fontsize=12)
# plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq2, chisqdof2, bin4, len(guesses2)), fontsize=12)
# plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq3, chisqdof3, bin4, len(guesses3)), fontsize=12)
# plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq4, chisqdof4, bin4, len(guesses4)), fontsize=12)
plt.legend(fontsize=12)

# plt.hist(time3, bins=bins[1], range=(min(time3), max(time3)), histtype='step', color='C0', label='e-')
# plt.hist(time2, bins=bins[0], range=(min(time2), max(time2)), histtype='step', color='C1', label='had')
# plt.hist(time5, bins=bins[2], range=(min(time5), max(time5)), histtype='step', color='C2', label='Sum')
# plt.plot(bins5, fit_function_gaussian(bins4, *fitparams1), color='C3', label='Double Gaussian Fit')
# # plt.plot(bins5, gaussian(bins4, *fitparams2), color='C3', label='Single Gaussian Fit')
# # plt.plot(bins5, fit_function(bins4, *fitparams3), color='C3', label='Double Bi_gaussian Fit')
# # plt.plot(bins5, bi_gaussian(bins4, *fitparams4), color='C3', label='Single Bi_gaussian Fit')
# plt.xlim(tmin5, tmax5)
# plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq1, chisqdof1, bin5, len(guesses1)), fontsize=12)
# # plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq2, chisqdof2, bin5, len(guesses2)), fontsize=12)
# # plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq3, chisqdof3, bin5, len(guesses3)), fontsize=12)
# # plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq4, chisqdof4, bin5, len(guesses4)), fontsize=12)
# plt.legend(fontsize=12)


plt.show()
