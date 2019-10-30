import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare

neutrino_energy = 100

fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_{}TeV377_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/had_{}TeV411_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
# fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_{}TeV534_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)
# fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)

PMT_ID = [1, 1, 1]
layerID, columnID, cellID, time, x, y, z, energy = get_data(fname, PMT_ID)
print(len(x))

layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(fname2, PMT_ID)
layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(fname1, PMT_ID)


PMT_pos = PMT_ID_to_pos(PMT_ID)
d = distance_to_vertex('A', PMT_pos)
theta = angle_to_vertex('A', PMT_pos)


# n, bins, patches = plt.hist(time1, bins=30, range=(104,134))
# # fitparams, fitcov = curve_fit(gaussian, bins[:-1], n, p0=[115, 15, 530])
# # fitparams, fitcov = curve_fit(bi_gaussian, bins[:-1], n, p0=[115, 15, 530, 0.3])
# fitparams, fitcov = curve_fit(fit_function_gaussian, bins[:-1], n, p0=[110, 125, 15, 20, 850, 400])
# # fitparams, fitcov = curve_fit(fit_function, bins[:-1], n, p0=[110, 125, 15, 20, 850, 400, 0.3, 0.3])
# plt.close()
#
# # chisq, p = chisquare(n, gaussian(bins[:-1], *fitparams), ddof=3)
# # chisq, p = chisquare(n, bi_gaussian(bins[:-1], *fitparams), ddof=4)
# chisq, p = chisquare(n, fit_function_gaussian(bins[:-1], *fitparams), ddof=6)
# # chisq, p = chisquare(n, fit_function(bins[:-1], *fitparams), ddof=8)
# # chisqdof = chisq/3
# # chisqdof = chisq/4
# chisqdof = chisq/6
# # chisqdof = chisq/8
# print(fitparams)
# print(p)
select = time < 1000
time = time[select]
print(max(time), min(time))
print(np.std(time), np.mean(time))
print(get_peak(time, 1))
print(get_FWHM(time, 1))


# plt.figure(figsize=(10,7))
# plt.hist(time, bins=30, range=(104,134), histtype='step', color='C1', label='{} TeV tau'.format(neutrino_energy))
# # plt.hist(time2, bins=30, range=(104,134), histtype='step', color='C0', label='{} TeV had'.format(neutrino_energy))
# # plt.hist(time1, bins=30, range=(104,134), histtype='step', color='C2', label='Sum')
# #
# # plt.plot(bins, gaussian(bins, fitparams[0], fitparams[2], fitparams[4]))
# # plt.plot(bins, gaussian(bins, fitparams[1], fitparams[3], fitparams[5]))
#
# # plt.plot(bins, bi_gaussian(bins, fitparams[1], fitparams[3], fitparams[5], fitparams[7]))
# # plt.plot(bins, bi_gaussian(bins, fitparams[0], fitparams[2], fitparams[4], fitparams[6]))
#
# # plt.plot(bins, gaussian(bins, *fitparams), label='Gaussian Fit, pos={:.2f}\nFWHM={:.2f}, Amp={:.2f}'.format(*fitparams))
# # plt.plot(bins, bi_gaussian(bins, *fitparams), label='Bi_gaussian Fit, pos={:.2f}\nFWHM={:.2f}, Amp={:.2f}\nr={:.2f}'.format(*fitparams))
# # plt.plot(bins, fit_function_gaussian(bins, *fitparams))
# # plt.plot(bins, fit_function(bins, *fitparams))#, label='Bi_gaussian Fit, pos={:.2f}\nFWHM={:.2f}, Amp={:.2f}\nr={:.2f}'.format(*fitparams))
#
# # plt.title('PMT_ID = [{}, {}, {}], distance = {:.4f} m\nangle={:.4f} degrees, chi^2/dof = {:.2f}'
# #            .format(PMT_ID[0], PMT_ID[1], PMT_ID[2], d, theta, chisqdof), fontsize=16)
# plt.title('PMT_ID = [{}, {}, {}], distance = {:.4f} m\nangle={:.4f} degrees'
#            .format(PMT_ID[0], PMT_ID[1], PMT_ID[2], d, theta), fontsize=16)
# plt.xlabel("Time [ns]", fontsize=14)
# plt.ylabel("Photon Count", fontsize=14)
# plt.legend(fontsize=14)
# # plt.xlim(104,134)
# #plt.savefig('tau_neutrino_100TeV.png')
# plt.show()
