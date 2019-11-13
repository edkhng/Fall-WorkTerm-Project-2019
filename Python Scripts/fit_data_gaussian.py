import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import timeit

neutrino_energy = 10

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_{}TeV423_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/had_{}TeV234_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_{}TeV163_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)

fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)
fname5 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)

PMT_IDs = [[1, 1, 0], [1, 1, 1], [1, 1, 5], [1, 1, 7]]

for PMT_ID in PMT_IDs:

    a = PMT_ID[0]
    b = PMT_ID[1]

    layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(fname1, PMT_ID)
    layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(fname2, PMT_ID)
    layerID3, columnID3, cellID3, time3, x3, y3, z3, energy3 = get_data(fname3, PMT_ID)
    layerID4, columnID4, cellID4, time4, x4, y4, z4, energy4 = get_data(fname4, PMT_ID)
    layerID5, columnID5, cellID5, time5, x5, y5, z5, energy5 = get_data(fname5, PMT_ID)

    PMT_pos = PMT_ID_to_pos(PMT_ID)
    dx, dy, dz = seperation_vector('A', PMT_pos)
    d = distance_to_vertex(dx, dy, dz)
    theta = angle_to_vertex('A', dx, dy, dz)

    N = min(len(time1), len(time2), len(time3))
    time_range = get_range(time4)
    bin_size = get_bin_size(time_range, N)

    bins_a = int((max(time4) - min(time4))/bin_size)
    bins_b = int((max(time5) - min(time5))/bin_size)

    range_a = (min(time4), max(time4))
    range_b = (min(time5), max(time5))

    n1, bins1 = np.histogram(time1, bins=bins_a, range=range_a)
    n2a, bins2a = np.histogram(time2, bins=bins_a, range=range_a)
    n2b, bins2b = np.histogram(time2, bins=bins_b, range=range_b)
    n3, bins3 = np.histogram(time3, bins=bins_b, range=range_b)
    n4, bins4 = np.histogram(time4, bins=bins_a, range=range_a)
    n5, bins5 = np.histogram(time5, bins=bins_b, range=range_b)

    peak_height1, peak_pos1 = get_peak(n1, bins1)
    peak_height2a, peak_pos2a = get_peak(n2a, bins2a)
    peak_height2b, peak_pos2b = get_peak(n2b, bins2b)
    peak_height3, peak_pos3 = get_peak(n3, bins3)
    peak_height4, peak_pos4 = get_peak(n4, bins4)
    peak_height5, peak_pos5 = get_peak(n5, bins5)

    FWHM1 = get_FWHM(n1, bins1)
    FWHM2a = get_FWHM(n2a, bins2a)
    FWHM2b = get_FWHM(n2b, bins2b)
    FWHM3 = get_FWHM(n3, bins3)
    FWHM4 = get_FWHM(n4, bins4)
    FWHM5 = get_FWHM(n5, bins5)

    tmin4, tmax4 = get_range_fit(n4, bins4, peak_height4)
    tmin5, tmax5 = get_range_fit(n5, bins5, peak_height5)

    N4 = range_size(time4, tmin4, tmax4)
    N5 = range_size(time5, tmin5, tmax5)

    bin4 = int((tmax4-tmin4)/bin_size)
    bin5 = int((tmax5-tmin5)/bin_size)

    # print(tmin4, tmax4)
    # print(bin4, bin5)

    # print("\nThe number of bins: [{}, {}, {}, {}, {}]".format(bin1, bin2, bin3, bin4, bin5))

    n1, bins1 = np.histogram(time1, bins=bin4, range=(tmin4, tmax4))
    n2a, bins2a = np.histogram(time2, bins=bin4, range=(tmin4, tmax4))
    n2b, bins2b = np.histogram(time2, bins=bin5, range=(tmin5, tmax5))
    n3, bins3 = np.histogram(time3, bins=bin5, range=(tmin5, tmax5))
    n4, bins4 = np.histogram(time4, bins=bin4, range=(tmin4, tmax4))
    n5, bins5 = np.histogram(time5, bins=bin5, range=(tmin5, tmax5))

    guesses1 = [peak_pos1, peak_pos2a, FWHM1, FWHM2a, peak_height1, peak_height2a]
    guesses2 = [peak_pos3, peak_pos2b, FWHM3, FWHM2b, peak_height3, peak_height2b]
    guesses3 = [peak_pos4, FWHM4, peak_height4]
    guesses4 = [peak_pos5, FWHM5, peak_height5]

    fitparams1, fitcov1 = curve_fit(fit_function_gaussian, bins4[:-1], n4, bounds=(0., np.inf), p0=guesses1, maxfev=100000)
    fitparams2, fitcov2 = curve_fit(fit_function_gaussian, bins5[:-1], n5, bounds=(0., np.inf), p0=guesses2, maxfev=100000)
    fitparams3, fitcov3 = curve_fit(gaussian, bins4[:-1], n4, bounds=(0., np.inf), p0=guesses3, maxfev=100000)
    fitparams4, fitcov4 = curve_fit(gaussian, bins5[:-1], n5, bounds=(0., np.inf), p0=guesses4, maxfev=100000)

    # fitparams1_error = np.sqrt(np.diag(fitcov1))
    # fitparams2_error = np.sqrt(np.diag(fitcov2))

    dof1 = bin4 - len(guesses1)
    dof2 = bin5 - len(guesses2)
    dof3 = bin4 - len(guesses3)
    dof4 = bin5 - len(guesses4)

    chisq1, p1 = chisquare(n4, fit_function_gaussian(bins4[:-1], *fitparams1), ddof=dof1)
    chisq3, p3 = chisquare(n5, fit_function_gaussian(bins5[:-1], *fitparams2), ddof=dof2)
    chisq2, p2 = chisquare(n4, gaussian(bins4[:-1], *fitparams3), ddof=dof3)
    chisq4, p4 = chisquare(n5, gaussian(bins5[:-1], *fitparams4), ddof=dof4)

    chisqdof1 = chisq1/dof1
    chisqdof2 = chisq2/dof2
    chisqdof3 = chisq3/dof3
    chisqdof4 = chisq4/dof4

    pos1a, pos1b, wid1a, wid1b, amp1a, amp1b = fitparams1
    pos2a, pos2b, wid2a, wid2b, amp2a, amp2b = fitparams2

    f = open('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/fitting_params_gaussian_tau_all.csv', 'a+')
    # if PMT_ID == [0,0,0]:
    #     f.write('# energy, PMT_ID, pos1, pos2, wid1, wid2, amp1, amp2, chi^2, bins, params\n')
    f.write('{}, {}-{}-{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(neutrino_energy, *PMT_ID, *fitparams1, chisq1, bin4, len(guesses1)))
    # f.write('# pos1, pos2, wid1, wid2, amp1, amp2, upos1, upos2, uwid1, uwid2, uamp1, uamp2')
    # f.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}'.format(*fitparams1, *fitparams1_error))

    f = open('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/fitting_params_gaussian_e_all.csv', 'a+')
    # if PMT_ID == [0,0,0]:
    #     f.write('# energy, PMT_ID, pos1, pos2, wid1, wid2, amp1, amp2, chi^2, bins, params\n')
    f.write('{}, {}-{}-{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(neutrino_energy, *PMT_ID, *fitparams2, chisq2, bin5, len(guesses2)))
    # f.write('# pos1, pos2, wid1, wid2, amp1, amp2, upos1, upos2, uwid1, uwid2, uamp1, uamp2')
    # f.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}'.format(*fitparams2, *fitparams2_error))

    #
    # f = open('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/fitting_params.txt'.format(neutrino_energy), 'a+')
    # f.write("#####################################################################\n")
    # f.write("\nGaussian Fits; PMT_ID: [{}, {}, {}]\n".format(*PMT_ID))
    # f.write("\n{} TeV Tau Neutrino\n".format(2*neutrino_energy))
    # f.write("Gaussian 1:\npos1 = {:.2f}\nwid1 = {:.2f}\namp1 = {:.2f}\n".format(pos1a, wid1a, amp1a))
    # f.write("Gaussian 2:\npos2 = {:.2f}\nwid2 = {:.2f}\namp2 = {:.2f}\n".format(pos1b, wid1b, amp1b))
    # f.write("\n{} TeV Electron Neutrino\n".format(2*neutrino_energy))
    # f.write("Gaussian 1:\npos1 = {:.2f}\nwid1 = {:.2f}\namp1 = {:.2f}\n".format(pos2a, wid2a, amp2a))
    # f.write("Gaussian 2:\npos2 = {:.2f}\nwid2 = {:.2f}\namp2 = {:.2f}\n".format(pos2b, wid2b, amp2b))
    # f.close()

    plt.figure(figsize=(12,10))
    plt.suptitle('{} TeV Time Residuals with Gaussian Fits;\nPMT_ID = [{}, {}, {}], distance = {:.4f} m, angle={:.4f} degrees'.format(neutrino_energy, *PMT_ID, d, theta), fontsize=14)

    plt.subplot(321)
    plt.hist(time1, bins=bins_a, range=range_a, histtype='step', color='C0', label='tau')
    plt.hist(time2, bins=bins_a, range=range_a, histtype='step', color='C1', label='had')
    plt.hist(time4, bins=bins_a, range=range_a, histtype='step', color='C2', label='Sum')
    tmin_a, tmax_a = get_range_plot(n4, bins4)
    plt.xlim(tmin_a, tmax_a)
    plt.legend(fontsize=10)

    plt.subplot(322)
    plt.hist(time3, bins=bins_b, range=range_b, histtype='step', color='C0', label='tau')
    plt.hist(time2, bins=bins_b, range=range_b, histtype='step', color='C1', label='had')
    plt.hist(time5, bins=bins_b, range=range_b, histtype='step', color='C2', label='Sum')
    tmin_b, tmax_b = get_range_plot(n5, bins5)
    plt.xlim(tmin_b, tmax_b)
    plt.legend(fontsize=10)

    plt.subplot(323)
    plt.plot(bins4, gaussian(bins4, pos1a, wid1a, amp1a), label='Gaussian 1')
    plt.plot(bins4, gaussian(bins4, pos1b, wid1b, amp1b), label='Gaussian 2')
    plt.hist(time4, bins=bins_a, range=range_a, histtype='step', color='C2', label='Sum')
    plt.plot(bins4, fit_function_gaussian(bins4, *fitparams1), color='C3', label='Double Gaussian Fit')
    plt.xlim(tmin4, tmax4)
    plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq1, chisqdof1, bin4, len(guesses1)), fontsize=11)
    plt.legend(fontsize=10)

    plt.subplot(324)
    plt.plot(bins5, gaussian(bins5, pos2a, wid2a, amp2a), label='Gaussian 1')
    plt.plot(bins5, gaussian(bins5, pos2b, wid2b, amp2b), label='Gaussian 2')
    plt.hist(time5, bins=bins_b, range=range_b, histtype='step', color='C2', label='Sum')
    plt.plot(bins5, fit_function_gaussian(bins5, *fitparams2), color='C3', label='Double Gaussian Fit')
    plt.xlim(tmin5, tmax5)
    plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq3, chisqdof3, bin5, len(guesses2)), fontsize=11)
    plt.legend(fontsize=10)

    plt.subplot(325)
    plt.hist(time4, bins=bins_a, range=range_a, histtype='step', color='C2', label='Sum')
    plt.plot(bins4, gaussian(bins4, *fitparams3), color='C3', label='Single Gaussian Fit')
    plt.xlim(tmin4, tmax4)
    plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq2, chisqdof2, bin4, len(guesses3)), fontsize=11)
    plt.legend(fontsize=10)

    plt.subplot(326)
    plt.hist(time5, bins=bins_b, range=range_b, histtype='step', color='C2', label='Sum')
    plt.plot(bins5, gaussian(bins5, *fitparams4), color='C3', label='Single Gaussian Fit')
    plt.xlim(tmin5, tmax5)
    plt.title('chi^2 = {:.2f}, chi^2/dof = {:.2f}, bins = {}, params = {}'.format(chisq4, chisqdof4, bin5, len(guesses4)), fontsize=11)
    plt.legend(fontsize=10)

    plt.subplots_adjust(top=0.90, bottom=0.06, left=0.08, right=0.93, hspace=0.35, wspace=0.15)
    plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/{}_TeV_string{}{}/{}_TeV_PMT_ID_{}{}{}_gaussian_fits.png'.format(neutrino_energy, neutrino_energy, a, b, neutrino_energy, *PMT_ID))
# plt.show()
