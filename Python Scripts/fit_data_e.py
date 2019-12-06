#!/usr/bin/env python
"""
This script fits the time residuals of the e- events for the given
sensors with a double bifrucated gaussian and also a single bifrucated
gaussian. The user will need to input the PMT_IDs of the sensors they want
fitted. The script produces a plot of the time residuals, and a plot of each
fit for the e- event. It also records the fitting parameters, chi^2
value and p-values in a csv. The script takes two parameters
the energy and event ID.
"""

import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import sys
import os

energy = int(sys.argv[1])
event_ID = int(sys.argv[2])

fname1 = '{} TeV Tau/tau_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, event_ID)
fname2 = '{} TeV Tau/{}_TeV_tau_decay_times.csv'.format(energy, energy)

data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')

hits = data2[:,1]
hit = hits[event_ID]
decay_times = data2[:,2]
decay_time = decay_times[event_ID]

PMT_IDs = [[1, 1, 0], [1,1,1], [1,1,2], [1,1,5], [1,1,6], [1,1,7], [1,0,0], [1,0,1], [1,0,2], [1,0,3], [1,0,6], [0,1,0], [0,1,1], [0,1,2], [0,1,3], [0,1,6], [0,0,0], [0,0,1], [0,0,2], [0,0,3], [0,0,4]]

for PMT_ID in PMT_IDs:

    a = PMT_ID[0]
    b = PMT_ID[1]

    fname3 = '{} TeV Final Results/base_e_had_event/string{}{}/{}TeV_e_had_merge_bins_PMT_ID{}{}{}.csv'.format(energy, a, b, energy, *PMT_ID)
    try:
        os.mkdir('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/final_plots'.format(energy, event_ID, decay_time, hit, a, b))
    except FileExistsError:
        pass
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

    PMT_pos = PMT_ID_to_pos(PMT_ID)
    dx, dy, dz = seperation_vector(PMT_pos)
    d = distance_to_vertex(dx, dy, dz)
    theta = angle_to_vertex(dx, dy, dz)

    n2 = had_n
    n3 = e_n
    n5 = e_n + had_n

    peak_height2, peak_pos2 = get_peak(n2, bins)
    peak_height3, peak_pos3 = get_peak(n3, bins)
    peak_height5, peak_pos5 = get_peak(n5, bins)

    FWHM2 = get_FWHM(n2, bins)
    FWHM3 = get_FWHM(n3, bins)
    FWHM5 = get_FWHM(n5, bins)

    cutoff = 4

    tmin5, tmax5 = get_range_fit(n5, bins, peak_height5, cutoff)

    bins_mid = bins[:-1] + bin_size/2

    range5 = (bins >= tmin5) & (bins <= tmax5)

    bins5 = bins[range5]

    bins5_mid = bins5[:-1] + bin_size/2

    if len(bins5_mid) < 9:
        continue

    range5 = (bins_mid >= tmin5) & (bins_mid <= tmax5)

    n5 = n5[range5]

    guesses2 = [peak_pos3, peak_pos2, FWHM3, FWHM2, peak_height3, peak_height2, 0.3, 0.3]
    guesses4 = [peak_pos5, FWHM5, peak_height5, 0.3]

    sigma5 = np.sqrt(n5)

    fitparams2, fitcov2 = curve_fit(fit_function, bins5_mid, n5, sigma=sigma5, bounds=(0., np.inf), p0=guesses2, maxfev=100000)
    fitparams4, fitcov4 = curve_fit(bi_gaussian, bins5_mid, n5, sigma=sigma5, bounds=(0., np.inf), p0=guesses4, maxfev=100000)

    fitparams2_error = np.sqrt(np.diag(fitcov2))
    fitparams4_error = np.sqrt(np.diag(fitcov4))

    dof2 = len(bins5_mid) - len(guesses2)
    dof4 = len(bins5_mid) - len(guesses4)

    # dx = bin_size/100
    # f_exp2 = []
    # f_exp4 = []

    #
    # for i in range(len(bins5)-1):
    #     x = np.arange(bins5[i], bins5[i+1], dx)
    #     y2 = fit_function(x, *fitparams2)
    #     y4 = bi_gaussian(x, *fitparams4)
    #     area2 = np.trapz(y2, x)
    #     area4 = np.trapz(y4, x)
    #     f_exp2.append(area2)
    #     f_exp4.append(area4)

    chisq2_bincenter, p2_bincenter = chisquare(n5, fit_function(bins5_mid, *fitparams2), ddof=len(guesses2)-1)
    chisq4_bincenter, p4_bincenter = chisquare(n5, bi_gaussian(bins5_mid, *fitparams4), ddof=len(guesses4)-1)

    # chisq2_area, p2_area = chisquare(n5, f_exp2, ddof=len(guesses2)-1)
    # chisq4_area, p4_area = chisquare(n5, f_exp4, ddof=len(guesses4)-1)

    chisqdof2 = chisq2_bincenter/dof2
    chisqdof4 = chisq4_bincenter/dof4

    pos2a, pos2b, wid2a, wid2b, amp2a, amp2b, r2a, r2b = fitparams2

    f = open('{} TeV Final Results/fitting_params_double_bi-gaussian_e.csv'.format(energy), 'a')
    # f.write('# eventID, decay_time, hits, layerID, columnID, cellID, pos1, pos2, wid1, wid2, amp1, amp2, r1, r2, upos1, upos2, uwid1, uwid2, uamp1, uamp2, ur1, ur2, chi^2, bins, params, p-value\n')
    f.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(event_ID, decay_time, hit, *PMT_ID, *fitparams2, *fitparams2_error, chisq2_bincenter, len(bins5_mid), len(guesses2), p2_bincenter))
    f.close()

    f = open('{} TeV Final Results/fitting_params_single_bi-gaussian_e.csv'.format(energy), 'a')
    # f.write('# eventID, decay_time, hits, layerID, columnID, cellID, pos, wid, amp, r, upos, uwid, uamp, ur, chi^2, bins, params, p-value\n')
    f.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(event_ID, decay_time, hit, *PMT_ID, *fitparams4, *fitparams4_error, chisq4_bincenter, len(bins5_mid), len(guesses4), p4_bincenter))
    f.close()

    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=(11,7), gridspec_kw={'height_ratios': [2, 1]})
    fig.suptitle('{} TeV Residuals with Bi_gaussian Fits; tau_decay_time = {:.2f} ns\nevent_ID: {}, PMT_ID = [{}, {}, {}], distance = {:.4f} m'.format(energy, decay_time, event_ID, *PMT_ID, d), fontsize=14)

    ax0.step(bins5[:-1], n5, where='post', color='C2', label='Data')
    # ax1.step(bins5[:-1], f_exp2, where='post', color='C1', label='Area under curve')
    ax0.step(bins5[:-1], fit_function(bins5_mid, *fitparams2), where='post', color='C3', label='Fit')
    ax0.set_title('Double Bi_gaussian Fit e-', fontsize=11)
    ax0.legend(fontsize=10)
    ax0.set_xlim(tmin5, tmax5)

    ax2.plot([tmin5,tmax5], [0,0], color='C2')
    # ax3.step(bins5[:-1], n5 - f_exp2, where='post', color='C1', label='Area under curve')
    ax2.step(bins5[:-1], n5 - fit_function(bins5_mid, *fitparams2), where='post', color='C3')
    ax2.set_title('chi^2 = {:.2f}, bins = {}, p = {:.2e}'.format(chisq2_bincenter, len(bins5_mid), p2_bincenter), fontsize=11)
    ax2.set_xlim(tmin5, tmax5)



    ax1.step(bins5[:-1], n5, where='post', color='C2', label='Data')
    # ax5.step(bins5[:-1], f_exp4, where='post', color='C1', label='Area under curve')
    ax1.step(bins5[:-1], bi_gaussian(bins5_mid, *fitparams4), where='post', color='C3', label='Fit')
    ax1.set_title('Single Bi_gaussian Fit e-', fontsize=11)
    ax1.legend(fontsize=10)
    ax1.set_xlim(tmin5, tmax5)

    ax3.plot([tmin5,tmax5], [0,0], color='C2')
    # ax7.step(bins5[:-1], n5 - f_exp4, where='post', color='C1', label='Area under curve')
    ax3.step(bins5[:-1], n5 - bi_gaussian(bins5_mid, *fitparams4), where='post', color='C3')
    ax3.set_title('chi^2 = {:.2f}, bins = {}, p = {:.2e}'.format(chisq4_bincenter, len(bins5_mid), p4_bincenter), fontsize=11)
    ax3.set_xlim(tmin5, tmax5)

    fig.subplots_adjust(top=0.87, bottom=0.08, left=0.08, right=0.96, hspace=0.27, wspace=0.15)
    plt.savefig('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/final_plots/bi_gaussian_fit_PMT_ID{}{}{}.png'.format(energy, event_ID, decay_time, hit, a, b, *PMT_ID))
    plt.close()



    n5 = np.insert(n5, 0, n5[0])

    plt.figure(figsize=(7,9))
    plt.suptitle('{} TeV Time Residuals with Bi_gaussian Fits; tau_decay_time = {:.2f} ns\nevent_ID: {}, PMT_ID = [{}, {}, {}], distance = {:.4f} m'.format(energy, decay_time, event_ID, *PMT_ID, d), fontsize=12)


    plt.subplot(211)
    plt.plot(bins5_mid, bi_gaussian(bins5_mid, pos2a, wid2a, amp2a, r2a), label='Bi_gaussian 1')
    plt.plot(bins5_mid, bi_gaussian(bins5_mid, pos2b, wid2b, amp2b, r2b), label='Bi_gaussian 2')
    plt.step(bins5, n5, color='C2', label='Data')
    plt.plot(bins5_mid, fit_function(bins5_mid, *fitparams2), color='C3', label='Fit')
    plt.xlim(tmin5, tmax5)
    plt.ylim(0)
    plt.xlabel('time [ns]')
    plt.ylabel('bin count')
    plt.title('Double Bi_gaussian Fit e-\nchi^2 = {:.2f}, bins = {}, p = {:.2e}'.format(chisq2_bincenter, len(bins5_mid), p2_bincenter), fontsize=11)
    plt.legend(fontsize=10)



    plt.subplot(212)
    plt.step(bins5, n5, color='C2', label='Sum')
    plt.plot(bins5_mid, bi_gaussian(bins5_mid, *fitparams4), color='C3', label='Fit')
    plt.xlim(tmin5, tmax5)
    plt.ylim(0)
    plt.xlabel('time [ns]')
    plt.ylabel('bin count')
    plt.title('Single Bi_gaussian Fit e-\nchi^2 = {:.2f}, bins = {}, p = {:.2e}'.format(chisq4_bincenter, len(bins5_mid), p4_bincenter), fontsize=11)
    plt.legend(fontsize=10)

    plt.subplots_adjust(top=0.87, bottom=0.08, left=0.12, right=0.90, hspace=0.33)
    plt.savefig('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/final_plots/bi_gaussian_curve_fit_PMT_ID{}{}{}.png'.format(energy, event_ID, decay_time, hit, a, b, *PMT_ID))
    plt.close()

    # plt.show()
