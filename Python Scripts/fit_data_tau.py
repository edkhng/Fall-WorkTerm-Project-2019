#!/usr/bin/env python
"""
This script fits the time residuals of the tau events for the given
sensors with a double bifrucated gaussian and also a single bifrucated
gaussian. The user will need to input the PMT_IDs of the sensors they want
fitted. The script produces a plot of the time residuals, and a plot of each
fit for the tau event. It also records the fitting parameters, chi^2
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

    had_n, bins = data3[:,1], data3[:,2]
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
    had_n = np.random.poisson(had_n)

    PMT_pos = PMT_ID_to_pos(PMT_ID)
    dx, dy, dz = seperation_vector(PMT_pos)
    d = distance_to_vertex(dx, dy, dz)
    theta = angle_to_vertex(dx, dy, dz)

    n1, bins1 = np.histogram(time1, bins=bins)
    n2 = had_n
    n4 = n1 + had_n

    peak_height1, peak_pos1 = get_peak(n1, bins)
    peak_height2, peak_pos2 = get_peak(n2, bins)
    peak_height4, peak_pos4 = get_peak(n4, bins)

    FWHM1 = get_FWHM(n1, bins)
    FWHM2 = get_FWHM(n2, bins)
    FWHM4 = get_FWHM(n4, bins)

    cutoff = 4

    tmin4, tmax4 = get_range_fit(n4, bins, peak_height4, cutoff)

    bins_mid = bins[:-1] + bin_size/2

    range4 = (bins >= tmin4) & (bins <= tmax4)

    bins4 = bins[range4]

    bins4_mid = bins4[:-1] + bin_size/2

    if len(bins4_mid) < 9:
        continue

    range4 = (bins_mid >= tmin4) & (bins_mid <= tmax4)

    n4 = n4[range4]

    guesses1 = [peak_pos1, peak_pos2, FWHM1, FWHM2, peak_height1, peak_height2, 0.3, 0.3]
    guesses3 = [peak_pos4, FWHM4, peak_height4, 0.3]

    sigma4 = np.sqrt(n4)

    fitparams1, fitcov1 = curve_fit(fit_function, bins4_mid, n4, sigma=sigma4, bounds=(0., np.inf), p0=guesses1, maxfev=100000)
    fitparams3, fitcov3 = curve_fit(bi_gaussian, bins4_mid, n4, sigma=sigma4, bounds=(0., np.inf), p0=guesses3, maxfev=100000)

    fitparams1_error = np.sqrt(np.diag(fitcov1))
    fitparams3_error = np.sqrt(np.diag(fitcov3))

    dof1 = len(bins4_mid) - len(guesses1)
    dof3 = len(bins4_mid) - len(guesses3)

    # dx = bin_size/100
    # f_exp1 = []
    # f_exp3 = []
    # for i in range(len(bins4)-1):
    #     x = np.arange(bins4[i], bins4[i+1], dx)
    #     y1 = fit_function(x, *fitparams1)
    #     y3 = bi_gaussian(x, *fitparams3)
    #     area1 = np.trapz(y1, x)
    #     area3 = np.trapz(y3, x)
    #     f_exp1.append(area1)
    #     f_exp3.append(area3)
    #


    chisq1_bincenter, p1_bincenter = chisquare(n4, fit_function(bins4_mid, *fitparams1), ddof=len(guesses1)-1)
    chisq3_bincenter, p3_bincenter = chisquare(n4, bi_gaussian(bins4_mid, *fitparams3), ddof=len(guesses3)-1)

    # chisq1_area, p1_area = chisquare(n4, f_exp1, ddof=len(guesses1)-1)
    # chisq2_area, p2_area = chisquare(n5, f_exp2, ddof=len(guesses2)-1)
    # chisq3_area, p3_area = chisquare(n4, f_exp3, ddof=len(guesses3)-1)
    # chisq4_area, p4_area = chisquare(n5, f_exp4, ddof=len(guesses4)-1)

    chisqdof1 = chisq1_bincenter/dof1
    chisqdof3 = chisq3_bincenter/dof3

    pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b = fitparams1

    f = open('{} TeV Final Results/fitting_params_double_bi-gaussian_tau.csv'.format(energy), 'a')
    # f.write('# eventID, decay_time, hits, layerID, columnID, cellID, pos1, pos2, wid1, wid2, amp1, amp2, r1, r2, upos1, upos2, uwid1, uwid2, uamp1, uamp2, ur1, ur2, chi^2, bins, params, p-value\n')
    f.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(event_ID, decay_time, hit, *PMT_ID, *fitparams1, *fitparams1_error, chisq1_bincenter, len(bins4_mid), len(guesses1), p1_bincenter))
    f.close()

    f = open('{} TeV Final Results/fitting_params_single_bi-gaussian_tau.csv'.format(energy), 'a')
    # f.write('# eventID, decay_time, hits, layerID, columnID, cellID, pos, wid, amp, r, upos, uwid, uamp, ur, chi^2, bins, params, p-value\n')
    f.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(event_ID, decay_time, hit, *PMT_ID, *fitparams3, *fitparams3_error, chisq3_bincenter, len(bins4_mid), len(guesses3), p3_bincenter))
    f.close()

    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=(11,7), gridspec_kw={'height_ratios': [2, 1]})
    fig.suptitle('{} TeV Residuals with Bi_gaussian Fits; tau_decay_time = {:.2f} ns\nevent_ID: {}, PMT_ID = [{}, {}, {}], distance = {:.4f} m'.format(energy, decay_time, event_ID, *PMT_ID, d), fontsize=14)

    ax0.step(bins4[:-1], n4, where='post', color='C2', label='Data')
    # ax0.step(bins4[:-1], f_exp1, where='post', color='C1', label='Area under curve')
    ax0.step(bins4[:-1], fit_function(bins4_mid, *fitparams1), where='post', color='C3', label='Fit')
    ax0.set_title('Double Bi_gaussian Fit Tau', fontsize=11)
    ax0.legend(fontsize=10)
    ax0.set_xlim(tmin4, tmax4)
    ax0.set_ylabel('bincount')

    ax2.plot([tmin4,tmax4], [0,0], color='C2')
    # ax2.step(bins4[:-1], n4 - f_exp1, where='post', color='C1', label='Area under curve')
    ax2.step(bins4[:-1], n4 - fit_function(bins4_mid, *fitparams1), where='post', color='C3')
    ax2.set_title('chi^2 = {:.2f}, bins = {}, p = {:.2e}'.format(chisq1_bincenter, len(bins4_mid), p1_bincenter), fontsize=11)
    ax2.set_xlim(tmin4, tmax4)
    ax2.set_xlabel('time [ns]')

    ax1.step(bins4[:-1], n4, where='post', color='C2', label='Data')
    # ax4.step(bins4[:-1], f_exp3, where='post', color='C1', label='Area under curve')
    ax1.step(bins4[:-1], bi_gaussian(bins4_mid, *fitparams3), where='post', color='C3', label='Fit')
    ax1.set_title('Single Bi_gaussian Fit Tau', fontsize=11)
    ax1.legend(fontsize=10)
    ax1.set_xlim(tmin4, tmax4)
    ax1.set_ylabel('bincount')

    ax3.plot([tmin4,tmax4], [0,0], color='C2')
    # ax6.step(bins4[:-1], n4 - f_exp3, where='post', color='C1', label='Area under curve')
    ax3.step(bins4[:-1], n4 - bi_gaussian(bins4_mid, *fitparams3), where='post', color='C3')
    ax3.set_title('chi^2 = {:.2f}, bins = {}, p = {:.2e}'.format(chisq3_bincenter, len(bins4_mid), p3_bincenter), fontsize=11)
    ax3.set_xlim(tmin4, tmax4)
    ax3.set_xlabel('time [ns]')

    fig.subplots_adjust(top=0.87, bottom=0.08, left=0.08, right=0.96, hspace=0.27, wspace=0.15)
    plt.savefig('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/final_plots/bi_gaussian_fit_PMT_ID{}{}{}.png'.format(energy, event_ID, decay_time, hit, a, b, *PMT_ID))
    plt.close()

    n4 = np.insert(n4, 0, n4[0])

    plt.figure(figsize=(7,9))
    plt.suptitle('{} TeV Time Residuals with Bi_gaussian Fits; tau_decay_time = {:.2f} ns\nevent_ID: {}, PMT_ID = [{}, {}, {}], distance = {:.4f} m'.format(energy, decay_time, event_ID, *PMT_ID, d), fontsize=12)

    plt.subplot(211)
    plt.plot(bins4_mid, bi_gaussian(bins4_mid, pos1a, wid1a, amp1a, r1a), label='Bi_gaussian 1')
    plt.plot(bins4_mid, bi_gaussian(bins4_mid, pos1b, wid1b, amp1b, r1b), label='Bi_gaussian 2')
    plt.step(bins4, n4, color='C2', label='Data')
    plt.plot(bins4_mid, fit_function(bins4_mid, *fitparams1), color='C3', label='Fit')
    plt.xlim(tmin4, tmax4)
    plt.xlabel('time [ns]')
    plt.ylabel('bin count')
    plt.ylim(0)
    plt.title('Double Bi_gaussian Fit Tau-\nchi^2 = {:.2f}, bins = {}, p = {:.2e}'.format(chisq1_bincenter, len(bins4_mid), p1_bincenter), fontsize=11)
    plt.legend(fontsize=10)

    plt.subplot(212)
    plt.step(bins4, n4, color='C2', label='Sum')
    plt.plot(bins4_mid, bi_gaussian(bins4_mid, *fitparams3), color='C3', label='Fit')
    plt.xlim(tmin4, tmax4)
    plt.ylim(0)
    plt.xlabel('time [ns]')
    plt.ylabel('bin count')
    plt.title('Single Bi_gaussian Fit Tau\nchi^2 = {:.2f}, bins = {}, p = {:.2e}'.format(chisq3_bincenter, len(bins4_mid), p3_bincenter), fontsize=11)
    plt.legend(fontsize=10)


    plt.subplots_adjust(top=0.87, bottom=0.08, left=0.12, right=0.90, hspace=0.33)
    plt.savefig('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/final_plots/bi_gaussian_curve_fit_PMT_ID{}{}{}.png'.format(energy, event_ID, decay_time, hit, a, b, *PMT_ID))
    plt.close()

    # plt.show()
