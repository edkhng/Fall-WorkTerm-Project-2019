#!/usr/bin/env python
"""
This script goes through every single sensor and plots the time residuals
for both the tau and e- events. The plots are grouped by strings, where
a single plot has the time residuals of all the sensors on that string.
The script also records the number of hits by each component, the approximate
time range, peak height, location and FWHM in a text file. The script
takes two arguments, the energy and the event ID. 
"""

import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
import timeit
import os
import sys


start = timeit.default_timer()

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

os.mkdir('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}'.format(energy, event_ID, decay_time, hit))

for k in range(9):
    a = k % 3
    if k <= 2:
        b = 0
    elif k <= 5:
        b = 1
    else:
        b = 2

    fig1, ax1 = plt.subplots(nrows=3, ncols=4, figsize=(18,10))
    fig2, ax2 = plt.subplots(nrows=3, ncols=4, figsize=(18,10))

    os.mkdir('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}'.format(energy, event_ID, decay_time, hit, a, b))

    f = open('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/string{}{}_output.txt'.format(energy, event_ID, decay_time, hit, a, b, a, b), 'w+')
    f.write("Decay time = {} ns\n".format(decay_time))
    f.close()
    for i in range(12):

        PMT_ID = [a, b, i]

        fname3 = '{} TeV Final Results/base_e_had_event/string{}{}/{}TeV_e_had_merge_bins_PMT_ID{}{}{}.csv'.format(energy, a, b, energy, *PMT_ID)

        try:
            data3 = np.loadtxt(fname3, delimiter=',', comments='#')
        except OSError:
            continue

        e_n, had_n, bins = data3[:,0], data3[:,1], data3[:,2]
        layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(data1, PMT_ID)

        bin_size = bins[1] - bins[0]
        last_bin = bins[-1] + bin_size
        bins = np.append(bins, last_bin)

        PMT_pos = PMT_ID_to_pos(PMT_ID)
        dx, dy, dz = seperation_vector(PMT_pos)
        d = distance_to_vertex(dx, dy, dz)
        theta = angle_to_vertex(dx, dy, dz)

        if i <= 3:
            l = 0
        elif i <= 7:
            l = 1
        else:
            l = 2

        np.random.seed(event_ID)
        # e- and had are smooth averages of 100 events, alter by drawing from a poisson
        # distribution with lambda equal to n
        e_n = np.random.poisson(e_n)
        had_n = np.random.poisson(had_n)

        n1, bins1 = np.histogram(time1, bins=bins)
        ax1[l][i % 4].step(bins[:-1], n1, color='C0', label='tau')
        ax1[l][i % 4].step(bins[:-1], had_n, color='C1', label='had')
        ax1[l][i % 4].step(bins[:-1], n1+had_n, color='C2', label='Sum')
        tmin_a, tmax_a = get_range_plot(n1+had_n, bins)
        ax1[l][i % 4].set_xlim(tmin_a, tmax_a)
        ax1[l][i % 4].set_ylim(0)
        if i == 0:
            ax1[l][i % 4].legend(fontsize=11)
        ax1[l][i % 4].set_title('PMT_ID=[{},{},{}], distance={:.2f} m\nangle={:.2f} degrees, bin_size={:.0f} ns'
                   .format(*PMT_ID, d, theta, bin_size), fontsize=11)


        ax2[l][i % 4].step(bins[:-1], e_n, color='C0', label='e-')
        ax2[l][i % 4].step(bins[:-1], had_n, color='C1', label='had')
        ax2[l][i % 4].step(bins[:-1], e_n+had_n, color='C2', label='Sum')
        tmin_b, tmax_b = get_range_plot(e_n+had_n, bins)
        ax2[l][i % 4].set_xlim(tmin_b, tmax_b)
        ax2[l][i % 4].set_ylim(0)
        if i == 0:
            ax2[l][i % 4].legend(fontsize=11)
        ax2[l][i % 4].set_title('PMT_ID=[{},{},{}], distance={:.2f} m\nangle={:.2f} degrees, bin_size={:.0f} ns'
                   .format(*PMT_ID, d, theta, bin_size), fontsize=11)

        f = open('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/string{}{}_output.txt'.format(energy, event_ID, decay_time, hit, a, b, a, b), 'a')
        f.write("#####################################################################\n")
        f.write("\nPMT_ID: [{}, {}, {}]\n".format(*PMT_ID))
        f.write("Total time range tau neutrino: {:.2f}-{:.2f} ns\n".format(tmin_a, tmax_a))
        f.write("Total time range electron neutrino: {:.2f}-{:.2f} ns\n".format(tmin_b, tmax_b))

        peak_height1, peak_pos1 = get_peak(n1, bins)
        peak_height2, peak_pos2 = get_peak(had_n, bins)
        peak_height3, peak_pos3 = get_peak(e_n, bins)
        peak_height4, peak_pos4 = get_peak(had_n+n1, bins)
        peak_height5, peak_pos5 = get_peak(e_n+had_n, bins)

        FWHM1 = get_FWHM(n1, bins)
        FWHM2 = get_FWHM(had_n, bins)
        FWHM3 = get_FWHM(e_n, bins)
        FWHM4 = get_FWHM(had_n+n1, bins)
        FWHM5 = get_FWHM(e_n+had_n, bins)

        f.write("\nTau; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height1, peak_pos1, FWHM1))
        f.write("Had; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height2, peak_pos2, FWHM2))
        f.write("e-; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height3, peak_pos3, FWHM3))
        f.write("Total tau neutrino; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height4, peak_pos4, FWHM4))
        f.write("Total electron neutrino; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height5, peak_pos5, FWHM5))
        f.close()

    fig1.suptitle('Time Residuals {} TeV Tau Neutrino; String_{}{}, eventID: {}, decay_time={:.2f} ns, hits={}'.format(2*energy, a, b, event_ID, decay_time, hit), fontsize=14)
    fig2.suptitle('Time Residuals {} TeV Electron Neutrino; String_{}{}, eventID: {}, hits={}'.format(2*energy, a, b, event_ID, hit), fontsize=14)
    fig1.subplots_adjust(top=0.90, bottom=0.04, left=0.06, right=0.93, hspace=0.3, wspace=0.2)
    fig2.subplots_adjust(top=0.90, bottom=0.04, left=0.06, right=0.93, hspace=0.3, wspace=0.2)
    fig1.savefig('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/tau_string{}{}_output.png'.format(energy, event_ID, decay_time, hit, a, b, a, b))
    fig2.savefig('{} TeV Final Results/tau{}_decay_time={:.2f}_hits={:.2e}/string{}{}/e_string{}{}_output.png'.format(energy, event_ID, decay_time, hit, a, b, a, b))
    plt.close()

stop = timeit.default_timer()
print('runtime=', stop-start)
# plt.show()
