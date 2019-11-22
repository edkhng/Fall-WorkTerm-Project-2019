import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
import timeit
import os

start = timeit.default_timer()

energy = 50

fname1 = '{} TeV Electron/e_{}TeV_merge.csv'.format(energy, energy)
fname2 = '{} TeV Had/had_{}TeV_merge.csv'.format(energy, energy)

data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')

os.mkdir('{} TeV Final Results/base_e_had_event'.format(energy))

for k in range(9):
    a = k % 3
    if k <= 2:
        b = 0
    elif k <= 5:
        b = 1
    else:
        b = 2

    fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(18,10))
    os.mkdir('{} TeV Final Results/base_e_had_event/string{}{}'.format(energy, a, b))
    for i in range(12):

        PMT_ID = [a, b, i]

        layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(data1, PMT_ID)
        layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(data2, PMT_ID)

        PMT_pos = PMT_ID_to_pos(PMT_ID)
        dx, dy, dz = seperation_vector(PMT_pos)
        d = distance_to_vertex(dx, dy, dz)
        theta = angle_to_vertex(dx, dy, dz)

        N = len(time2)/100
        time_range = get_range(time2)
        if N < 20 and time_range > 20:
            continue
        bin_size = get_bin_size(time_range, N)

        bins = int((max(time2) - min(time2))/bin_size)

        range_a = (min(time2), max(time2))

        if i <= 3:
            l = 0
        elif i <= 7:
            l = 1
        else:
            l = 2

        n1, bins1 = np.histogram(time1, bins=bins, range=range_a)
        n2, bins2 = np.histogram(time2, bins=bins, range=range_a)

        n1 = n1/100
        n2 = n2/100

        ax[l][i % 4].step(bins1[:-1], n1, color='C0', label='e-')
        ax[l][i % 4].step(bins1[:-1], n2, color='C1', label='had')
        ax[l][i % 4].step(bins1[:-1], (n1+n2), color='C2', label='Sum')
        tmin, tmax = get_range_plot(n2, bins2)
        ax[l][i % 4].set_xlim(tmin, tmax)
        ax[l][i % 4].set_ylim(0)
        if i == 0:
            ax[l][i % 4].legend(fontsize=11)
        ax[l][i % 4].set_title('PMT_ID=[{},{},{}], distance={:.2f} m\nangle={:.2f} degrees, bin_size={} ns'
                   .format(*PMT_ID, d, theta, bin_size), fontsize=11)

        header = 'e_n, had_n, bins'
        np.savetxt('{} TeV Final Results/base_e_had_event/string{}{}/{}TeV_e_had_merge_bins_PMT_ID{}{}{}.csv'.format(energy, a, b, energy, a, b, i), np.transpose([n1,n2,bins2[:-1]]), fmt='%.4f', delimiter=',', header=header, comments='#')

    fig.suptitle('Time Residuals {} Average of all e- and had; String_{}{}'.format(energy, a,b), fontsize=14)
    fig.subplots_adjust(top=0.90, bottom=0.04, left=0.06, right=0.93, hspace=0.3, wspace=0.2)
    fig.savefig('{} TeV Final Results/base_e_had_event/string{}{}/{}TeV_e_had_merge_string{}{}.png'.format(energy, a, b, energy, a, b))

stop = timeit.default_timer()
print('runtime=', stop-start)
