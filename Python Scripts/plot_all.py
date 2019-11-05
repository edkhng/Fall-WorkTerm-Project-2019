import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import os
import timeit

start = timeit.default_timer()

neutrino_energy = 100

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_{}TeV377_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/had_{}TeV411_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_{}TeV371_nt_Ntuple.csv'.format(neutrino_energy, neutrino_energy)

fname4 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)
fname5 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_had_merge_{}_TeV.csv'.format(neutrino_energy, neutrino_energy)

# for k in range(9):
#     a = k % 3
#     if k <= 2:
#         b = 0
#     elif k <= 5:
#         b = 1
#     else:
#         b = 2

a = 1
b = 1

#plt.figure(figsize=(14,10))
fig1, ax1 = plt.subplots(nrows=3, ncols=4, figsize=(18,10))
fig2, ax2 = plt.subplots(nrows=3, ncols=4, figsize=(18,10))

# os.mkdir('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/string{}{}'.format(neutrino_energy, a, b))
for i in range(12):

    PMT_ID = [a, b, i]
    layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(fname1, PMT_ID)
    layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(fname2, PMT_ID)
    layerID3, columnID3, cellID3, time3, x3, y3, z3, energy3 = get_data(fname3, PMT_ID)
    layerID4, columnID4, cellID4, time4, x4, y4, z4, energy4 = get_data(fname4, PMT_ID)
    layerID5, columnID5, cellID5, time5, x5, y5, z5, energy5 = get_data(fname5, PMT_ID)

    f = open('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/string{}{}/{}_TeV_string{}{}_output.txt'.format(neutrino_energy, a, b, neutrino_energy, a, b), 'a+')
    f.write("########################################################\n")
    f.write("\nPMT_ID: [{}, {}, {}]\n".format(*PMT_ID))
    f.write("\nBefore Removing Tails\n")
    f.write("Number of hits {} TeV Tau: {}\n".format(neutrino_energy, len(x1)))
    f.write("Number of hits {} TeV Had: {}\n".format(neutrino_energy, len(x2)))
    f.write("Number of hits {} TeV e-: {}\n".format(neutrino_energy, len(x3)))
    f.write("Total number of hits tau neutrino {} TeV: {}\n".format(neutrino_energy, len(x4)))
    f.write("Total number of hits electron neutrino {} TeV: {}\n".format(neutrino_energy, len(x5)))

    PMT_pos = PMT_ID_to_pos(PMT_ID)
    dx, dy, dz = seperation_vector('A', PMT_pos)
    d = distance_to_vertex(dx, dy, dz)
    theta = angle_to_vertex('A', dx, dy, dz)

    time1 = clean_data(time1, d)
    time2 = clean_data(time2, d)
    time3 = clean_data(time3, d)
    time4 = clean_data(time4, d)
    time5 = clean_data(time5, d)

    f.write("\nAfter Removing Tails\n")
    f.write("\nNumber of hits {} TeV Tau: {}\n".format(neutrino_energy, len(time1)))
    f.write("Number of hits {} TeV Had: {}\n".format(neutrino_energy, len(time2)))
    f.write("Number of hits {} TeV e-: {}\n".format(neutrino_energy, len(time3)))
    f.write("Total number of hits tau neutrino {} TeV: {}\n".format(2*neutrino_energy, len(time4)))
    f.write("Total number of hits electron neutrino {} TeV: {}\n".format(2*neutrino_energy, len(time5)))

    if len(time1) < 50:
        continue
    elif len(time1) < 200:
        bin_size = 20
    elif len(time1) < 400:
        bin_size = 10
    elif len(time1) < 700:
        bin_size = 4
    elif len(time1) < 1400:
        bin_size = 2
    else:
        bin_size = 1

    bins = []
    time = [time1, time2, time3, time4, time5]
    for j in range(5):
        bin = int((max(time[j]) - min(time[j])) / bin_size)
        bins.append(bin)


    if i <= 3:
        l = 0
    elif i <= 7:
        l = 1
    else:
        l = 2

    n1, bins1, patches1 = ax1[l][i % 4].hist(time1, bins=bins[0], range=(min(time1), max(time1)), histtype='step', color='C0', label='tau')
    n2, bins2, patches2 = ax1[l][i % 4].hist(time2, bins=bins[1], range=(min(time2), max(time2)), histtype='step', color='C1', label='had')
    n4, bins4, patches4 = ax1[l][i % 4].hist(time4, bins=bins[3], range=(min(time4), max(time4)), histtype='step', color='C2', label='Sum')
    ax1[l][i % 4].set_xlim(min(time4), max(time4))
    ax1[l][i % 4].legend(fontsize=11)
    ax1[l][i % 4].set_title('PMT_ID=[{},{},{}], distance={:.2f} m\nangle={:.2f} degrees, bin_size={} ns'
               .format(*PMT_ID, d, theta, bin_size), fontsize=11)

    n3, bins3, patches3 = ax2[l][i % 4].hist(time3, bins=bins[2], range=(min(time3), max(time3)), histtype='step', color='C0', label='e-')
    n2, bins2, patches2 = ax2[l][i % 4].hist(time2, bins=bins[1], range=(min(time2), max(time2)), histtype='step', color='C1', label='had')
    n5, bins5, patches5 = ax2[l][i % 4].hist(time5, bins=bins[4], range=(min(time5), max(time5)), histtype='step', color='C2', label='Sum')

    ax2[l][i % 4].set_xlim(min(time5), max(time5))
    ax2[l][i % 4].legend(fontsize=11)
    ax2[l][i % 4].set_title('PMT_ID=[{},{},{}], distance={:.2f} m\nangle={:.2f} degrees, bin_size={} ns'
               .format(*PMT_ID, d, theta, bin_size), fontsize=11)


    f.write("\nTau time range: {:.2f}-{:.2f} ns\n".format(min(time1), max(time1)))
    f.write("Had time range: {:.2f}-{:.2f} ns\n".format(min(time2), max(time2)))
    f.write("e- time range: {:.2f}-{:.2f} ns\n".format(min(time3), max(time3)))
    f.write("Total time range tau neutrino: {:.2f}-{:.2f} ns\n".format(min(time4), max(time4)))
    f.write("Total time range electron neutrino: {:.2f}-{:.2f} ns\n".format(min(time5), max(time5)))

    peak_height1, peak_pos1 = get_peak(n1, bins1)
    peak_height2, peak_pos2 = get_peak(n2, bins2)
    peak_height3, peak_pos3 = get_peak(n3, bins3)
    peak_height4, peak_pos4 = get_peak(n4, bins4)
    peak_height5, peak_pos5 = get_peak(n5, bins5)

    FWHM1 = get_FWHM(n1, bins1)
    FWHM2 = get_FWHM(n2, bins2)
    FWHM3 = get_FWHM(n3, bins3)
    FWHM4 = get_FWHM(n4, bins4)
    FWHM5 = get_FWHM(n5, bins5)

    f.write("\nTau; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height1, peak_pos1, FWHM1))
    f.write("Had; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height2, peak_pos2, FWHM2))
    f.write("e-; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height3, peak_pos3, FWHM3))
    f.write("Total tau neutrino; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height4, peak_pos4, FWHM4))
    f.write("Total electron neutrino; peak_height, peak_pos, FWHM: {} photons, {} ns, {} ns\n".format(peak_height5, peak_pos5, FWHM5))
    f.close()

#fig1.suptitle('Time Residuals {} TeV Tau Neutrino; String_{}{}'.format(2*neutrino_energy, a,b), fontsize=14)
#fig2.suptitle('Time Residuals {} TeV Electron Neutrino; String_{}{}'.format(2*neutrino_energy, a,b), fontsize=14)
fig1.subplots_adjust(top=0.94, bottom=0.06, left=0.06, right=0.93, hspace=0.3, wspace=0.2)
fig2.subplots_adjust(top=0.94, bottom=0.06, left=0.06, right=0.93, hspace=0.3, wspace=0.2)
fig1.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/string{}{}/{}_TeV_tau_string{}{}_output.png'.format(neutrino_energy, a, b, neutrino_energy, a, b))
fig2.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/string{}{}/{}_TeV_e_string{}{}_output.png'.format(neutrino_energy, a, b, neutrino_energy, a, b))
stop = timeit.default_timer()
print('runtime=', stop-start)
plt.show()
