import numpy as np
import matplotlib.pyplot as plt
import analysis_functions as af
from scipy.optimize import curve_fit

simID = 1
energy = 1
type = 2

#fname1 = '{} TeV/tau_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)
fname1 = '{} TeV/e_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)
fname2 = '{} TeV/had_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)
#fname = '{}TeV_{}x{}/muon{}_nt_Ntuple.csv'.format(energy, size, size, simID)

#fname3 = '{} TeV/tau_had_{}.csv'.format(energy, simID)
fname3 = '{} TeV/e_had_{}.csv'.format(energy, simID)


#print(len(data))

# select specific PMT by inputing numbers for layerID, columnID, and cellID
PMT_ID = [1, 1, 3]

#layerID, columnID, cellID, time, x, y, z, energy = af.get_data(fname, PMT_ID)
layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = af.get_data(fname1, PMT_ID)
layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = af.get_data(fname2, PMT_ID)
layerID3, columnID3, cellID3, time3, x3, y3, z3, energy3 = af.get_data(fname3, PMT_ID)

PMT_pos = af.PMT_ID_to_pos(PMT_ID, 150)
distance = af.distance_to_vertex(simID, PMT_pos)
angle = af.angle_to_vertex(simID, PMT_pos)

guesses1 = [137, 30, 500, 0.3]
guesses2 = [132, 25, 820, 0.1]
guesses3 = [130, 20, 1000, 0.3]
guesses4 = [213, 217, 15, 15, 700, 450, 0.25, 0.35]

tmin = 168
tmax = 175
bin_size = 0.25  # ns
dt = tmax - tmin
bins = int(dt/bin_size)

plt.figure(figsize=(8,5))

#n, bins, patches = plt.hist(time, bins=bins, range=(tmin, tmax), histtype='step', color='C1', label='muon')
n1, bins1, patches1 = plt.hist(time1, bins=bins, range=(tmin, tmax), histtype='step', color='C1', label='e-')
# n1, bins1, patches1 = plt.hist(time1, bins=bins, range=(tmin, tmax), histtype='step', color='C1', label='tau')
n2, bins2, patches2 = plt.hist(time2, bins=bins, range=(tmin, tmax), histtype='step', color='C3', label='had')
n3, bins3, patches3 = plt.hist(time3, bins=bins, range=(tmin, tmax), histtype='step', color='C0', label='Sum')

t_fit = np.linspace(tmin, tmax, bins)
t_plot = np.linspace(tmin, tmax, 4000)

# fitparams1, fitcov1 = curve_fit(af.bi_gaussian, t_fit, n1, p0=guesses1)
# fitparams2, fitcov2 = curve_fit(af.bi_gaussian, t_fit, n2, p0=guesses2)
# fitparams3, fitcov3 = curve_fit(af.bi_gaussian, t_fit, n3, p0=guesses3)
# fitparams4, fitcov4 = curve_fit(af.fit_function, t_fit, n3, p0=guesses4, maxfev=1000000)


# print(fitparams1)
# print(fitparams2)
# print(fitparams3)
# print(fitparams4)


# plt.plot(t_plot, af.bi_gaussian(t_plot, *fitparams1))
# plt.plot(t_plot, af.bi_gaussian(t_plot, *fitparams2))
# plt.plot(t_plot, af.bi_gaussian(t_plot, *fitparams4[[0,2,4,6]]))
# plt.plot(t_plot, af.bi_gaussian(t_plot, *fitparams4[[1,3,5,7]]))
# plt.plot(t_plot, af.fit_function(t_plot, *fitparams4))
# plt.plot(t_plot, af.bi_gaussian(t_plot, *fitparams1)+af.bi_gaussian(t_plot, *fitparams2))

# from scipy.stats import chisquare

# chisq, p = chisquare(n3, af.fit_function(t_fit, fitparams1[0], fitparams2[0], fitparams1[1], fitparams2[1], fitparams1[2], fitparams2[2], fitparams1[3], fitparams2[3]), ddof=8)
# chisq, p = chisquare(n3, af.fit_function(t_fit, *fitparams4), ddof=8)
# print(p)

# plt.title('Chiquare/dof: {:.4f}, PMT: {}, distance = {:.2f} m, angle = {:.2f} deg'.format(p, PMT_ID, distance, angle), fontsize=14)
plt.title('PMT: {}, distance = {:.2f} m, angle = {:.2f} deg'.format(PMT_ID, distance, angle), fontsize=14)
plt.xlabel("Time [ns]", fontsize=14)
plt.ylabel("Photon Count", fontsize=14)
plt.tick_params(axis="x", labelsize=12)
plt.tick_params(axis="y", labelsize=12)
plt.legend(fontsize=14)
plt.show()
