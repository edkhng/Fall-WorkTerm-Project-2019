"""
Fit the distribution of decay times to an exponential and return the mean.
The distribution and fit is then plotted. The bins and range will need
to be adjusted by the user.
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

energy = 100
fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/{}_TeV_tau_decay_times.csv'.format(energy, energy)
data = np.loadtxt(fname, delimiter=',')

event_ID = data[:,0]
hits = data[:,1]
decay_times = data[:,2]

def decay(x, mu):
    return mu*np.exp(-mu*x)

# adjust number of bins and range for plot
bin = 15
range = (0,55)

n, bins = np.histogram(decay_times, bins=bin, range=range, density=True)
sigma = np.sqrt(n)
fitparams, fitcov = curve_fit(decay, bins[:-1], n, sigma=sigma, p0=[0.1])

err = np.sqrt(np.diag(fitcov))
mean = 1/fitparams[0]
sigma_mean = 1/fitparams[0]*err[0]/fitparams[0]
print("Fit mean decay time: {} +/- {}".format(mean, sigma_mean))

# calculation for the mean decay time
E = energy*1e6  # energy in MeV
m0 = 1776.86  # rest mass of tau in MeV/c^2
gamma = np.sqrt(E**2 - m0**2)/m0  # not exact but pretty close
t = 2.9e-4*gamma  # 2.9e-4 is the mean decay time, result in ns
print("Expected mean decay time: {}".format(t))

plt.figure(figsize=(10,7))
plt.plot(bins, decay(bins, fitparams), label='fit')
plt.hist(decay_times, bins=bin, range=range, histtype='step', density=True, align='mid', label='simulation')
plt.xlabel('time [ns]')
plt.xlim(range)
plt.legend()
plt.title('mean decay time: {:.2f} +/- {:.2f} [ns]'.format(mean, sigma_mean))
# plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/decay_time_bin={:.2f}ns.png'.format(energy, energy, range[1]/bin))
plt.show()
