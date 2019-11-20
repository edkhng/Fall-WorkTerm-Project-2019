import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/50 TeV/50_TeV_tau_decay_times.csv'
data = np.loadtxt(fname, delimiter=',')
# data = sorted(data, key=lambda x: x[0])
# np.savetxt(fname, data, fmt=['%.0f', '%.0f', '%.4f'], delimiter=',')

event_ID = data[:,0]
hits = data[:,1]
decay_times = data[:,2]


# plt.figure(figsize=(10,7))
# plt.hist(hits, bins=40)
# plt.show()


def decay(x, mu):
    return mu*np.exp(-mu*x)

n, bins = np.histogram(decay_times, bins=10, range=(0,40), density=True)
sigma = 1/np.sqrt(n)
fitparams, fitcov = curve_fit(decay, bins[:-1], n, sigma=sigma, p0=[0.1])

print(fitparams)
print(np.sqrt(fitcov))

plt.figure(figsize=(10,7))
plt.plot(bins, decay(bins, fitparams), label='fit')
plt.hist(decay_times, bins=10, range=(0,40), histtype='step', density=True, align='mid', label='simulation')
plt.xlabel('time [ns]')
plt.xlim(0,40)
plt.legend()
plt.title('mean decay time: {:.2f} +/- {:.2f} [ns]'.format(float(1/fitparams), float(1/fitparams*np.sqrt(fitcov)/fitparams)))
plt.savefig('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/100 TeV/decay_time_bin=6ns.png')
plt.show()
