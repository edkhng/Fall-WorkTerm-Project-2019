import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import gaussian

# fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM1/bad/5000000_photons_violet_nt_Ntuple.csv'
# fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM1/data/test_with_pitchfork.csv'
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM1/data/photon_violet_merge_pitchfork.csv'
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM1/data/photon_uv2_pitchfork.csv'
fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM1/data/photon_blue0_pitchfork.csv'
# fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM3/photon_violet4_pitchfork.csv'
# fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM3/photon_blue2_pitchfork.csv'

data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')
data3 = np.loadtxt(fname3, delimiter=',', comments='#')
# layerID = data[:,0]
# columnID = data[:,1]
# cellID = data[:,2]
# time = data[:,3] - 241.489 # [ns]
# select = (layerID == 0) & (cellID == 2)
# print(len(time))
# time = time[select]

time1 = data1 - 241.489 # [ns]
time2 = data2 - 241.489 # [ns]
time3 = data3 - 241.489 # [ns]
# time1 = data1 - 391.915 # [ns]
# time2 = data2 - 391.915 # [ns]

print(len(time1))
print(len(time2))
print(len(time3))

def my_gaussian(x, mu, sigma):
    return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-(x-mu)**2/(2*sigma**2))

mu = 0
sigma = 10/(2*np.sqrt(2*np.log(2)))

n1, bins1, patches1 = plt.hist(time1, bins=300, range=(0,300), density=True)
n2, bins2, patches2 = plt.hist(time2, bins=300, range=(0,300), density=True)
n3, bins3, patches3 = plt.hist(time3, bins=300, range=(0,300), density=True)
plt.close()

xx1 = np.linspace(-149, 150, len(n1))
pulse1 = my_gaussian(xx1, 0, sigma)
tt1 = np.convolve(n1, pulse1)
tt1 / sum(tt1)

xx2 = np.linspace(-149, 150, len(n2))
pulse2 = my_gaussian(xx2, 0, sigma)
tt2 = np.convolve(n2, pulse2)
tt2 / sum(tt2)

xx3 = np.linspace(-149, 150, len(n3))
pulse3 = my_gaussian(xx3, 0, sigma)
tt3 = np.convolve(n3, pulse3)
tt3 / sum(tt3)
# fname = "c:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM1/data/['P2'],['SDOM1'],up,violet,['20V'],['2500Hz'].csv"
# data = np.loadtxt(fname, delimiter=',', comments='#')
# t = data[0,:]
# x = data[1,:]
# background = np.mean(x[0:50])
# print(background)
# x -= background
# x /= sum(x)

fname1 = "c:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM1/data/['P2'],['SDOM1'],up,uv,['20V'],['2500Hz'],corrected.csv"
# fname1 = "c:/Users/Edmond Ng/Documents/WorkTerm 2/Data/optical_properties_analysis_mine/P2_sDOM3/['P2'],['SDOM3'],up,violet,['20V'],['2500Hz'],corrected.csv"
data1 = np.loadtxt(fname1, delimiter=',', comments='#')
t1 = data1[0,:]
x1 = data1[1,:]
x1 /= sum(x1)

a = np.linspace(0, 299, 300)
aa = np.linspace(-150, 448, 599)

plt.figure()
# plt.step(a, n, where='post', label='Sim Raw Data')
# plt.plot(xx, pulse, label='Pulse Shape')
plt.step(aa, tt3, where='pre', lw=1, label='Simulation uv')
# plt.step(aa, tt1, where='pre', lw=1, label='Simulation blue')
# plt.step(aa, tt2, where='pre', lw=1, label='Simulation violet')
# plt.plot(aa, tt, label='Simulation')
plt.step(t1, x1, where='post', lw=1, label='STRAW Data corrected uv')
# plt.plot(t1, x1, lw=1, label='STRAW Data corrected blue')
# plt.hist(time3, bins=300, range=(0,300), histtype='step', density = True, label='uv')
# plt.hist(time2, bins=300, range=(0,300), histtype='step', density = True, label='violet')
# plt.hist(time1, bins=300, range=(0,300), histtype='step', density = True, label='blue')
# plt.step(t, x, where='post', lw=1, label='STRAW Data')
plt.title("Photon Time Distribution, P2_sDOM1")
plt.xlim(-100,200)
plt.ylim(0.00003, 0.1)
plt.yscale('log')
plt.xlabel('Time [ns]')
plt.ylabel('Normalized Photon Count')
plt.legend()
plt.show()

