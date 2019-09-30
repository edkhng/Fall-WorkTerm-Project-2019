'''
Merge the tau and had data together and the e- and had data
'''
import numpy as np

simID = 1
energy = 1
# size = 150
fname1 = '{} TeV/tau_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)
fname2 = '{} TeV/e_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)
fname3 = '{} TeV/had_{}TeV{}_nt_Ntuple.csv'.format(energy, energy, simID)


data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')
data3 = np.loadtxt(fname3, delimiter=',', comments='#')
# data4 = np.loadtxt(fname1, delimiter=',', comments='#')

# data5 = np.concatenate((data1, data2, data3), axis=0)

data13 = np.append(data1, data3, axis=0)
data23 = np.append(data2, data3, axis=0)

np.savetxt('{} TeV/tau_had_{}.csv'.format(energy, simID), data13, delimiter=',', comments='#')
np.savetxt('{} TeV/e_had_{}.csv'.format(energy, simID), data23, delimiter=',', comments='#')

# np.savetxt('{}TeV_{}x{}/tau{}_nt_Ntuple_merge.csv'.format(energy, size, size, simID), data5, delimiter=',', comments='#')
