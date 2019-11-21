'''
Merge the tau and had data together and the e- and had data
'''
import numpy as np

energy = 10

fname1 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_{}TeV423_nt_Ntuple.csv'.format(energy, energy)
fname3 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/had_{}TeV234_nt_Ntuple.csv'.format(energy, energy)
fname2 = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_{}TeV163_nt_Ntuple.csv'.format(energy, energy)


data1 = np.loadtxt(fname1, delimiter=',', comments='#')
data2 = np.loadtxt(fname2, delimiter=',', comments='#')
data3 = np.loadtxt(fname3, delimiter=',', comments='#')
# data4 = np.loadtxt(fname1, delimiter=',', comments='#')

# data5 = np.concatenate((data1, data2, data3), axis=0)

data13 = np.append(data1, data3, axis=0)
data23 = np.append(data2, data3, axis=0)

np.savetxt('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/tau_had_merge_{}_TeV.csv'.format(energy, energy), data13, delimiter=',', comments='#')
np.savetxt('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/e_had_merge_{}_TeV.csv'.format(energy, energy), data23, delimiter=',', comments='#')

# np.savetxt('{}TeV_{}x{}/tau{}_nt_Ntuple_merge.csv'.format(energy, size, size, simID), data5, delimiter=',', comments='#')
