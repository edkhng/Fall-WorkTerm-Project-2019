import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
import timeit

fname1 = '50 TeV Electron/e_50TeV_merge.csv'
fname2 = '50 TeV Had/had_50TeV_merge.csv'
start = timeit.default_timer()
layerID1, columnID1, cellID1, time1, x1, y1, z1, energy1 = get_data(fname1)
layerID2, columnID2, cellID2, time2, x2, y2, z2, energy2 = get_data(fname2)
stop = timeit.default_timer()
print('runtime=', stop-start)

PMT_ID = [1, 1, 1]

select1 = (layerID1 == PMT_ID[0]) & (columnID1 == PMT_ID[1]) & (cellID1 == PMT_ID[2])

layerID1 = layerID1[select1]
columnID1 = columnID1[select1]
cellID1 = cellID1[select1]
time1 = time1[select1]

select2 = (layerID2 == PMT_ID[0]) & (columnID2 == PMT_ID[1]) & (cellID2 == PMT_ID[2])

layerID2 = layerID2[select2]
columnID2 = columnID2[select2]
cellID2 = cellID2[select2]
time2 = time2[select2]

N = min(len(time1)/87, len(time2)/99)
time_range = get_range(time2)
bin_size = get_bin_size(time_range, N)

bins_a = int((max(time2) - min(time2))/bin_size)
range_a = (min(time2), max(time2))

n1/87, bins1, patches1 = plt.hist(time1, bins=bins_a, range=range_a, histtype='step', color='C0', label='e-')
n2/99, bins2, patches2 = plt.hist(time2, bins=bins_a, range=range_a, histtype='step', color='C1', label='had')
tmin_a, tmax_a = get_range_plot(n4, bins4)
plt.xlim(tmin_a, tmax_a)
plt.xlabel('time [ns]')
plt.ylabel('photon count')
plt.legend()
plt.savefig('test_e_had.png')
