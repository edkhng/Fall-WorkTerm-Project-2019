import numpy as np
import matplotlib.pyplot as plt
from analysis_functions import *
import timeit

fname1 = '50 TeV Electron/e_50TeV_merge.csv'
fname2 = '50 TeV Had/had_50TeV_merge.csv'

N = min(len(time1)/100, len(time2)/100)
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
