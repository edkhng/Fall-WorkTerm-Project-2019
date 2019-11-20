# import numpy as np
# import os

fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/100 TeV/tau_100TeV_b.out'

decay_times = []
event_IDs = []
hits = []

f = open(fname, 'r')
for line in f.readlines():
    words = line.split()
    for word in words:
        if word == ';decay_time:':
            decay_time = words[words.index(word) + 1]
            decay_times.append(float(decay_time))

        if word == 'ID:':
            event_ID = words[words.index(word) + 1]
            event_IDs.append(float(event_ID))

        if word == 'Hit:':
            hit = words[words.index(word) + 1]
            hits.append(float(hit))

f = open('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/100 TeV/100_TeV_tau_decay_times.csv', 'a+')
f.write('#Event ID,Hits,Decay Time\n')
for i in range(len(event_IDs)):
    f.write('{},{},{}\n'.format(event_IDs[i],hits[i],decay_times[i]))
f.close()

# import matplotlib.pyplot as plt
# print(max(decay_times))
# print(len(decay_times))
# print(decay_times)
# plt.hist(decay_times, bins=15, range=(0,40))
# plt.show()
# print(len(decay_times))
# print(decay_times)
# print(len(event_IDs))
# # print(event_IDs)
# print(len(hits))
# print(hits)
