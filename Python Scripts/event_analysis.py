"""
This script analyzes the csv file from read_output.py. It outputs
the number of events in the specified hits range and decay time range,
and also plots out a breakdown of the number of hits by decay modes.
"""
import csv
energy = 50

hits = []
names = []
ID = []
decay_time = []
fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/{}_TeV_tau_decay_times_with_particles.csv'.format(energy, energy)
with open(fname,'r') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for i,row in enumerate(csv_reader):
        if i == 0:
            continue
        hits.append(float(row[1]))
        ID.append(float(row[0]))
        decay_time.append(float(row[2]))
        names.append(row[3])

hits = np.asarray(hits)
decay_time = np.asarray(decay_time)
ID = np.asarray(ID)

decay_time_range = (0, 7)
hits_range = (40000, 200000)
f = open('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Plots/key.txt'.format(energy), 'a')
f.write('\n############################################################################')
f.write("\nNumber of events with decay time: {} - {} ns\n".format(*decay_time_range))
f.write(str(len(ID[(decay_time > decay_time_range[0]) & (decay_time < decay_time_range[1])])))
f.write("\n\nNumber of events with hits: {} - {}\n".format(*hits_range))
f.write(str(len(ID[(hits > hits_range[0]) & (hits < hits_range[1])])))

select = (decay_time > decay_time_range[0]) & (decay_time < decay_time_range[1])
hits = hits[select]
ID = ID[select]

select = (hits > hits_range[0]) & (hits < hits_range[1])
hits = hits[select]
ID = ID[select]

f.write("\n\nNumber of hits with decay time {} - {} ns and hits {} - {}\n".format(*decay_time_range, *hits_range))
f.write(str(len(ID)))
f.write("\n\nEvent ID:\n")
f.write(str(ID))
f.close()

# Counting the decay modes
e = pi = mu = 0
for particle in names:
    if particle == 'e-':
        e += 1
    elif particle == 'mu-':
        mu += 1
    else:
        pi += 1

print("e- avg: {}".format(e))
print("pi avg: {}".format(pi))
print("muon avg: {}".format(mu))

pi = []
e = []
mu = []

for i in range(len(hits) - 1):
    if names[i+1] == 'e-':
        e.append(float(ID[i+1]))
    elif names[i+1] == 'mu-':
        mu.append(float(ID[i+1]))
    else:
        pi.append(float(ID[i+1]))

print(mu)
print(e)
print(pi)
import matplotlib.pyplot as plt

plt.figure()
plt.hist(e, bins=10, range=(0,358000), histtype='step', label='e-')
plt.hist(pi, bins=10, range=(0,358000), histtype='step', label='pions')
plt.hist(mu, bins=10, range=(0,358000), histtype='step', label='mu-')
plt.hist(mu+pi+e, bins=10, range=(0,358000), histtype='step', label='sum')
plt.title('Hits distribution')
plt.xlabel('hits')
plt.ylabel('frequency')
plt.legend()
plt.show()
