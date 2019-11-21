"""
The output file for the tau simulation includes the decay time.
This script scans through the output file and records the decay
time for each event and the number of hits, then saves it to a csv.
"""
energy = 10
fname = 'C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/Outputs/tau_{}TeV.out'.format(energy, energy)

decay_times = []
event_IDs = []
hits = []

f = open(fname, 'r')
for line in f.readlines():
    # scan through every line in the output file and break it down into each word
    words = line.split()
    for word in words:
        if word == ';decay_time:':
            # scan each word for "decay_time", when found save the next word
            # which will be the decay time
            decay_time = words[words.index(word) + 1]
            decay_times.append(float(decay_time))

        if word == 'ID:':
            # scan each word for "ID", when found save the next word
            # which will be the event ID number
            event_ID = words[words.index(word) + 1]
            event_IDs.append(float(event_ID))

        if word == 'Hit:':
            # scan each word for "Hit", when found save the next word
            # which will be the number of hits for that event
            hit = words[words.index(word) + 1]
            hits.append(float(hit))

# save the results to a csv file
f = open('C:/Users/Edmond Ng/Documents/WorkTerm 2/Data/{} TeV/{}_TeV_tau_decay_times.csv', 'a+').format(energy, energy)
f.write('#Event ID,Hits,Decay Time\n')
for i in range(len(event_IDs)):
    # lists are out of order, so order by the ID number
    index = event_IDs.index(i)
    f.write('{},{},{}\n'.format(event_IDs[index],hits[index],decay_times[index]))
f.close()
