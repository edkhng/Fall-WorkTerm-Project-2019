"""
All the useful reoccuring functions in one place.
"""

import numpy as np
import matplotlib.pyplot as plt

def get_data(data, PMT_ID = ''):
    """Obtains the data from the csv file. If PMT_ID is specified then
       the data for that specific PMT is returned if then all the data
       returned."""
    layerID = data[:,0]
    columnID = data[:,1]
    cellID = data[:,2]
    time = data[:,3]  # [ns]
    x = data[:,4]/1000  # [m]
    y = data[:,5]/1000  # [m]
    z = data[:,6]/1000  # [m]
    energy = data[:,7]
    if PMT_ID:
        select = (layerID == PMT_ID[0]) & (columnID == PMT_ID[1]) & (cellID == PMT_ID[2]) & (time < 900)

        layerID = layerID[select]
        columnID = columnID[select]
        cellID = cellID[select]
        time = time[select]
        energy = energy[select]
        x = x[select]
        y = y[select]
        z = z[select]
    else:
        return layerID, columnID, cellID, time, x, y, z, energy

    return layerID, columnID, cellID, time, x, y, z, energy


def find_good_fits(data1, data2):
    """Keep ones with amplitude ratio less than 4, a FWHM ratio
     less than 4, and a time difference less than 100 ns. Also
     remove fits with inf chi^2."""
    time_diff = abs(data1[:,7] - data1[:,6])
    amp_ratio = data1[:,10]/data1[:,11]
    wid_ratio = data1[:,8]/data1[:,9]
    chisq1 = data1[:,22]
    chisq2 = data2[:,14]
    clean_data1 = []
    clean_data2 = []
    select1 = time_diff < 100
    select2 = (amp_ratio < 4) & (amp_ratio > 1/4)
    select3 = (wid_ratio < 4) & (wid_ratio > 1/4)
    select4 = ~np.isinf(chisq1)
    select5 = ~np.isinf(chisq2)
    for i in range(26):
        clean_data1.append(data1[:,i][select1 & select2 & select3 & select4 & select5])

    for i in range(18):
        clean_data2.append(data2[:,i][select1 & select2 & select3 & select4 & select5])
    return clean_data1, clean_data2



def get_fit_data(data, hits_range=(0,400000), decay_time_range=(0,70)):
    """Keep only the data within the time range and hits range"""
    if len(data) == 26:
        event_ID1, decay_time1, hits1, layerID1, columnID1, cellID1, pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b, upos1a, upos1b, uwid1a, uwid1b, uamp1a, uamp1b, ur1a, ur1b, chisq1, bins1, params1, p1 = data

        select1 = (decay_time1 < decay_time_range[1]) & (decay_time1 > decay_time_range[0])

        event_ID1, decay_time1, hits1, layerID1, columnID1, cellID1, pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b, upos1a, upos1b, uwid1a, uwid1b, uamp1a, uamp1b, ur1a, ur1b, chisq1, bins1, params1, p1 = event_ID1[select1], decay_time1[select1], hits1[select1], layerID1[select1], columnID1[select1], cellID1[select1], pos1a[select1], pos1b[select1], wid1a[select1], wid1b[select1], amp1a[select1], amp1b[select1], r1a[select1], r1b[select1], upos1a[select1], upos1b[select1], uwid1a[select1], uwid1b[select1], uamp1a[select1], uamp1b[select1], ur1a[select1], ur1b[select1], chisq1[select1], bins1[select1], params1[select1], p1[select1]

        select1 = (hits1 < hits_range[1]) & (hits1 > hits_range[0])

        event_ID1, decay_time1, hits1, layerID1, columnID1, cellID1, pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b, upos1a, upos1b, uwid1a, uwid1b, uamp1a, uamp1b, ur1a, ur1b, chisq1, bins1, params1, p1 = event_ID1[select1], decay_time1[select1], hits1[select1], layerID1[select1], columnID1[select1], cellID1[select1], pos1a[select1], pos1b[select1], wid1a[select1], wid1b[select1], amp1a[select1], amp1b[select1], r1a[select1], r1b[select1], upos1a[select1], upos1b[select1], uwid1a[select1], uwid1b[select1], uamp1a[select1], uamp1b[select1], ur1a[select1], ur1b[select1], chisq1[select1], bins1[select1], params1[select1], p1[select1]

        return event_ID1, decay_time1, hits1, layerID1, columnID1, cellID1, pos1a, pos1b, wid1a, wid1b, amp1a, amp1b, r1a, r1b, upos1a, upos1b, uwid1a, uwid1b, uamp1a, uamp1b, ur1a, ur1b, chisq1, bins1, params1, p1

    else:
        event_ID4, decay_time4, hits4, layerID4, columnID4, cellID4, pos4, wid4, amp4, r4, upos4, uwid4, uamp4, ur4, chisq4, bins4, params4, p4 = data

        select1 = (decay_time4 < decay_time_range[1]) & (decay_time4 > decay_time_range[0])

        event_ID4, decay_time4, hits4, layerID4, columnID4, cellID4, pos4, wid4, amp4, r4, upos4, uwid4, uamp4, ur4, chisq4, bins4, params4, p4 = event_ID4[select1], decay_time4[select1], hits4[select1], layerID4[select1], columnID4[select1], cellID4[select1], pos4[select1], wid4[select1], amp4[select1], r4[select1], upos4[select1], uwid4[select1], uamp4[select1], ur4[select1], chisq4[select1], bins4[select1], params4[select1], p4[select1]

        select1 = (hits4 < hits_range[1]) & (hits4 > hits_range[0])

        event_ID4, decay_time4, hits4, layerID4, columnID4, cellID4, pos4, wid4, amp4, r4, upos4, uwid4, uamp4, ur4, chisq4, bins4, params4, p4 = event_ID4[select1], decay_time4[select1], hits4[select1], layerID4[select1], columnID4[select1], cellID4[select1], pos4[select1], wid4[select1], amp4[select1], r4[select1], upos4[select1], uwid4[select1], uamp4[select1], ur4[select1], chisq4[select1], bins4[select1], params4[select1], p4[select1]

        return event_ID4, decay_time4, hits4, layerID4, columnID4, cellID4, pos4, wid4, amp4, r4, upos4, uwid4, uamp4, ur4, chisq4, bins4, params4, p4

def fit_function(x, pos1, pos2, wid1, wid2, amp1, amp2, r1, r2):
    """Fit function for two merged bifurcated gaussians."""
    bi_gaussian1 = bi_gaussian(x, pos1, wid1, amp1, r1)
    bi_gaussian2 = bi_gaussian(x, pos2, wid2, amp2, r2)
    return bi_gaussian1 + bi_gaussian2


def fit_function_gaussian(x, pos1, pos2, wid1, wid2, amp1, amp2):
    """Fit function for two merged gaussians."""
    gaussian1 = gaussian(x, pos1, wid1, amp1)
    gaussian2 = gaussian(x, pos2, wid2, amp2)
    return gaussian1 + gaussian2


def bi_gaussian(x, pos, wid, amp, r):
    """Function for bifurcated Gaussian centered at pos, where wid is the
       FWHM and r is the ratio of widths of the right side to left side."""

    mask = x < pos
    if r != 0:
        # seperate gaussian with different widths to the right and left of the center
        y1 = gaussian(x[mask], pos, r*wid/(r+1), amp)
        y2 = gaussian(x[~mask], pos, wid/(r+1), amp)
        return np.append(y1, y2)

    else:
        # normal gaussian
        return gaussian(x, pos, wid, amp)


def gaussian(x, pos, wid, amp):
    """Function for Gaussian centered at pos, where wid is the FWHM
       and amp is the amplitude."""
    return amp*np.exp(-4*np.log(2)*((x-pos)/wid)**2)


def get_peak(n, bins):
    """Find the approximate location and height of the peak of the distribution."""
    peak_height = max(n)
    if len(bins[np.where(n==max(n))]) != 0:
        peak_pos = np.mean(bins[np.where(n==max(n))])
    else:
        peak_pos = float(bins[np.where(n==max(n))])
    return peak_height, peak_pos


def get_range(time):
    """Get time range of first 80% of data"""
    time.sort()
    N = len(time)
    time = time[0:-int(N*0.2)]
    time_range = max(time) - min(time)
    return time_range


def get_range_plot(n, bins):
    """Return a tmin and tmax that covers an appropriate amount of the data for plotting."""
    peak_height = get_peak(n, bins)[0]
    cutoff = max(int(peak_height/10), 2)
    in_range = n >= cutoff
    tmin = bins[0]
    bins = bins[1:]
    bins = bins[in_range]
    tmax = bins[-1]
    return tmin, tmax


def get_FWHM(n, bins):
    """Find the FWHM of the distribution."""
    peak_height = max(n)
    half_max = int(peak_height/2)
    in_range = n >= half_max
    bins = bins[:-1]
    bins = bins[in_range]
    FWHM = max(bins) - min(bins)
    return FWHM


def get_range_fit(n, bins, peak_height, cutoff=3):
    """Return a tmin and tmax that covers an appropriate amount of the data for fitting."""
    third_max = int(peak_height/cutoff)
    in_range = n >= third_max
    bins = bins[:-1]
    bins = bins[in_range]
    tmin = min(bins)
    tmax = max(bins)
    return tmin, tmax

def get_bin_size(time_range, N):

    if N/time_range < 1:
        if time_range < 56:
            bin_size = 4
        else:
            bin_size = 8
    elif N/time_range < 3:
        if time_range < 40:
            bin_size = 2
        else:
            bin_size = 4
    elif N/time_range < 6:
        if time_range < 80:
            bin_size = 2
        else:
            bin_size = 4
    elif N/time_range < 10:
        if time_range < 60:
            bin_size = 1
        else:
            bin_size = 2
    elif N/time_range < 20:
        if time_range < 80:
            bin_size = 1
        else:
            bin_size = 2
    else:
        bin_size = 1

    return bin_size


def range_size(time, tmin, tmax):
    Ntot = len(time)
    select = (time >= tmin) & (time < tmax)
    time = time[select]
    N = len(time)
    return N/Ntot


def PMT_ID_to_pos(PMT_ID):
    """Function to obtain the cartesian coordinates of the PMT given the PMT_ID,
       the size of the detector and the type (i.e., number of layers, columns and cells)."""
    x = 50*PMT_ID[0] - 50
    y = 50*PMT_ID[1] - 50
    z = 12.5*PMT_ID[2] - 68.75
    return x, y, z


def seperation_vector(PMT_pos):
    """Returns the seperation vector between the vertex and the PMT."""
    vertex = [-10, -20, -62.5]
    dx = PMT_pos[0] - vertex[0]
    dy = PMT_pos[1] - vertex[1]
    dz = PMT_pos[2] - vertex[2]
    return dx, dy, dz


def distance_to_vertex(dx, dy, dz):
    """Returns the seperation distance between the vertex and the PMT."""
    return np.sqrt(dx**2 + dy**2 + dz**2)


def angle_to_vertex(dx, dy, dz):
    """Returns the angle between the vertex and the PMT."""
    angle = np.arctan((np.sqrt(dx**2 + dy**2)) / dz)

    # convert angle to degrees
    angle = angle*180/np.pi
    return angle
