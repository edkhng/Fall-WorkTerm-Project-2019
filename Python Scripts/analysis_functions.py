import numpy as np
import matplotlib.pyplot as plt

def get_data(fname, PMT_ID = ''):
    """Obtains the data from the csv file. If PMT_ID is specified then
       the data for that specific PMT is returned if then all the data
       returned."""
    data = np.loadtxt(fname, delimiter=',', comments='#')
    layerID = data[:,0]
    columnID = data[:,1]
    cellID = data[:,2]
    time = data[:,3]  # [ns]
    x = data[:,4]/1000  # [m]
    y = data[:,5]/1000  # [m]
    z = data[:,6]/1000  # [m]
    energy = data[:,7]
    if PMT_ID:
        select = (layerID == PMT_ID[0]) & (columnID == PMT_ID[1]) & (cellID == PMT_ID[2])

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
    return abs(amp)*np.exp(-4*np.log(2)*((x-pos)/(wid))**2)


def get_peak(n, bins):
    """Find the approximate location and height of the peak of the distribution."""
    peak_height = max(n)
    if len(bins[np.where(n==max(n))]) != 0:
        peak_pos = np.mean(bins[np.where(n==max(n))])
    else:
        peak_pos = float(bins[np.where(n==max(n))])
    return peak_height, peak_pos


def get_range(time):
    """Get time range of first 85% of data"""
    time.sort()
    N = len(time)
    time = time[0:-int(N*0.2)]
    time_range = max(time) - min(time)
    return time_range


def get_range_plot(n, bins):
# def get_range_plot(time):
    """Return a tmin and tmax that covers an appropriate amount of the data for plotting."""
    peak_height = get_peak(n, bins)[0]
    cutoff = max(int(peak_height/10), 4)
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


def get_range_fit(n, bins, peak_height):
    """Return a tmin and tmax that covers an appropriate amount of the data for fitting."""
    third_max = int(peak_height/3)
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


def seperation_vector(vertex_type, PMT_pos):
    """Returns the seperation vector between the vertex and the PMT."""
    if vertex_type == 'A':
        vertex = [-10, -20, -62.5]
    elif vertex_type == 'B':
        vertex = [-20, -62.5, -10]
    dx = vertex[0] - PMT_pos[0]
    dy = vertex[1] - PMT_pos[1]
    dz = vertex[2] - PMT_pos[2]
    return dx, dy, dz


def distance_to_vertex(dx, dy, dz):
    """Returns the seperation distance between the vertex and the PMT."""
    return np.sqrt(dx**2 + dy**2 + dz**2)


def angle_to_vertex(vertex_type, dx, dy, dz):
    """Returns the angle between the vertex and the PMT."""

    if vertex_type == 'A':
        angle = np.arctan((np.sqrt(dx**2 + dy**2)) / dz)
    elif vertex_type == 'B':
        angle = np.arctan((np.sqrt(dx**2 + dz**2)) / dy)

    # convert angle to degrees
    angle = angle*180/np.pi
    return angle
