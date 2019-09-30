import numpy as np

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
    return amp*np.exp(-4*np.log(2)*((x-pos)/(wid))**2)


def PMT_ID_to_pos(PMT_ID, size, type=0):
    """Function to obtain the cartesian coordinates of the PMT given the PMT_ID,
       the size of the detector and the type (i.e., number of layers, columns and cells)."""
    if size == 150:
        x = 50*PMT_ID[0] - 50
        y = 50*PMT_ID[1] - 50
        z = 12.5*PMT_ID[2] - 68.75

    elif size == 300:
        # type 1 is configuration with 3 layers, 3 columns and 12 cells
        if type == 1:
            x = 100*PMT_ID[0] - 100
            y = 100*PMT_ID[1] - 100
            z = 25*PMT_ID[2] - 137.5
        # type 2 is configuration with 6 layers, 6 columns and 24 cells
        elif type ==2:
            x = 50*PMT_ID[0] - 125
            y = 50*PMT_ID[1] - 125
            z = 12.5*PMT_ID[2] - 143.75
    return x, y, z


def seperation_vector(vertex_type, PMT_pos):
    """Returns the seperation vector between the vertex and the PMT."""
    if vertex_type == 1:
        vertex = [-20, -10, -62.5]
    elif vertex_type == 2:
        vertex = [-20, -62.5, -10]
    dx = vertex[0] - PMT_pos[0]
    dy = vertex[1] - PMT_pos[1]
    dz = vertex[2] - PMT_pos[2]
    return dx, dy, dz


def distance_to_vertex(vertex_type, PMT_pos):
    """Returns the seperation distance between the vertex and the PMT."""
    dx, dy, dz = seperation_vector(vertex_type, PMT_pos)
    return np.sqrt(dx**2 + dy**2 + dz**2)


def angle_to_vertex(vertex_type, PMT_pos):
    """Returns the angle between the vertex and the PMT."""
    dx, dy, dz = seperation_vector(vertex_type, PMT_pos)

    if vertex_type == 1:
        angle = np.arctan((np.sqrt(dx**2 + dy**2)) / dz)
    elif vertex_type == 2:
        angle = np.arctan((np.sqrt(dx**2 + dz**2)) / dy)

    # convert angle to degrees
    angle = angle*180/np.pi
    return angle
