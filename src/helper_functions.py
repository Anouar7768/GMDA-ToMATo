import numpy as np
import re

from tqdm import tqdm


def read_xyz(filename, n_rows, n_columns):
    """
    Get coordinates from filename and return a vectorset with all the
    coordinates, in XYZ format.
    Parameters
    ----------
    filename : string
        Filename to read
    Returns
    -------
    atoms : list
        List of atomic types
    V : array
        (N,3) where N is number of atoms
    """

    f = open(filename, 'r')
    V = list()

    # Use the number of rows to not read beyond the end of a file
    for lines_read, line in enumerate(f):

        if lines_read == n_rows:
            break

        numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', line)
        numbers = [float(number) for number in numbers]

        if len(numbers) == n_columns:
            V.append(np.array(numbers))
        else:
            exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read + 2))

    f.close()
    V = np.array(V)
    return V


def get_conformations(raw):
    """
    Transform the raw matrix into conformation matrix
    :param raw: raw matrix
    :return: conformations : list of conformations
    """
    conformations = []
    for k in range(int(len(raw) / 10)):
        conformations.append([raw[k * 10:(k + 1) * 10]])
    return conformations


def get_sample_conformation(conformations, step):
    """
    Extract a sample of conformations from all conformations
    :param conformations: the total list of conformations
    :param step: step at which we bypass conformations
    :return: conform_sample: the extracted sample of conformations
    """
    conform_sample = []
    for k in range(len(conformations)):
        if k % step == 0:
            conform_sample.append(conformations[k])
    return conform_sample


def get_dihedral_sample(dihedral, step):
    """
    Extract a sample from dihedral dataset, to show the final clusters
    :param dihedral: original dihedral dataset
    :param step: step at which we bypass conformations
    :return: dihedral_sample: the extracted sample from dihedral dataset
    """
    dihedral_sample = []
    for k in range(len(dihedral)):
        if k % step == 0:
            dihedral_sample.append(dihedral[k])
    dihedral_sample = np.array(dihedral_sample)
    return dihedral_sample


def compute_RMSD_matrix(conform_sample):
    """
    Compute the RMSD matrix for a sample of conformations
    :param conform_sample: the sample of conformations
    :return: RMSD_m_sample: the associated RMSD matrix
    """
    n = len(conform_sample)
    RMSD_m_sample = np.zeros((n, n))
    for k in tqdm(range(n)):
        for i in range(n):
            RMSD_m_sample[k, i] = quaternion_rmsd(np.array(conform_sample[k][0]), np.array(conform_sample[i][0]))
    return RMSD_m_sample


def quaternion_rotate(X, Y):
    """
    Calculate the rotation
    Parameters
    ----------
    X : array
        (N,D) matrix, where N is the number of points and D is dimension.
    Y: array
        (N,D) matrix, where N is the number of points and D is dimension.
    Returns
    -------
    rot : matrix
        Rotation matrix (D,D)
    """
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    Q = np.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = np.asarray([np.dot(Q[k].T, W[k]) for k in range(N)])
    W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    A = np.sum(Qt_dot_W, axis=0)
    eigen = np.linalg.eigh(A)
    r = eigen[1][:, eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot


def quaternion_transform(r):
    """
    Get optimal rotation
    note: translation will be zero when the centroids of each molecule are the
    same
    """
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3, :3]
    return rot


def makeW(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    W = np.asarray([
        [r4, r3, -r2, r1],
        [-r3, r4, r1, r2],
        [r2, -r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return W


def makeQ(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    Q = np.asarray([
        [r4, -r3, r2, r1],
        [r3, r4, -r1, r2],
        [-r2, r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return Q


def centroid(X):
    """
    Calculate the centroid from a vectorset X.
    https://en.wikipedia.org/wiki/Centroid
    Centroid is the mean position of all the points in all of the coordinate
    directions.
    C = sum(X)/len(X)
    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    C : float
        centeroid
    """
    C = X.mean(axis=0)
    return C


def quaternion_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    based on doi:10.1016/1049-9660(91)90036-O
    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    rmsd : float
    """
    rot = quaternion_rotate(P, Q)
    P = np.dot(P, rot)
    return rmsd(P, Q)


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    rmsd : float
        Root-mean-square deviation
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
    return np.sqrt(rmsd / N)
