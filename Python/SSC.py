"""
@author: garg26

Driver module for sparse subspace clustering in Python.

@modified: fcseidl

Added documentation and a wrapper function which allows clustering to be 
performed in a single line of code.
"""

import numpy as np
from sklearn import cluster

from DataProjection import DataProjection
from BuildAdjacency import BuildAdjacency
from OutlierDetection import OutlierDetection
from BestMap import BestMap
from SpectralClustering import SpectralClustering
from SparseCoefRecovery import SparseCoefRecovery

from KSubspaces import KSubspaces


def sparseSubspaceClustering(
        X, n_clusters=8, ground_truth=None, remove_outliers=False, affine=False, OptM='Lasso', lam=1e-2):
    """
    Wrapper function which removes outliers then performs sparse subspace 
    clustering, added by Frank Seidl.

    Parameters
    ----------
    X : array, shape (N, D)
        Data matrix of N samples, D features.
    n_clusters : int, optional
        Number of clusters to form. The default is 8.
    ground_truth : array, shape (N,), optional
        True labels of each datapoint. Not required for clustering.
    remove_outliers : bool, optional
        Whether or not to perform outlier detection. Default is False.
    affine : bool, optional
        Whether to use affine (rather than linear) subspaces. 
        The default is False.
    OptM : string, optional
        Type of optimization, {'L1Perfect','L1Noise','Lasso','L1ED'}.
        The default is 'Lasso'.
    lam : float, optional
        Regularizartion parameter of LASSO, typically between 0.001 and
        0.1 or the noise level for 'L1Noise'. The default is 1e-2.

    Raises
    ------
    Exception
        If data has too many outliers.

    Returns
    -------
    labels : 1D array
        Cluster membership of each non-outlier sample sample.
    ground_truth : 1D array
        True label of each non-outlier. Nonsense if ground_truth parameter was 
        not provided.

    """
    X = X.T
    D, N = X.shape
    Fail = False
    # build dummy ground truth if needed
    if ground_truth is None:
        ground_truth = np.ones(N)
        for i in range(n_clusters):
            ground_truth[i] = i + 1
    # express data as sparse combinations of each other
    C = SparseCoefRecovery(X, cst=affine, Opt=OptM, lmbda=lam)
    # Make small values 0
    eps = np.finfo(float).eps
    C[np.abs(C) < eps] = 0
    # detect and remove outliers?
    if remove_outliers:
        C, ground_truth, OutlierIndx, Fail = OutlierDetection(C, ground_truth)
    if not Fail:
        adj = BuildAdjacency(C, K=0)  # K=0 to keep all edges
        # sklearn implementation of spectral clustering rather than JHU port
        labels = cluster.spectral_clustering(adj, n_clusters)
    else:
        raise Exception('Too many outliers detected for sparse subspace clustering')
    return labels, ground_truth


def SSC_test():
    """Basic test to check SSC."""
    
    # generating data

    D = 30  # Dimension of ambient space
    n = 2  # Number of subspaces
    d = 2
    N1 = 50
    N2 = 50  # N1 and N2: number of points in subspace 1 and 2
    # Generating N1 points in a d dim. subspace
    X1 = np.random.randn(D, d) * np.random.randn(d, N1)
    # Generating N2 points in a d dim. subspace
    X2 = np.random.randn(D, d) * np.random.randn(d, N2)
    X = np.concatenate((X1, X2), axis=1)

    # Generating the ground-truth for evaluating clustering results
    s = np.concatenate((1 * np.ones([1, N1]), 2 * np.ones([1, N2])), axis=1)
    r = 0  # Enter the projection dimension e.g. r = d*n, enter r = 0 to not project
    Cst = 0  # Enter 1 to use the additional affine constraint sum(c) == 1
    OptM = 'L1Noise'  # OptM can be {'L1Perfect','L1Noise','Lasso','L1ED'}
    lmbda = 0.001  # Regularization parameter in 'Lasso' or the noise level for 'L1Noise'
    # Number of top coefficients to build the similarity graph, enter K=0 for using the whole coefficients
    K = 0 #max(d1, d2)
    if Cst == 1:
        K = d + 1  # For affine subspaces, the number of coefficients to pick is dimension + 1

    Xp = DataProjection(X, r, 'NormalProj')
    
    # testing clustering
    
    '''
    # using wrapper function
    Grps, sc = sparseSubspaceClustering(
        Xp.T, n, ground_truth=s, affine=Cst, OptM=OptM, lam=lmbda)
    Grps = BestMap(sc, Grps)
    Missrate = float(np.sum(sc != Grps)) / sc.size
    print("\n\nMisclassification rate: {:.4f} %\n\n".format(Missrate * 100))
    '''
    
    '''
    # calling internal tools
    CMat = SparseCoefRecovery(Xp, Cst, OptM, lmbda)
    # Make small values 0
    eps = np.finfo(float).eps
    CMat[np.abs(CMat) < eps] = 0

    CMatC, sc, OutlierIndx, Fail = OutlierDetection(CMat, s)

    if Fail == False:
        CKSym = BuildAdjacency(CMatC, K)
        Grps = SpectralClustering(CKSym, n)
        Grps = BestMap(sc, Grps)
        Missrate = float(np.sum(sc != Grps)) / sc.size
        print("\n\nMisclassification rate: {:.4f} %\n\n".format(Missrate * 100))
    else:
        print("Something failed")
    '''
    
    labels, _, __ = KSubspaces(Xp.T, n, d)
    labels = BestMap(s, labels)
    Missrate = float(np.sum(s != labels)) / s.size
    print("\n\nMisclassification rate: {:.4f} %\n\n".format(Missrate * 100))

if __name__ == "__main__":
    SSC_test()
