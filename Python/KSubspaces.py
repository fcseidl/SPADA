#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 16:10:20 2020

@author: fcseidl

k-subspaces clustering algorithm generalizes k-means to cluster points near 
hyperplanes rather than centroids.
"""

import numpy as np
from sklearn.decomposition import PCA
    

def assignment_step(X, bases, origins):
    """
    Assign each datapoint to the nearest subspace defined by the given bases
    (and, in the case of affine subspaces, origins).

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        Data matrix, assignment is performed on rows.
    bases : array, shape (k, d, n_features)
        The rows of bases[i] are an orthogonal basis for the ith subspace.
    origins : array, shape (k, n_features), optional
        origins[i] is a point in the ith subspace.

    Returns
    -------
    tot_res : float
        Sum of residuals of points to nearest subspaces.
    labels : array, shape (n_samples,)
        Clustering labels of each sample, in [0, k).

    """
    tot_res = 0
    labels = []
    k, d, n_feat = bases.shape
    
    for x in X: # for each datapoint
        lbl = 0
        residual = np.linalg.lstsq(bases[0].T, x - origins[0], rcond=None)[1]
        for i in range(1, k):   # for each subspace
            res = np.linalg.lstsq(bases[i].T, x - origins[i], rcond=None)[1]
            if res < residual:
                residual = res
                lbl = i
        labels.append(lbl)
        tot_res += residual
        
    return tot_res, np.array(labels)


def update_step(X, k, d, labels):
    """
    Update subspaces to best-fit clusters.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        Data matrix.
    k : int
        Number of clusters
    d : int
        Dimension of subspaces.
    labels : array, shape (n_samples,)
        Clustering labels of each sample, in [0, k).

    Returns
    -------
    bases : array, shape (k, d, n_features)
        The rows of bases[i] are an orthogonal basis for the ith subspace.
    origins : array, shape (k, n_features)
        origins[i] is a point in the ith subspace.

    """
    n_samp, n_feat = X.shape
    
    bases = []
    origins = []
    
    for i in range(k):
        idxs = np.argwhere(labels == i).flatten()  # indices of cluster k
        if len(set(idxs)) != 0:
            cluster = X[idxs]
            pca = PCA(n_components=d).fit(cluster)
            bases.append(pca.components_)
            origins.append(pca.mean_)
        else:
            bases.append(np.zeros(d, n_feat))
            origins.append(np.zeros(n_feat))
    
    return np.array(bases), np.array(origins)
            

def KSubspaces(X, k, d, max_iter=1000, n_restarts=10):
    """
    Perform k-subspaces clustering, a generalization of k-means. 
    Initialization subspaces are chosen to contain random groups of points.

    Parameters
    ----------
    X : array, shape (n_samples, n_features)
        Data matrix, clustering is performed on rows.
    k : int
        Number of clusters. Must be >= 2.
    d : int
        Dimension of subspaces. Must satisfy 0 <= (d + 1) * k < n_samples.
    max_iter : int, optional
        Maximum number of iterations per restart. The default is 1000.
    n_restarts : int, optional
        Number of times to restart the algorithm with new initial conditions.

    Returns
    -------
    labels : array, shape (n_samples,)
        Clustering labels of each sample, in [0, k).
    bases : array, shape (k, d, n_features)
        The rows of bases[i] are an orthogonal basis for the ith subspace.
    origins : array, shape (k, n_features)
        origins[i] is a point in the ith subspace.

    """
    n_samp, n_feat = X.shape
    
    # check for bad input
    assert(k >= 2)
    assert(0 <= d and (d + 1) * k < n_samp)
    
    res_best = sum([ np.linalg.norm(x) for x in X ])
    labels_best = bases_best = origins_best = None
    
    for _ in range(n_restarts):
        labels = np.random.randint(k, size=n_samp) # initialize
        labels_prev = labels
        for it in range(max_iter):
            bases, origins = update_step(X, k, d, labels)
            res, labels = assignment_step(X, bases, origins)
            if (labels == labels_prev).all():   # converged
                break
            labels_prev = labels
        if res < res_best:
            labels_best = labels
            bases_best = bases
            origins_best = origins
    
    return labels_best, bases_best, origins_best
