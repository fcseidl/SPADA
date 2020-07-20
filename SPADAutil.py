#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 23:28:59 2020

@author: fcseidl

Various auxiliary methods for SPADA research.
"""

import numpy as np
from scipy.optimize import linprog, minimize
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def bestSilhouetteKMeans(X, max_n_clusters=10):
    """
    Perform k-means clustering using average silhouette method to determine 
    optimal number of clusters.

    Parameters
    ----------
    X : array, (n_samples, n_features)
        Data matrix.
    max_n_clusters : int, optional
        Maximum number of clusters to try.

    Returns
    -------
    n_clusters : int
        Number of clusters used.
    cluster_centers : array, shape (n_clusters, n_features)
        Coordinates of cluster centers.
    labels : array, shape (n_samples,)
        Labels of each point.
    
    """
    # TODO: this doesn't work so well
    
    # determine centers, labels, and silhouette scores for each k.
    candidate_centers = []
    candidate_labels = []
    scores = []
    for k in range(2, max_n_clusters):
        kmeans = KMeans(n_clusters=k).fit(X)
        candidate_centers.append(kmeans.cluster_centers_)
        candidate_labels.append(kmeans.labels_)
        scores.append(silhouette_score(X, kmeans.labels_))
    
    # determine k with highest silhouette score.
    k_best = 2
    score_best = scores[0]
    for k in range(3, max_n_clusters):
        if scores[k - 2] > score_best:
            k_best = k
            score_best = scores[k - 2]
    
    return k_best, candidate_centers[k_best - 2], candidate_labels[k_best - 2]


def dropoutRate(Yn):
    """
    Compute fraction of dropouts for a gene.
    """
    return 1 - (np.count_nonzero(Yn) / Yn.shape[0])


def uniformFromUnitSphere(n, k=None):
    """
    Sample a uniform distribution over the n-dimensional unit sphere.

    Parameters
    ----------
    n : int
        number of dimensions
    k : int
        desired number of samples
        
    Returns
    -------
    k samples from a uniform distribution over the n-dimensional unit 
    sphere, with shape (n, k)
    
    """
    if k is None:
        x = np.random.randn(n)
        x /= np.linalg.norm(x)
        return x
    return np.array([ uniformFromUnitSphere(n) for _ in range(k) ])


def scaledSolidAngle(Y):
    """
    Approximate solid angle of cone(Y) for a matrix Y >= 0. Result is divided 
    by the solid angle of the positive orthant in R^N, where Y has shape (N, L)
    
    """
    N, L = Y.shape
    n_points = int(1e4)  # TODO: magic number
    count = 0
    for _ in range(n_points):
        x = uniformFromUnitSphere(N)
        c = np.zeros(L)
        lp = linprog(c, A_eq=Y, b_eq=x)
        if lp.success: count += 1
    return count / n_points


def inferBackground(X):
    """
    Guess a background profile which underlies each sample in normalized data 
    matrix.

    Parameters
    ----------
    X : array (n_samples, n_features)
        Data matrix.

    Returns
    -------
    b : array (n_features,)
        Background vector b which averages the L1 normalized columns of X.
    norm : float
        Frobenius norm of normalized X minus [b, ..., b].
        
    """
    # preprocess X
    M = X.shape[1]
    X = normalize(X, norm='l1', axis=0)
    b = np.mean(X, axis=1)
    norm = np.sqrt(
                sum(
                    [ np.linalg.norm(X[:, m] - b) ** 2
                     for m in range(M) ]
                    )
                )
    return b, norm
    
    
    
