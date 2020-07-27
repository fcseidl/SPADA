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
    centers = labels = []
    k_best = 1
    score_best = -1
    for k in range(2, max_n_clusters):
        kmeans = KMeans(n_clusters=k).fit(X)
        score = silhouette_score(X, kmeans.labels_)
        if score > score_best:
            k_best = k
            score_best = score
            centers = kmeans.cluster_centers_
            labels = kmeans.labels_
    return k_best, centers, labels


def dropoutRate(Yn):
    """
    Compute fraction of dropouts for a gene.
    """
    return 1 - (np.count_nonzero(Yn) / Yn.shape[0])
    
    
