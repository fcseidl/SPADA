#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 23:28:59 2020

@author: fcseidl

Various auxiliary methods for SPADA research.
"""

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


def kMeansClustering(X, k):
    """Wrapper function for sklearn.cluster.KMeans class."""
    return KMeans(n_clusters=k).fit_predict(X)


def maxSilhouetteClusters(X, algorithm, max_n_clusters=10):
    """
    Call a clustering algorithm on a dataset trying a range of possible 
    numbers of clusters, and determine which number of clusters gives the 
    highest average silhouette score.

    Parameters
    ----------
    X : array, (n_samples, n_features)
        Data matrix.
    algorithm : callable
        Clustering algorithm to use. algorithm(X, k) returns an array 
        containing labels dividing the samples into k clusters.
    max_n_clusters : int, optional
        Try forming every number of clusters in [2, max_n_clusters). Return 
        the clustering with the best silhouette score. The default is 10.

    Returns
    -------
    n_clusters : int
        Number of clusters used.
    labels : array, shape (n_samples,)
        Labels of each point.

    """
    print("Clustering data with 2 thru", max_n_clusters, "clusters...")
    n_clusters = 2
    labels_best = algorithm(X, 2)
    score_best = silhouette_score(X, labels_best)
    for k in range(3, max_n_clusters):
        labels = algorithm(X, k)
        score = silhouette_score(X, labels)
        if score > score_best:
            n_clusters = k
            labels_best = labels
            score_best = score
    print(n_clusters, "clusters chosen to maximize average silhouette score")
    return n_clusters, labels_best
    

def dropoutRate(Yn):
    """
    Compute fraction of dropouts for a gene.
    """
    return 1 - (np.count_nonzero(Yn) / Yn.shape[0])
    
    
