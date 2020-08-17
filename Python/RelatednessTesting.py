#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:08:48 2020

@author: fcseidl

Tests to decide whether RNA-seq datasets are joint.
"""

import numpy as np
import preprocessing
import SPADAutil as util
from scipy.optimize import nnls
from sklearn.cluster import KMeans

import matplotlib.pyplot as plt


def residualsToCone(X, Y):
    """
    Compute bulk profiles' average normalized l2 distance to conical hull.

    Parameters
    ----------
    X : array, shape (N, M)
        bulk data matrix
    Y : array, shape (N, L)
        single-cell data matrix

    Returns
    -------
    Average over j of min ||X[:, j] - b|| / ||X[:, j]|| for b in cone(Y). 
    l2 norm is used.
    
    """
    M = X.shape[1]
    normalized_residuals = []
    for j in range(M):
        _, norm = nnls(Y, X[:, j])
        normalized_residuals.append(norm / np.linalg.norm(X[:, j]))
    return sum(normalized_residuals) / M


def pvalue(X, Y):
    """
    Compute a bound for the probability under the null hypothesis that 
    bulk and scRNA-seq data are totally unrelated.

    Parameters
    ----------
    X : array, shape (N, M)
        bulk data matrix
    Y : array, shape (N, L)
        single-cell data matrix

    Returns
    -------
    Chebyshev bound for probability under null hypothesis.
    
    """
    true_residuals = residualsToCone(X, Y)   # true order
    permuted_residuals = []
    for n in range(100):     # TODO: magic number here...
        X_permuted = np.random.permutation(X)
        permuted_residuals.append(residualsToCone(X_permuted, Y))
    exp = np.mean(permuted_residuals)
    a = np.abs(exp - true_residuals)
    var = np.var(permuted_residuals)
    return var / (a ** 2)   # Chebyshev bound


def identifyJointDatasets(bulkfile, scfile, delim=',', quiet=False):
    """
    Read bulk and scRNA-seq data from csv files and calculate p-value.
    
    REQUIRES: rows correspond to samples, columns to genes. The first 
    row and column contain gene/sample names or numbers. Bulk and single-cell 
    files are preprocessed to contain the same genes in the same order.
    The data are raw counts, not log-transformed.

    Parameters
    ----------
    bulkfile : path to bulk data file
    scfile : path to single-cell data file
    delim : character, optional
        Separating character in file format. The default is ','.
    quiet : bool, optional
        True prevents additional messages from printing. The default is False.

    Returns
    -------
    Chebyshev bound for probability that conic residuals are below their 
    observed values under null hypothesis.
    
    """
    # conditional print function
    def printIf(*args):
        if not quiet:
            print(*args)
        
    printIf("---Loading datasets---")
    
    printIf("Reading data matrices X and Y...")
    X = preprocessing.csvToMatrix(bulkfile, delim)
    Y = preprocessing.csvToMatrix(scfile, delim)
    
    N = X.shape[0]
    assert(Y.shape[0] == N)
    printIf("Number of genes:", N)
    printIf("Number of bulk samples:", X.shape[1])
    printIf("Number of single cells:", Y.shape[1])
    
    printIf("---Performing hypothesis testing---")
    
    printIf("Removing genes with low variance...")
    def lowVariance(Yn):
        return np.var(Yn) < 80  # TODO: magic number
    X, Y = preprocessing.removeRowsPred(X, Y, lowVariance)
    
    N = X.shape[0]
    assert(Y.shape[0] == N)
    printIf("Number of remaining genes:", N)
    
    printIf("Scaling genes by variance...")
    X, Y = preprocessing.scaleRowsByVariance(X, Y)
    
    printIf("Estimating bound for probability of residuals under null hypothesis...")
    p = pvalue(X, Y)
    printIf("p <=", p)
    
    return p


def clusterHeterogeneity(A, B, clusterer):
    """
    Assess relatedness of two data matrices using cluster heterogeneity 
    technique.

    Parameters
    ----------
    A : array, shape (n_feat, n_samp_A)
        First data matrix.
    B : array, shape (n_feat, n_samp_B)
        Second data matrix.
    clusterer : callable
        Function which implements a clustering algorithm. Receives data 
        matrix of shape (n_samples, n_features) as input, and returns 
        n_clusters, labels.

    Returns
    -------
    A relatedness score which is close to 0.25 if the clustering algorithm is
    unable to distinguish the two datasets, and close to 0 if it is.
    
    """
    n_feat, n_samp_A = A.shape
    n_samp_B = B.shape[1]
    assert(B.shape[0] == n_feat) # both datasets have same number of features
    AB = np.concatenate((A, B), axis=1) # joint data matrix
    
    # perform clustering (must transpose joint matrix)
    n_clusters, labels = clusterer(AB.T)
    
    # Identify cluster of each datapoint. 
    # CA[k] is the number of points from cluster k in A.
    # CB[k] is the number of points from cluster k in B. 
    # C[k] is the size of cluster k.
    CA = np.zeros(n_clusters, dtype=int)
    CB = np.zeros(n_clusters, dtype=int)
    C = np.zeros(n_clusters, dtype=int)
    for k in labels[:n_samp_A]:        # labels of samples of A
        CA[k] += 1
        C[k] += 1
    for k in labels[n_samp_A:]:       # labels of samples of B
        CB[k] += 1
        C[k] += 1
    
    # overall heterogeneity
    h = n_samp_A * n_samp_B / (n_samp_A + n_samp_B)**2
    
    # cluster heterogeneities
    H = [ CA[k] * CB[k] / C[k]**2 for k in range(n_clusters) ]
    
    # unified cluster heterogeneity score
    #h_star = 1 / (M + N) * sum([ C[k] * H[k] for k in range(K) ])
    h_star = sum(H) / n_clusters
    
    print("Unified cluster heterogeneity score =", h_star)
    print("Overall heterogeneity =", h)
    
    ind = np.arange(n_clusters) 
    clusters = [ k for k in range(n_clusters) ]
    width = 0.35       
    plt.bar(ind, CA, width, label='A')
    plt.bar(ind + width, CB, width,
        label='B')
    
    plt.ylabel('Number of points from dataset')
    plt.xlabel('Cluster')
    plt.title('Cluster composition')
    
    plt.xticks(ind + width / 2, clusters)
    plt.legend(loc='best')
    plt.show()
    
    return h_star / h
    
    
    
    
    