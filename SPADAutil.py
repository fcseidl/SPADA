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
    
    
    
