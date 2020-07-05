#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:08:48 2020

@author: fcseidl

Tests to decide whether bulk and scRNA-seq datasets are joint.
"""

import numpy as np
import preprocessing
from scipy.optimize import  nnls


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
    print("---Loading datasets---")
    
    print("Reading data matrices X and Y...")
    X = preprocessing.csvToMatrix(bulkfile)
    Y = preprocessing.csvToMatrix(scfile)
    
    N = X.shape[0]
    assert(Y.shape[0] == N)
    print("Number of genes:", N)
    print("Number of bulk samples:", X.shape[1])
    print("Number of single cells:", Y.shape[1])
    
    print("---Performing hypothesis testing---")
    
    print("Removing genes with low variance...")
    def lowVariance(Yn):
        return np.var(Yn) < 80  # TODO: magic number
    X, Y = preprocessing.removeRowsPred(X, Y, lowVariance)
    
    N = X.shape[0]
    assert(Y.shape[0] == N)
    print("Number of remaining genes:", N)
    
    print("Scaling genes by variance...")
    X, Y = preprocessing.scaleRowsByVariance(X, Y)
    
    print("Estimating bound for probability of residuals under null hypothesis...")
    p = pvalue(X, Y)
    print("p <=", p)
    
    return p

    