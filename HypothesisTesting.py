#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:08:48 2020

@author: fcseidl

Tests to decide whether bulk and scRNA-seq datasets are joint.
"""

import numpy as np
from sklearn.preprocessing import normalize
from sklearn.decomposition import FactorAnalysis
from scipy.spatial import ConvexHull
from scipy.optimize import linprog, nnls


def HullContainment(X, Y):
    """
    NOTE: false negatives are too likely with this method.
    
    Check that bulk and scRNA-seq data matrices are joint by verifying that 
    bulk profiles are inside convex hull of single-cell profiles.

    Parameters
    ----------
    X : array, shape (N, M)
        bulk data matrix
    Y : array, shape (N, L)
        single-cell data matrix

    Returns
    -------
    Whether normalized columns of X are convex combinations of normalized 
    columns of Y.
    """
    M = X.shape[1]
    # preprocessing: tranpose and normalize
    X = normalize(X.T, norm='l1')
    Y = normalize(Y.T, norm='l1')
    # convex hull should have only vertices from Y
    points = np.r_[X, Y]
    hull = ConvexHull(points, incremental=False, qhull_options='QJ QbB')
    vertices = np.unique(hull.simplices.ravel())
    return (vertices >= M).all()


def fractionInsideHull(X, Y):
    """
    NOTE: fractionInsideCone is more efficient and less likely to produce
    OptimizeWarnings.
    
    Count fraction of normalized bulk RNA-seq profiles which are convex 
    combinations of normalized single cell profiles.

    Parameters
    ----------
    X : array, shape (N, M)
        bulk data matrix
    Y : array, shape (N, L)
        single-cell data matrix

    Returns
    -------
    Fraction of normalized columns of X which are inside convex hull of 
    normalized Y columns.
    """
    M = X.shape[1]
    L = Y.shape[1]
    # preprocessing: normalization
    X = normalize(X, norm='l1', axis=0)
    Y = normalize(Y, norm='l1', axis=0)
    count = 0
    for j in range(M):  # see if X[:, j] is a convex combination of Y cols
        c = np.zeros(L)
        # add additional equality constraint to make combination convex
        # (this constraint is redundant)
        A_equation = np.r_[Y, np.ones((1, L))]
        b_equation = np.r_[X[:, j], np.ones(1)]
        lp = linprog(c, A_eq=A_equation, b_eq=b_equation)
        if lp.success: count += 1
    return count / M


def fractionInsideCone(X, Y):
    """
    Count fraction of bulk RNA-seq profiles which are conical combinations of
    single cell profiles.

    Parameters
    ----------
    X : array, shape (N, M)
        bulk data matrix
    Y : array, shape (N, L)
        single-cell data matrix

    Returns
    -------
    Fraction of columns of X which are inside conical hull of columns of Y.
    """
    M = X.shape[1]
    L = Y.shape[1]
    count = 0
    for j in range(M):  # see if X[:, j] is a conical combination of Y cols
        c = np.zeros(L)
        lp = linprog(c, A_eq=Y, b_eq=X[:, j])
        if lp.success: count += 1
    return count / M


def residualsToCone(X, Y):
    """
    Compute bulk profiles' average normalized distance to conical hull.

    Parameters
    ----------
    X : array, shape (N, M)
        bulk data matrix
    Y : array, shape (N, L)
        single-cell data matrix

    Returns
    -------
    Average over j of min ||X[:, j] - b|| / ||X[:, j]|| for b in cone(Y).
    """
    M = X.shape[1]
    normalized_residuals = []
    for j in range(M):
        _, norm = nnls(Y, X[:, j])
        normalized_residuals.append(norm / np.linalg.norm(X[:, j]))
    return sum(normalized_residuals) / M
        

if __name__ == "__main__":
    # test with 3 cell line mixture data
    if 1:
        import preprocessing
        np.random.seed(44)
        
        # local files, not available on other machines
        bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_MIX3cl_bulkESET.csv"
        scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_MIX3cl_scESET.csv"
        
        print("Reading full data matrices X and Y...")
        X = preprocessing.csvToMatrix(bulkfile)
        Y = preprocessing.csvToMatrix(scfile)
        
        N = X.shape[0]
        assert(Y.shape[0] == N)
        print("Number of genes:", N)
        print("Number of bulk samples:", X.shape[1])
        print("Number of single cells:", Y.shape[1])
        
        print("Shuffling X and Y together...")
        indices = np.arange(N)
        np.random.shuffle(indices)
        X = X[indices]
        Y = Y[indices]
        
        print("Constructing half data matrices X1, X2, Y1, Y2...")
        halfway = int(N / 2)
        X1 = X[:halfway]
        X2 = X[halfway:]
        Y1 = Y[:halfway]
        Y2 = Y[halfway:]
    
        print("Are X1 and Y1 joint? Expect 0, receive",
              residualsToCone(X1, Y1))
        print("Are X1 and Y2 joint? Expect > 0, receive",
              residualsToCone(X1, Y2))
        print("Are X2 and Y1 joint? Expect > 0, receive",
              residualsToCone(X2, Y1))
        print("Are X2 and Y2 joint? Expect 0, receive",
              residualsToCone(X2, Y2))
    
    # test with simulated data
    if 0:
        import simulations as sims
        #np.random.seed(34)
        
        N = 273
        M = 72
        L = 213
        K = 3
        lam = 0.1
        alpha = [ 2e4, 1e4, 7e4 ]
        #alpha = [ 1e3 for _ in range(K) ]  # assumes symmetric Dirichlet prior
        A1 = sims.randomA(N, K)
        A2 = sims.randomA(N, K)
        X1 = sims.bulk(N, M, K, alpha, A=A1)
        X2 = sims.bulk(N, M, K, alpha, A=A2)
        Y1 = sims.true_single_cell(N, L, K, alpha, A=A1)
        sims.doubleExpDropouts(Y1, lam)
        
        print("using residuals of true data:")
        print('Are X1 and Y1 joint? Expect 0, receive', 
              residualsToCone(X1, Y1))
        print('Are X2 and Y1 joint? Expect > 0, receive', 
              residualsToCone(X2, Y1))
        
        '''
        print("using conical hull of true data:")
        print('Are X1 and Y1 joint? Expect 1.0, receive', 
              fractionInsideCone(X1, Y1))
        print('Are X2 and Y1 joint? Expect 0.0, receive', 
              fractionInsideCone(X2, Y1))
        
        model = FactorAnalysis(n_components=K)
        join1 = np.concatenate((X1, Y1), axis=1).T
        Z1 = model.fit_transform(join1).T
        join2 = np.concatenate((X2, Y1), axis=1).T
        Z2 = model.fit_transform(join2).T
        
        print("using conical hull of dimension-reduced data:")
        print('Are X1 and Y1 joint? Expect 1.0, receive', 
              fractionInsideCone(Z1[:, :M], Z1[:, M:]))
        print('Are X2 and Y1 joint? Expect 0.0, receive', 
              fractionInsideCone(Z2[:, :M], Z2[:, M:]))
        '''
    
    