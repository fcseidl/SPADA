#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:08:48 2020

@author: fcseidl

Tests to decide whether bulk and scRNA-seq datasets are joint.
"""

import numpy as np
from sklearn.preprocessing import normalize
from scipy.spatial import ConvexHull
from scipy.optimize import linprog


def HullContainment(X, Y):
    """
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


def fractionInsideHull(X, Y):
    """
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


if __name__ == "__main__":
    import simulations as sims
    
    N = 100
    M = 60
    L = 800
    K = 4
    alpha = [ 1e4 for _ in range(K) ]  # assumes symmetric Dirichlet prior
    A1 = sims.randomA(N, K)
    A2 = sims.randomA(N, K)
    X1 = sims.bulk(N, M, K, alpha, A=A1)
    X2 = sims.bulk(N, M, K, alpha, A=A2)
    Y1 = sims.true_single_cell(N, L, K, alpha, A=A1)
    
    print("using convex hull:")
    print('Are X1 and Y1 joint? Expect 1.0, receive', 
          fractionInsideHull(X1, Y1))
    print('Are X2 and Y1 joint? Expect 0.0, receive', 
          fractionInsideHull(X2, Y1))
    
    print("using conical hull:")
    print('Are X1 and Y1 joint? Expect 1.0, receive', 
          fractionInsideCone(X1, Y1))
    print('Are X2 and Y1 joint? Expect 0.0, receive', 
          fractionInsideCone(X2, Y1))
    
    '''
    print('Are X1 and Y1 joint? Expect True, receive', 
          HullContainment(X1, Y1))
    print('Are X2 and Y1 joint? Expect False, receive', 
          HullContainment(X2, Y1))
    '''
    
    