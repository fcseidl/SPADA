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
from scipy.optimize import linprog, nnls
from sklearn.decomposition import PCA, FactorAnalysis


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

    