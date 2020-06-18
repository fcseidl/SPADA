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


if __name__ == "__main__":
    import simulations as sims
    
    N = 15
    M = 8
    L = 80
    K = 4
    alpha = [ 10 for _ in range(K) ]
    A1 = sims.randomA(N, K)
    A2 = sims.randomA(N, K)
    X1 = sims.bulk(N, M, K, alpha, A=A1)
    X2 = sims.bulk(N, M, K, alpha, A=A2)
    Y1 = sims.true_single_cell(N, L, K, alpha, A=A1)
    
    print('Are X1 and Y1 joint? Expect True, receive', 
          HullContainment(X1, Y1))
    print('Are X2 and Y1 joint? Expect False, receive', 
          HullContainment(X2, Y1))
    
    