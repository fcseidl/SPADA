#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 23:28:59 2020

@author: fcseidl

Various auxiliary methods for SPADA research.
"""

import numpy as np
from scipy.optimize import linprog
   

def marker_quality(A):
    """
    Parameters
    ----------
    A : array
        Signature matrix for N genes in K cell types.

    Returns
    -------
    mq : array
        N x K matrix whose n,k entry is marker quality of nth gene for kth 
        cell type
    """
    mq = np.zeros(A.shape)
    for n in range(A.shape[0]):
        s = sum(A[n])
        mq[n] = [ A[n][k] / s for k in range(A.shape[1]) ]
    return mq


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
    n_points = int(1e3)  # TODO: magic number
    count = 0
    for _ in range(n_points):
        x = uniformFromUnitSphere(N)
        c = np.zeros(L)
        lp = linprog(c, A_eq=Y, b_eq=x)
        if lp.success: count += 1
    return count / n_points
    
    
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    # Plot description of variances of gene expressions in real scRNA-seq data
    if 1:
        import preprocessing
        
        scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_islets_sc.csv"
        
        print("Reading single-cell data matrix Y...")
        Y = preprocessing.csvToMatrix(scfile)
        
        '''
        n = 400
        indices = set([ np.random.randint(Y.shape[0]) for _ in range(n) ])
        '''
        indices = np.arange(Y.shape[0])
        
        print("Computing variances for", len(indices), "random genes..." )
        V = [ np.var(Y[l]) for l in indices ]
        V.sort()
        
        print("Plotting variances...")
        plt.plot(np.arange(len(V)), V)
        plt.xlabel("gene")
        plt.ylabel("variance")
    
    # assess solid angle of original and PCA cones for simulated single-cell 
    # data
    if 0:
        import simulations as sims
        from sklearn.decomposition import PCA
        
        print("Generating joint datasets X1, Y1, and X2, Y2...")
        X, Y = sims.simulateURSMdata()
        
        print("Performing PCA on single-cell data...")
        D = 206
        pca = PCA(n_components=D)
        Y_hat = pca.fit_transform(Y.T).T
        
        print("Scaled solid angle of original cone ~=", scaledSolidAngle(Y))
        print("SSA of cone in PCA feature space ~=", scaledSolidAngle(Y_hat))
    
    # sample from 2d unit sphere first quadrant
    if 0:
        points = uniformFromUnitSphere(2, 200)
        
        # reflect to positive orthant
        points = np.abs(points)
        
        # 2d case
        plt.scatter(
            points[:, 0],
            points[:, 1],
            )
    
    # sample from 3d unit sphere positive orthant
    if 0:
        points = uniformFromUnitSphere(3, 200)
        
        # reflect to positive orthant
        points = np.abs(points)
        
        # 3d case
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2]
            )