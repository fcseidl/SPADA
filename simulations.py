#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 12:00:08 2020

@author: fcseidl

Simulate RNA-seq data for SPADA lab research.
"""

import numpy as np
import statistics as stats
from sklearn.preprocessing import normalize


def bulk(N, M, K, depth, alpha, A=None):
    """
    Simulate bulk RNA-seq data according to Dirichlet-multinomial model from 
    URSM.
    
    Parameters
    ----------
    N : number of genes
    M : number of samples
    K : number of cell types
    depth : sequencing depth. Assumed to be constant across samples.
    alpha : parameter for Dirichlet distribution
    A : array, optional
        N x K normalized signature matrix. Random if not supplied.

    Returns
    -------
    X : array
        N x M data matrix giving expression level of each gene in each sample.
    """
    if A is None:
        A = np.random.rand(N, K)
        for k in range(K):
            A[:, k] /= sum(A[:, k])
    X = []
    for j in range(M):
        wj = np.random.dirichlet(alpha)
        xj = np.random.multinomial(depth, A.dot(wj))
        X.append(xj)
    return np.array(X).T


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


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    np.set_printoptions(precision=3)
    
    N = 100
    M = 40
    K = 8
    depth = 600
    alpha = [ 1 for _ in range(K) ]
    #alpha = [9, 16, 0.3, 6]
    
    A = np.random.rand(N, K)
    for k in range(K):
        A[:, k] /= sum(A[:, k])
    
    mq = marker_quality(A)
    X = bulk(N, M, K, depth, alpha, A=A)
    
    plt.scatter(
        [ max(mqn) for mqn in mq ],
        [ np.var(Xn) for Xn in X ]
        )
    plt.xlabel('marker quality')
    plt.ylabel('variance across samples')
    
    