#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 12:00:08 2020

@author: fcseidl

Simulate RNA-seq data for SPADA lab research.
"""

import numpy as np
from sklearn.preprocessing import normalize


eps = 1e-6


def categorical(p):
    """
    A categorical random variable.

    Parameters
    ----------
    p : sequence of nonnegative floats, sum(p) = 1
        p[i] is the probability of returning i.

    Returns
    -------
    A single sample from categorical distribution.
    """
    counts = np.random.multinomial(1, p)
    return list(counts).index(1)
    


def randomA(N, K):
    """
    Generate uniform random N x K signature matrix.

    Parameters
    ----------
    N : number of genes
    K : number of cell types

    Returns
    -------
    A : array
        signature matrix
    """
    A = np.random.rand(N, K)
    A += eps
    A = normalize(A, norm='l1', axis=0)
    return A
    


def bulk(N, M, K, alpha, A=None):
    """
    Simulate bulk RNA-seq data according to Dirichlet-multinomial model from 
    URSM. No noise is induced.
    
    Parameters
    ----------
    N : number of genes
    M : number of samples
    K : number of cell types
    alpha : sequence of floats, length K
        parameter for Dirichlet distribution
    A : array, optional
        N x K normalized signature matrix. Random if not supplied.

    Returns
    -------
    X : array
        N x M data matrix giving expression level of each gene in each sample.
    """
    if A is None:
        A = randomA(N, K)
    X = []
    for j in range(M):
        wj = np.random.dirichlet(alpha)
        depth = np.random.poisson(50 * N)   # sequencing depth
        xj = np.random.multinomial(depth, A.dot(wj))
        X.append(xj)
    return np.array(X).T


def true_single_cell(N, L, K, alpha, A=None):
    """
    Simulate scRNA-seq data according to multinomial model from URSM. No 
    noise or dropout events are induced.

    Parameters
    ----------
    N : number of genes
    L : number of genes
    K : number of cell types
    alpha : sequence of floats, length K
        proportions of each cell type, can be parameter for bulk dirichlet 
        distribution
    A : array, optional
        N x K normalized signature matrix. Random if not supplied.

    Returns
    -------
    Y : array
        N x L data matrix giving expression level of each gene in each cell.
    """
    if A is None:
        A = randomA()
    proportions = np.array(alpha) / sum(alpha)
    Y = []
    for l in range(L):
        Gl = categorical(proportions)  # cell type
        depth = np.random.negative_binomial(2, 
                                            1 - 2 * N / (2 * N + 2))
        yj = np.random.multinomial(depth, A[:, Gl])
        Y.append(yj)
    return np.array(Y).T


def doubleExpDropouts(Y, lam):
    """
    Induce dropouts according to decaying squared exponenetial model as in 
    ZIFA.

    Parameters
    ----------
    Y : array
        N x L true expression levels for single cell data
    lam : float
        double exponential decay coefficient
    
    Effects
    -------
    Randomly zero elements of Y with probability 
    double exponential in the log of true expression.
    """
    N = Y.shape[0]
    L = Y.shape[1]
    for n in range(N):
        for l in range(L):
            if Y[n, l] > 0:
                if np.random.rand() < np.e ** (-lam * np.log(Y[n, l]) ** 2):
                    Y[n, l] = 0
    

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
    
    N = 70
    M = 20
    L = 100
    K = 8
    lam = 0.1
    alpha = [ 1 for _ in range(K) ]
    #alpha = [9, 16, 0.3, 6]
    
    A = randomA(N, K)
    X = bulk(N, M, K, alpha, A=A)
    Y = true_single_cell(N, L, K, alpha, A=A)
    print("noiseless bulk data:\n", X)
    print("noiseless single-cell data, no dropouts:\n", Y)
    doubleExpDropouts(Y, lam)
    print("single-cell data after dropouts:\n", Y)
    
    '''
    mq = marker_quality(A)
    plt.scatter(
        [ max(mqn) for mqn in mq ],
        [ np.var(Xn) for Xn in X ]
        )
    plt.xlabel('marker quality')
    plt.ylabel('variance across samples')
    '''
    
    