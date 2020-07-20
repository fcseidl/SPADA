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
    A single sample from the categorical distribution.
    
    """
    counts = np.random.multinomial(1, p)
    return list(counts).index(1)
    

def randomA(N, K):
    """
    Generate random N x K signature matrix.

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


def randomalpha(K):
    """
    Generate random paramter for Dirichlet distribution.

    Parameters
    ----------
    K : int
        Number of cells types, length of alpha.

    Returns
    -------
    alpha : array, shape (K,)

    """
    return np.random.rand(K) * 1e5
    

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
        A = randomA(N, K)
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
                if np.random.rand() < np.e ** (-lam * (np.log(Y[n, l]) ** 2)):
                    Y[n, l] = 0


def simulateJointData(N=273, M=72, L=213, K=3, lam=0.1, alpha=None, A=None):
    """
    Create a simulated joint dataset. Default settings resemble the real 
    dataset used in the URSM paper.

    Parameters
    ----------
    N : int, optional
        Number of genes. The default is 273.
    M : int, optional
        Number of bulk samples. The default is 72.
    L : int, optional
        Number of single cells. The default is 213.
    K : int, optional
        Number of cell types. The default is 3.
    lam : float, optional
        Paramter to double exponential dropout model. The default is 0.1.
    alpha : array, shape (K,) optional
        Parameter for Dirichlet distribution. Randomly generated if not given.
    A : array, shape (N, K), optional
        Signature matrix. Randomly generated if not given.

    Returns
    -------
    X : array, shape (N, M)
        bulk data
    Y : array shape (N, L)
        single-cell data
        
    """
    if A is None:
        A = randomA(N, K)
    else:
        N, K = A.shape
    if alpha is None:
        alpha = randomalpha(K)
    X = bulk(N, M, K, alpha, A=A)
    Y = true_single_cell(N, L, K, alpha, A=A)
    doubleExpDropouts(Y, lam)
    return X, Y


def g_star(x, c):
    """non-linear monotone tranformation from MMC paper"""
    return 1 / (1 + np.e ** (-c * x))
    