#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 23:28:59 2020

@author: fcseidl

Various auxiliary methods for SPADA research.
"""

import numpy as np


def ZIFApreprocessing(Y):
    """
    Prepare scRNA-seq data matrix for ZIFA dimensionality reduction.

    Parameters
    ----------
    Y : array
        N x L single-cell data matrix

    Returns
    -------
    Y with zero columns removed. If Y contains raw counts, then Y is also 
    log-transformed.
    """
    # remove zero cols
    L = L0 = Y.shape[1]
    l = 0
    while l < L:
        if Y[:, l].sum() < 1e-6:
            Y[:, l] = Y[:, -1]  # replace zero col with last col
            Y = Y[:, :-1]       # cut off duplicate last col
            L -= 1
        else:
            l += 1
    if L < L0:
        print("Removed", L0 - L, "all-zero columns from data.")
        
    # log-transform Y if data is raw
    if (Y - np.array(Y, dtype='int32')).sum() < 1e-6:
        print("Data was entirely integral. Assumed this was because",
              "log-transformation had not yet been applied and",
              "log-transformed data.")
        Y = np.log(Y + 1)   # lazy avoid division by zero
    
    return Y
        

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
    