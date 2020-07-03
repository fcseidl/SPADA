#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 23:28:59 2020

@author: fcseidl

Various auxiliary methods for SPADA research.
"""

import numpy as np


def recoverMatrix(Y, MY):
    """
    Recover the matrix of a linear transformation from a spanning set and its 
    image.

    Parameters
    ----------
    Y : array, shape (N, L)
        Columns span R^N
    MY : array, shape (D, L) 
        Image of Y under unkown matrix M of shape (D, N)

    Returns
    -------
    The unknown matrix M.
    """
    
        

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