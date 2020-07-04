#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:58:09 2020

@author: fcseidl

Preprocessing bulk and scRNA-seq data.
"""

import numpy as np
import csv


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
    while l < L:    # TODO: terrible slow way too do this
        if Y[:, l].sum() < 1e-6:
            Y = np.delete(Y, l, axis=1)
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


def csvToMatrix(filename, delim=','):
    """
    Convert a csv table to a 2d numpy array by removing first row and col
    (containing names), and converting to floats.

    Parameters
    ----------
    filename : string
        name of csv file
    delim : character, optional
        delimiting character, comma by default

    Returns
    -------
    A 2d numpy array containing (unlabeled) data from file.
    """
    with open(filename) as readfile:
        read_csv = csv.reader(readfile, delimiter=delim)
        result = np.array([ row for row in read_csv ])
    result = result[1:, 1:]
    result = result.astype(np.float)
    return result


def removeRowsPred(X, Y, pred):
    """
    Remove rows of bulk and scRNA-seq data matrices according to a predicate.

    Parameters
    ----------
    X : array, shape (N, M)
        bulk data matrix
    Y : array, shape (N, L)
        single-cell data matrix
    pred : function
        Predicate to indicate that a row should be removed

    Returns
    -------
    X, Y modified to exclude rows n where pred(Y[n]) is True.
    """
    N = Y.shape[0]
    indices = []
    for n in range(N):
        if pred(Y[n]): indices.append(n)
    X = np.delete(X, indices, axis=0)
    Y = np.delete(Y, indices, axis=0)
    return X, Y


def scaleRowsByVariance(X, Y):
    """
    Scale rows of bulk and scRNA-seq data matrices to weight genes by variance.

    Parameters
    ----------
    X : array, shape (N, M)
        bulk data matrix
    Y : array, shape (N, L)
        single-cell data matrix

    Returns
    -------
    X, Y, where each row is scaled in proportion to the original variance of 
    that gene in Y.
    """
    # L-inifinity normalized variances
    var = np.array([ np.var(Yn) for Yn in Y ])
    var /= max(var)
    
    # multiply X, Y by diag(var) TODO: could be done faster...
    V = np.diag(var)
    X = V.dot(X)
    Y = V.dot(Y)
    
    return X, Y
            

def findBlockExpressingNGenes(data, n_samps, n_expressed):
    """
    Search data for contiguous block of n_samps samples which express exactly 
    n_expressed genes collectively.

    Parameters
    ----------
    data : array, shape (N, M)
        Data matrix for N genes and M samples.
    n_samps : int
        Desired number of samples
    n_expressed : TYPE
       Desired number of genes expressed in samples.

    Returns
    -------
    If a satisfactory block exists, the index of the first sample in the block,
    otherwise -1.
    
    NOTE: slow, complexity O(N * M * n_samps).
    """
    N, M = data.shape
    assert n_samps <= M and n_expressed <= N
    for begin in range(M - n_samps + 1):
        block = data[:, begin:begin + n_samps]
        nonzero = block[np.any(block, axis=1)]
        if nonzero.shape[0] == n_expressed:
            return begin
    return -1

