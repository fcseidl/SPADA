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
    Prepare scRNA-seq data matrix for ZIFA dimensionality reduction. Note that 
    the preprocessed data matrix must still be transposed for ZIFA to give a 
    meaningful result.

    Parameters
    ----------
    Y : array
        N x L single-cell data matrix

    Returns
    -------
    Y : array
        Y with zero rows removed. If Y contains raw counts, then Y is also 
        log-transformed.
    
    """
    print("Preprocessing data for ZIFA dimensionality reduction...")
    # remove zero cols
    N = Y.shape[0]
    n = 0
    remove = []
    for n in range(N):
        if Y[n].sum() < 1e-6:
            remove.append(n)
    Y = np.delete(Y, remove, axis=0)
    if len(remove) > 0:
        print("Removed", len(remove), "all-zero rows from data.")
        
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
    A 2d numpy float array containing data from file, with no title row or 
    column.
    
    """
    with open(filename) as readfile:
        read_csv = csv.reader(readfile, delimiter=delim)
        result = np.array([ row for row in read_csv ])
    result = result[1:, 1:]
    result = result.astype(np.float)
    return result


def csvToNumpy(filename, delim=','):
    """
    Read a rectangular csv file into a numpy array.

    Parameters
    ----------
    filename : string
        name of csv file
    delim : character, optional
        delimiting character, comma by default

    Returns
    -------
    A 2d numpy string array containing file data.

    """
    with open(filename) as readfile:
        read_csv = csv.reader(readfile, delimiter=delim)
        result = np.array([ row for row in read_csv ])
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
    print("Removed", len(indices), "out of", N, "rows.")
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

