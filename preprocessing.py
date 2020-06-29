#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:58:09 2020

@author: fcseidl

Preprocessing bulk and scRNA-seq data.
"""

import numpy as np
import csv


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


if __name__ == "__main__":
    bulkfile = (
    "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_MIX3cl_bulkESET.csv")
    scfile = (
    "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_MIX3cl__scESET.csv")
      
    # simple test for findBlockExpressingNGenes()
    if 0:
        data = np.array([[0, 1, 0, 2, 0, 3],
                         [0, 0, 1, 2, 0, 3],
                         [0, 1, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 3]])
        print(findBlockExpressingNGenes(data, 1, 3))
    
    # try to find 72 12-13 PCW fetal prefrontal cortical cells in Kang data
    if 0:
        filename = (
            "/Users/fcseidl/Downloads/GSE25219_DABG_pvalue.csv")
        with open(filename) as readfile:
            read_tsv = csv.reader(readfile, delimiter=",")
            data = np.array([ row for row in read_tsv ])
        data = data[1:, 1:]             # trim first row and col
        data = data.astype(np.float)    # convert to floats
        print(findBlockExpressingNGenes(data, 72, 16947))
    
    # begin preprocessing of Camp data by ignoring non-fetal cells 
    if 0:
        filename = (
            "/Users/fcseidl/Downloads/GSE75140_hOrg.fetal.master.data.frame.txt")
        print("Reading data from file...")
        with open(filename) as readfile:
            read_tsv = csv.reader(readfile, delimiter="\t")
            data = np.array([ row for row in read_tsv ])
        print(data)
        print("Removed rows of data not from fetal cells:")
        data = data[-226:]
        print(data)
