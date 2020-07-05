#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 16:07:51 2020

@author: fcseidl

This module contains various tests.
"""

import numpy as np
from sklearn.decomposition import PCA, FactorAnalysis
import matplotlib.pyplot as plt
import csv

import preprocessing
import simulations as sims
import HypothesisTesting as ht
import SPADAutil as util

# how much variance is explained by principal components?
if 1:
    print("Dataset: pancreatic islets single cell")
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_sc.csv"
    
    print("Reading data matrix Y...")
    Y = preprocessing.csvToMatrix(scfile)
    
    print("Performing PCA...")
    pca = PCA(n_components=0.9)   # Explain at least 90% of variance
    pca.fit(Y.T)
    
    print("Inferred dimensionality:", pca.components_.shape[0])
    print("Explained percentages of variance:", pca.explained_variance_ratio)
    
# how much variance is explained by each gene?
if 0:
    print("test not implemented")

# hypothesis testing with 2 real datasets
if 0:
    print("Bulk: pancreatic islets")
    print("Single-cell: pancreatic islets")
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_sc.csv"
    ht.identifyJointDatasets(bulkfile, scfile)
    print()
    
    print("Bulk: 3cl mixture")
    print("Single-cell: 3cl mixture")
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_sc.csv"
    ht.identifyJointDatasets(bulkfile, scfile)
    print()
    
    print("Bulk: 3cl mixture")
    print("Single-cell: pancreatic islets")
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_islets_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_islets_sc.csv"
    ht.identifyJointDatasets(bulkfile, scfile)
    print()
    
    print("Bulk: pancreatic islets")
    print("Single-cell: 3cl mixture")
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_3cl_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_3cl_sc.csv"
    ht.identifyJointDatasets(bulkfile, scfile)
    print()

# hypothesis testing with simulated data
if 0:
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X1, Y1 = sims.simulateURSMdata()
    X2, Y2 = sims.simulateURSMdata()
    
    print('Are X1 and Y1 joint? Expect 0, receive', 
          hp.residualsToCone(X1, Y1))
    print('Are X2 and Y1 joint? Expect > 0, receive', 
          hp.residualsToCone(X2, Y1))
    
    print("Probability of X1 and Y1 under null hypothesis <=",
          hp.pvalue(X1, Y1))
    print("Probability of X2 and Y1 under null hypothesis <=",
          hp.pvalue(X2, Y1))

# hypothesis testing on dimensionality-reduced simulated data with PCA
if 0:
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X1, Y1 = sims.simulateURSMdata()
    X2, Y2 = sims.simulateURSMdata()
    
    print("Performing PCA on single-cell data...")
    D = 200
    pca = PCA(n_components=D)
    Y1_hat = pca.fit_transform(Y1.T).T
    
    print("Transforming bulk data into feature space...")
    X1_hat = pca.transform(X1.T).T
    X2_hat = pca.transform(X2.T).T
    
    print('Are X1 and Y1 joint? Expect 0, receive', 
          ht.residualsToCone(X1_hat, Y1_hat))
    print('Are X2 and Y1 joint? Expect > 0, receive', 
          ht.residualsToCone(X2_hat, Y1_hat))
    
    print("Probability of X1 and Y1 under null hypothesis <=",
          ht.pvalue(X1_hat, Y1_hat))
    print("Probability of X2 and Y1 under null hypothesis <=",
          ht.pvalue(X2_hat, Y1_hat))
    
# hypothesis testing on dimensionality-reduced simulated data with FA
if 0:
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X1, Y1 = sims.simulateURSMdata()
    X2, Y2 = sims.simulateURSMdata()
    
    print("Performing Factor Analysis on single-cell data...")
    D = 25
    fa = FactorAnalysis(n_components=D)
    Y1_hat = fa.fit_transform(Y1.T).T
    
    print("Transforming bulk data into feature space...")
    X1_hat = fa.transform(X1.T).T
    X2_hat = fa.transform(X2.T).T
    
    print('Are X1 and Y1 joint? Expect 0, receive', 
          ht.residualsToCone(X1_hat, Y1_hat))
    print('Are X2 and Y1 joint? Expect > 0, receive', 
          ht.residualsToCone(X2_hat, Y1_hat))
    
    print("Probability of X1 and Y1 under null hypothesis <=",
          ht.pvalue(X1_hat, Y1_hat))
    print("Probability of X2 and Y1 under null hypothesis <=",
          ht.pvalue(X2_hat, Y1_hat))

# hypothesis testing with two halves of a real dataset
if 0:
    np.random.seed(0)
    
    # NOTE: these are local files which are unavailable on other machines
    
    # 3 cell line mixture
    bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_bulk.csv"
    scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_3cl_sc.csv"
    
    # pancreatic islets data
    #bulkfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_bulk.csv"
    #scfile = "/Users/fcseidl/Documents/SPADA/SPADA/datasets/ssf_islets_sc.csv"
    
    print("Reading full data matrices X and Y...")
    X = preprocessing.csvToMatrix(bulkfile)
    Y = preprocessing.csvToMatrix(scfile)
    
    N = X.shape[0]
    assert(Y.shape[0] == N)
    print("Number of genes:", N)
    print("Number of bulk samples:", X.shape[1])
    print("Number of single cells:", Y.shape[1])
    
    print("Shuffling X alongside Y...")
    indices = np.arange(N)
    np.random.shuffle(indices)
    X = X[indices]
    Y = Y[indices]
    
    print("Constructing half data matrices X1, X2, Y1, Y2...")
    halfway = int(N / 2)
    X1 = X[:halfway]
    X2 = X[halfway:]
    Y1 = Y[:halfway]
    Y2 = Y[halfway:]

    print("Are X1 and Y1 joint? Expect 0, receive",
          ht.residualsToCone(X1, Y1))
    print("Are X1 and Y2 joint? Expect > 0, receive",
          ht.residualsToCone(X1, Y2))
    print("Are X2 and Y1 joint? Expect > 0, receive",
          ht.residualsToCone(X2, Y1))
    print("Are X2 and Y2 joint? Expect 0, receive",
          ht.residualsToCone(X2, Y2))
    
    print("Probability of X1 and Y1 under null hypothesis <=",
          ht.pvalue(X1, Y1))
    print("Probability of X1 and Y2 under null hypothesis <=",
          ht.pvalue(X1, Y2))
    print("Probability of X2 and Y1 under null hypothesis <=",
          ht.pvalue(X2, Y1))
    print("Probability of X2 and Y2 under null hypothesis <=",
          ht.pvalue(X2, Y2))

# perform ZIFA on simulated data
if 0:
    from ZIFA import ZIFA
    from preprocessing import ZIFApreprocessing
    from SPADAutil import marker_quality
    np.set_printoptions(precision=3)
    
    N = 70
    M = 20
    L = 100
    K = 4
    lam = 0.1
    alpha = [ 1 for _ in range(K) ]
    #alpha = [9, 16, 0.3, 6]
    
    A = sims.randomA(N, K)
    X = sims.bulk(N, M, K, alpha, A=A)
    Y = sims.true_single_cell(N, L, K, alpha, A=A)
    print("noiseless bulk data:\n", X)
    print("noiseless single-cell data, no dropouts:\n", Y)
    sims.doubleExpDropouts(Y, lam)
    print("single-cell data after dropouts:\n", Y)
    
    Y = ZIFApreprocessing(Y)
    Z, params = ZIFA.fitModel(Y, K)
    print("ZIFA estimated latent positions:\n", Z)
    
    '''
    mq = marker_quality(A)
    plt.scatter(
        [ max(mqn) for mqn in mq ],
        [ np.var(Xn) for Xn in X ]
        )
    plt.xlabel('marker quality')
    plt.ylabel('variance across samples')
    '''
    
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
        
        print("Scaled solid angle of original cone ~=", 
              util.scaledSolidAngle(Y))
        print("SSA of cone in PCA feature space ~=", 
              util.scaledSolidAngle(Y_hat))
    
    # sample from 2d unit sphere first quadrant
    if 0:
        points = util.uniformFromUnitSphere(2, 200)
        
        # reflect to positive orthant
        points = np.abs(points)
        
        # 2d case
        plt.scatter(
            points[:, 0],
            points[:, 1],
            )
    
    # sample from 3d unit sphere positive orthant
    if 0:
        points = util.uniformFromUnitSphere(3, 200)
        
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
    
# Plot description of variances of gene expressions in real scRNA-seq data
if 0:
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
    print("Generating joint datasets X1, Y1, and X2, Y2...")
    X, Y = sims.simulateURSMdata()
    
    print("Performing PCA on single-cell data...")
    D = 206
    pca = PCA(n_components=D)
    Y_hat = pca.fit_transform(Y.T).T
    
    print("Scaled solid angle of original cone ~=", 
          util.scaledSolidAngle(Y))
    print("SSA of cone in PCA feature space ~=", 
          util.scaledSolidAngle(Y_hat))

# sample from 2d unit sphere first quadrant
if 0:
    points = util.uniformFromUnitSphere(2, 200)
    
    # reflect to positive orthant
    points = np.abs(points)
    
    # 2d case
    plt.scatter(
        points[:, 0],
        points[:, 1],
        )

# sample from 3d unit sphere positive orthant
if 0:
    points = util.uniformFromUnitSphere(3, 200)
    
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

# simple test for findBlockExpressingNGenes()
if 0:
    data = np.array([[0, 1, 0, 2, 0, 3],
                     [0, 0, 1, 2, 0, 3],
                     [0, 1, 0, 0, 0, 0],
                     [0, 1, 0, 0, 0, 3]])
    print(preprocessing.findBlockExpressingNGenes(data, 1, 3))

# try to find 72 12-13 PCW fetal prefrontal cortical cells in Kang data
if 0:
    filename = (
        "/Users/fcseidl/Downloads/GSE25219_DABG_pvalue.csv")
    with open(filename) as readfile:
        read_tsv = csv.reader(readfile, delimiter=",")
        data = np.array([ row for row in read_tsv ])
    data = data[1:, 1:]             # trim first row and col
    data = data.astype(np.float)    # convert to floats
    print(preprocessing.findBlockExpressingNGenes(data, 72, 16947))

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
