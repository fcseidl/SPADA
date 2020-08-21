#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 21:57:40 2020

@author: fcseidl

Functions to generate data and visualizations for tech report.
"""

import numpy as np
from sklearn.preprocessing import normalize

from BestMap import BestMap

import preprocessing
import simulations as sims
import RelatednessTesting as rt
import SPADAutil as util


def mouse_brain_test(clusterer):
    n_clusters = [2, 5, 10, 20] # various numbers of clusters to use
    
    print("Loading Chen data matrix...")
    chen = preprocessing.csvToMatrix(
        "/Users/fcseidl/Documents/SPADA/RNAseq/mouse_brain/ssf_chen.tsv", 
        delim='\t')
    
    print("Loading Campbell data matrix...")
    campbell = preprocessing.csvToMatrix(
        "/Users/fcseidl/Documents/SPADA/RNAseq/mouse_brain/ssf_campbell.tsv", 
        delim='\t')
    
    print("Removing sparse genes from both datasets...")
    def sparse(Yn):
        return util.dropoutRate(Yn) > 0.8
    chen, campbell = preprocessing.removeRowsPred(chen, campbell, sparse)
    campbell, chen = preprocessing.removeRowsPred(campbell, chen, sparse)
    
    print("Randomly choosing subset of data...")
    L = 120
    chen = np.random.permutation(chen)
    campbell = np.random.permutation(campbell)
    chen = chen[:, :L]
    campbell = campbell[:, :L]
    
    print("Applying cosine normalization to samples...")
    chen = normalize(chen)
    campbell = normalize(campbell)
    
    for k in n_clusters:
        print("Using", k, "clusters...")
        
        print("Clustering Chen data alone...")
        chen_lbls = clusterer(chen.T, k)
        
        print("Clustering Campbell data alone...")
        campbell_lbls = clusterer(campbell.T, k)
        
        print("Clustering both datasets together...")
        joint = np.concatenate((chen, campbell), axis=1)
        joint_lbls = clusterer(joint.T, k)
        chen_joint_lbls = BestMap(chen_lbls, joint_lbls[:L])
        campbell_joint_lbls = BestMap(campbell_lbls, joint_lbls[L:])
        
        print("breakpoint passed")


if __name__ == "__main__":
    mouse_brain_test(util.kMeansClustering)
    
    