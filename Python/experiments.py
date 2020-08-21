#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 21:57:40 2020

@author: fcseidl

Functions to generate data and visualizations for tech report.
"""

import numpy as np
from sklearn.preprocessing import normalize
from sklearn.metrics import normalized_mutual_info_score
import matplotlib.pyplot as plt

from BestMap import BestMap
from SSC import sparseSubspaceClustering

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
    
    '''
    _, chen = sims.simulateJointData()
    _, campbell = sims.simulateJointData()
    '''
     
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
        
        print("NMI between chen separate and joint labels =",
              normalized_mutual_info_score(chen_lbls, joint_lbls[:L]))
        
        print("NMI between campbell separate and joint labels =",
              normalized_mutual_info_score(campbell_lbls, joint_lbls[L:]))
        
        print("Plotting cluster composition bar chart...")
        joint_distribution = np.zeros(k)
        chen_distribution = np.zeros((k, k))
        for c,j in zip(chen_lbls, joint_lbls[:L]):
            chen_distribution[int(c)][int(j)] += 1
            joint_distribution[j] += 1
        campbell_distribution = np.zeros((k, k))
        for c,j in zip(campbell_lbls, joint_lbls[L:]):
            campbell_distribution[int(c)][int(j)] += 1
            joint_distribution[j] += 1
        
        '''
        fig, ax = plt.subplots()
        
        size = 0.3
        cmap = plt.get_cmap("tab20c")
        #outer_colors = cmap(np.arange(k))
        #inner_colors = np.concatenate([ outer_colors for _ in range(k) ])
        
        ax.pie(joint_distribution, radius=1, #colors=outer_colors,
               wedgeprops=dict(width=size, edgecolor='w'))
        
        ax.pie(joint_distribution.flatten(), radius=1-size, #colors=inner_colors,
               wedgeprops=dict(width=size, edgecolor='w'))
        
        ax.set(aspect="equal", title='cluster composition')
        plt.show()
        '''
        
        index = np.arange(k)
        ticks = np.arange(k)
        width = 0.35
        y_offset = np.zeros(k)
        
        plt.bar(index + width, joint_distribution, width, label=("joint"))
        
        for c in range(k):
            plt.bar(index, chen_distribution[c], width, bottom=y_offset, 
                    label=("chen %d" % c))
            y_offset += chen_distribution[c]
            plt.bar(index, campbell_distribution[c], width, bottom=y_offset,
                    label=("campbell %d" % c))
            y_offset += campbell_distribution[c]
            
        plt.ylabel('number of points from dataset')
        plt.xlabel('joint cluster')
        plt.title('Cluster Composition')
        plt.xticks(index + width / 2, ticks)
        plt.legend(loc='best')
        plt.show()
        

if __name__ == "__main__":
    mouse_brain_test(sparseSubspaceClustering)
    #mouse_brain_test(util.kMeansClustering)
    
    