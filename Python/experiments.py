#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 21:57:40 2020

@author: fcseidl

Generate data and visualizations for tech report.
"""

import numpy as np
from sklearn.preprocessing import normalize
from sklearn.metrics import normalized_mutual_info_score
import matplotlib.pyplot as plt

from BestMap import BestMap
from SSC import sparseSubspaceClustering
from KSubspaces import KSubspaces

import preprocessing
import simulations as sims
import RelatednessTesting as rt
import SPADAutil as util


# data--note that paths are local and will not work on other machines
print("Loading Chen data matrix...")
chen = preprocessing.csvToMatrix(
    "/Users/fcseidl/Documents/SPADA/RNAseq/mouse_brain/ssf_chen.tsv", 
    delim='\t')

print("Loading Campbell data matrix...")
campbell = preprocessing.csvToMatrix(
    "/Users/fcseidl/Documents/SPADA/RNAseq/mouse_brain/ssf_campbell.tsv", 
    delim='\t')

print("Loading Segerstolpe data matrix...")
seger = np.genfromtxt(
    "/Users/fcseidl/Documents/SPADA/RNAseq/human_pancreas/ssf_segerstolpe.tsv", 
    delimiter='\t')
seger = seger[1:, 2:]

print("Loading Xin data matrix...")
xin = np.genfromtxt(
    "/Users/fcseidl/Documents/SPADA/RNAseq/human_pancreas/ssf_xin.tsv", 
    delimiter='\t')
xin = xin[1:, 1:]

# load cell labels
print("Loading Segerstolpe labels...")
s_lbl = np.genfromtxt(
    "/Users/fcseidl/Documents/SPADA/RNAseq/human_pancreas/Segerstolpe/labels.tsv", 
    delimiter='\t', dtype=str)
s_lbl = s_lbl[1:, -1].astype(str)

print("Loading Xin labels...")
x_lbl = np.genfromtxt(
    "/Users/fcseidl/Documents/SPADA/RNAseq/human_pancreas/Xin/labels.tsv", 
    delimiter='\t', dtype=str)
x_lbl = x_lbl[1:, -1].astype(str)


# preprocessing
print("Chen and Campbell preprocessing:")
print("Removing sparse genes from chen and campbell...")
def sparse(Yn):
    return util.dropoutRate(Yn) > 0.8
chen, campbell = preprocessing.removeRowsPred(chen, campbell, sparse)
campbell, chen = preprocessing.removeRowsPred(campbell, chen, sparse)

print("Randomly choosing subset of data...")
L = 500
ch_idx = np.random.permutation(np.arange(chen.shape[1]))[:L]
ca_idx = np.random.permutation(np.arange(campbell.shape[1]))[:L]
chen = chen[:, ch_idx]
campbell = campbell[:, ca_idx]

print("Applying cosine normalization to samples...")
chen = normalize(chen)
campbell = normalize(campbell)

print("Segerstolpe and Xin preprocessing:")

s_types = ["alpha cell", "beta cell", "delta cell", "gamma cell"]
x_types = ["alpha", "beta", "delta", "PP"]

print("Removing additional cell types...")
remove = []
for i in range(seger.shape[1]):
    if s_lbl[i] not in s_types:
        remove.append(i)
seger = np.delete(seger, remove, axis=1)
s_lbl = np.delete(s_lbl, remove)

remove = []
for i in range(xin.shape[1]):
    if x_lbl[i] not in x_types:
        remove.append(i)
xin = np.delete(xin, remove, axis=1)
x_lbl = np.delete(x_lbl, remove)

print("Removing sparse genes from both datasets...")
def sparse(Yn):
    return util.dropoutRate(Yn) > 0.7
seger, xin = preprocessing.removeRowsPred(seger, xin, sparse)
xin, seger = preprocessing.removeRowsPred(xin, seger, sparse)

print("Randomly choosing subset of data...")
L = 500
s_idx = np.random.permutation(np.arange(seger.shape[1]))[:L]
seger = seger[:, s_idx]
s_lbl = s_lbl[s_idx]
x_idx = np.random.permutation(np.arange(xin.shape[1]))[:L]
xin = xin[:, x_idx]
x_lbl = x_lbl[x_idx]

print("Applying cosine normalization to samples...")
seger = normalize(seger)
xin = normalize(xin)



def mouse_hypothalamus(clusterer):
    n_clusters = [2, 5, 10] # various numbers of clusters to use
    
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
        
        index = np.arange(k)
        ticks = np.arange(k)
        width = 0.4
        y_offset = np.zeros(k)
        
        #plt.bar(index + width, joint_distribution, width, label=("joint"))
        
        for c in range(k):
            plt.bar(index, chen_distribution[c], width, bottom=y_offset, 
                    label=("chen %d" % c))
            y_offset += chen_distribution[c]
            plt.bar(index, campbell_distribution[c], width, bottom=y_offset,
                    label=("campbell %d" % c))
            y_offset += campbell_distribution[c]
            
        plt.ylabel('number of points')
        plt.xlabel('joint cluster')
        plt.title('Cluster Composition')
        plt.xticks(index + width / 2, ticks)
        plt.legend(loc='best')
        plt.show()
        

def human_pancreatic_islets(clusterer):
    n_clusters = [4, 10] # various numbers of clusters to use
    
    for k in n_clusters:
        print("Clustering both datasets together into", k, "clusters...")
        joint = np.concatenate((seger, xin), axis=1)
        lbl = clusterer(joint.T, k)
        
        print("NMI between Segerstolpe true and predicted labels =",
              normalized_mutual_info_score(s_lbl, lbl[:L]))
        
        print("NMI between Xin true and predicted labels =",
              normalized_mutual_info_score(x_lbl, lbl[L:])) 
        
        print("Plotting cluster composition bar chart...")
        dist = np.zeros((k, 4))
        for i in range(L):
            # Segerstolpe data
            predicted = lbl[i]
            truth = s_lbl[i]
            dist[predicted, s_types.index(truth)] += 1
            # Xin data
            predicted = lbl[L + i]
            truth = x_lbl[i]
            dist[predicted, x_types.index(truth)] += 1
        
        index = np.arange(k)
        ticks = np.arange(k)
        width = 0.4
        y_offset = np.zeros(k)
        
        for j in range(4):
            plt.bar(index, dist[:, j], width, bottom=y_offset, 
                    label=x_types[j])
            y_offset += dist[:, j]
            
        plt.ylabel('number of points')
        plt.xlabel('cluster')
        plt.title('Cluster Composition')
        plt.xticks(index + width / 2, ticks)
        plt.legend(loc='best')
        plt.show()
        
        
# algorithms
kmeans = util.kMeansClustering
lSSC = sparseSubspaceClustering
aSSC = lambda X, k : sparseSubspaceClustering(X, k, affine=True)
ksubspaces = KSubspaces


# collect results
print("\n\n--------------------k-means on mice data-------------------------")
mouse_hypothalamus(kmeans)

print("\n\n------------------linear SSC on mice data------------------------")
mouse_hypothalamus(lSSC)

print("\n\n------------------affine SSC on mice data------------------------")
mouse_hypothalamus(aSSC)

print("\n\n------------------k-subspaces on mice data-----------------------")
mouse_hypothalamus(ksubspaces)

print("\n\n------------------k-means on pancreas data-----------------------")
human_pancreatic_islets(kmeans)

print("\n\n----------------linear SSC on pancreas data----------------------")
human_pancreatic_islets(lSSC)

print("\n\n----------------affine SSC on pancreas data----------------------")
human_pancreatic_islets(aSSC)

print("\n\n----------------k-subspaces on pancreas data---------------------")
human_pancreatic_islets(ksubspaces)
    
    