# SPADA
Code for SPADA lab genomics research on bulk and single-cell RNA-sequencing data.

The RelatednessTesting folder contains Python code which implements
-a joint bulk and scRNA-seq data simulation based on the URSM (https://projecteuclid.org/euclid.aoas/1520564486) and ZIFA (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0805-z) models
-a hypothesis test for whether a bulk and single-cell dataset are joint (come from similar tissue)
-a new algorithm called cluster heterogeneity for measuring the similarity of clustered datasets.

The Preprocessing folder contains C++ code for handling datasets stored in csv and tsv files.
