# Targeted co-expression networks precisely delineate the pathways associated with the APP protein within Alzheimer's brains

## Motivation
Gene co-expression networks (GCNs) on transcriptomes are too coarse grain models to specifically investigate associations of gene expression with phenotypes of interest. Targeted GCNs (TGCNs) are a new modelling approach we designed that creates more specific networks, focused on a phenotype, APP in this case. A TGCN tries to explain as much variation of the phenotype as possible within the sample while describing, functionally, the genes involved. We created [APP-TGCNs on the ROSMAP](https://htmlpreview.github.io/?https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_ROSMAP.html) gene expression profiles and replicated results on the [MSBB cohort](https://htmlpreview.github.io/?https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_MSBB.html).

## Results
The [ROSMAP APP-TGCN](https://htmlpreview.github.io/?https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_ROSMAP.html) has 25 modules 880 genes size. The hubs explain most of the gene expression variation in APP (R2 0.96). Within the hubs we find the ATP1B1 gene (leading a module enriched for dopaminergic neurons and GO terms like synaptic vesicle cycle). ATP1B1 replicates as hub in the [MSBB APP-TGCN](https://htmlpreview.github.io/?https://github.com/aliciagp/TGCN/blob/master/notebooks/validation.html). ABCD3 leads genes enriched for astrocytes and activity located mainly in the cell plasma membrane. FRZB is hub in a module enriched for endothelial and mural cells. NFASC module is enriched for oligodendrocytes and neuron development related terms, and replicates in the MSBB TGCN. In association to AD as a phenotype, the APP-TGCN shows APP involvement through the HNRNPA2B1 module (enriched with mRNA metabolism and splicing annotations), and ATP1B1 and NFASC module.

We show that APP expression in brain plays cell-specific key roles in several biological processes: neuron development through oligodendrocytes and in regulation of the synaptic vesicle cycle in dopaminergic neurons, which we found associated with AD pathobiology.


## How to use the TGCN R package

You can install the development version of the TGCN R package like this:

``` r
remotes::install_github("aliciagp/TGCN")
```

To create a targeted gene co-expression network, you just need to load the library and run this chunk of code:

``` r
library(TGCN)

testAllCutoffs(exprData=input,        # gene expression matrix
               target=toPredict,      # target to predict
               covs=covs,             # donors metadata
               train.split=0.7,       # proportion of samples use to train the models
               nfolds=5,              # number of folds for the cross-validation
               t=10,                  # number of LASSO runs
               path=tgcn_path,        # path where results will be saved
               targetName="APP",      # name of the target
               tissueName="ROSMAP",   # name of the tissue
               seed=1234,             # seed to ensure reproducibility
               cutoffs=10:1,          # cutoffs to be tested (number of times a gene was selected by LASSO)
               n=100,                 # the maximum module size
               m=10,                  # the number of genes added in each iteration if approach="enrichment" 
               s=10,                  # the minimum module size 
               minCor=0.3,            # the minimum correlation of the genes added to the modules
               maxTol=3,              # the maximum number of tries to get a better enrichment if approach=enrichment
               save=T,                # if save=T, results will be saved as separate files
               overwrite=T,           # if overwrite=T, modules annotation will be overwritten
               approach="enrichment", # the approach selected to complete the seed modules
               report=T)              # if report=T, an automated report will be created

```

Results will be saved into three different folders: hubGenes, Net and results. Results folders includes an html report of the TGCN created. Here you can find the reports for the creation and replication of the APP-TGCNs:

- [Notebook for the creation and results of the ROSMAP based APP-TGCN](https://htmlpreview.github.io/?https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_ROSMAP.html)

- [Notebook for the creation and results of the MSBB based APP-TGCN used for replication of results from the ROSMAP APP-TGCN](https://htmlpreview.github.io/?https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_MSBB.html)

- [Notebook for the replication of ROSMAP APP-TGCN on the MSBB APP-TGCN](https://htmlpreview.github.io/?https://github.com/aliciagp/TGCN/blob/master/notebooks/validation.html)


## Credits

The authors of the package are Alicia Gómez-Pascual, PhD student, Funded by Fundación Séneca and Juan A. Botía, Professor from the University of Murcia, Spain. This method was developed with the collaboration of Laura Ibañez from Washington University in Saint Louis, USA. The development of the software and some aspects of the algorithm design have been possible thanks to the previous work of Guillermo Rocamora at GOS Institute of Child Health, University College London, UK. We also thank Mina Ryten and Sonia García at GOS Institute of Child Health for their previous works.

 

