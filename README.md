# Targeted co-expression networks precisely delineate the pathways associated with the APP protein within Alzheimer's brains

Gene co-expression networks (GCNs) on transcriptomes are too coarse grain models to specifically investigate associations of gene expression with phenotypes of interest. Targeted GCNs (TGCNs) are a new modelling approach we designed that creates more specific networks, focused on a phenotype, APP in this case. A TGCN tries to explain as much variation of the phenotype as possible within the sample while describing, functionally, the genes involved. We created [APP-TGCNs on the ROSMAP](https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_ROSMAP.html) gene expression profiles and replicated results on the [MSBB cohort](https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_MSBB.html).

The [ROSMAP APP-TGCN](https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_ROSMAP.html) has 25 modules 880 genes size. The hubs explain most of the gene expression variation in APP (R2 0.96). Within the hubs we find the ATP1B1 gene (leading a module enriched for dopaminergic neurons and GO terms like synaptic vesicle cycle). ATP1B1 replicates as hub in the [MSBB APP-TGCN](https://github.com/aliciagp/TGCN/blob/master/notebooks/validation.html). ABCD3 leads genes enriched for astrocytes and activity located mainly in the cell plasma membrane. FRZB is hub in a module enriched for endothelial and mural cells. NFASC module is enriched for oligodendrocytes and neuron development related terms, and replicates in the MSBB TGCN. In association to AD as a phenotype, the APP-TGCN shows APP involvement through the HNRNPA2B1 module (enriched with mRNA metabolism and splicing annotations), and ATP1B1 and NFASC module.

We show that APP expression in brain plays cell-specific key roles in several biological processes: neuron development through oligodendrocytes and in regulation of the synaptic vesicle cycle in dopaminergic neurons, which we found associated with AD pathobiology.


Here you can find the notebooks for the creation and replication of the APP-TGCNs:

- [Notebook for the creation and results of the ROSMAP based APP-TGCN](https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_ROSMAP.html)

- [Notebook for the creation and results of the MSBB based APP-TGCN used for replication of results from the ROSMAP APP-TGCN](https://github.com/aliciagp/TGCN/blob/master/notebooks/template_APP_MSBB.html)

- [Notebook for the replication of ROSMAP APP-TGCN on the MSBB APP-TGCN](https://github.com/aliciagp/TGCN/blob/master/notebooks/validation.html)


