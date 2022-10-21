# LINDA ENCORE
Workflow of LINDA analyses over the ENCORE data

## Organization of the Workflow

This repository shows the steps needed to be performed by the user in order to run the LINDA analysis from the gene and transcript expression data obtained from [ENCORE project](https://www.encodeproject.org/encore-matrix/?type=Experiment&status=released&internal_tags=ENCORE). 

The analysis follows the workflow depicted in the following graph:

<p align="center">
    <img src="https://github.com/enio23/LINDA_ENCORE/blob/main/pipeline.jpeg" alt>
    <em> 1.The Gene and Transcript abundance data for two replicates (of the HepG2 and K562 cell-line) for the Knockdown (*KD*) and Control (*Ctrl*) (GEO:GSE88002) are downloaded and the expected counts are accesed. 2.Genes and Transcripts are filtered based on their expression. 3.With edgeR we perform differential analysis over Gene and Transcript abundance values for the KDvsCtrl comparison. 4. From the differential Gene Expression data we estimate TF activity scores by using the DoRothEA resource as well as a permutation test to estimate the significance of the activity scores (by randomizing Gene-To-TF relations 1000 times). 5. We filter DIGGER resource to only contain interactions between expressed genes. 6. We map the quantified Transcripts to Domain Pfam ID's and for each Domain we estimate an 'effect' score as the average logFC values of Transcript abundances and a 'significance' score after performing a Fissher aggregation of p-values of the mapping transcripts. 7. We give the three inputs to LINDA (TF's/Domains/DIGGER) in order to contextuaize regulatory signalling networks. </em>
</p>


Here are additionally stored the codes needed for the computational benchmarking of the obtained ENCORE networks to the [PerturbSeq study](https://gwps.wi.mit.edu/guide_mde) described in [Replogle, J.M., et.al.](https://www.sciencedirect.com/science/article/pii/S0092867422005979?via%3Dihub).

## How to run the workflow

Please test/run the workflow in the following order:

### Gene_Expression
Here are provided the scripts ([/Gene_Expression/analysis_script.R](https://github.com/enio23/LINDA_ENCORE/tree/main/Gene_Expression)) which have been used to estimate the differentially expressed genes (DGE's) for the KDvsCtrl comparison with [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) as well as TF activity estimation with [viper](https://www.bioconductor.org/packages/release/bioc/html/viper.html) and [DoRothEA](https://www.bioconductor.org/packages/release/bioc/html/viper.html). For more details, please see documentation provided in */Gene_Expression/README.md*.

### Transcript_Expression
Here are provided scripts ([/Transcript_Expression/analysis_script.R](https://github.com/enio23/LINDA_ENCORE/tree/main/Transcript_Expression)) which have been used to assign an Exon ID to all the genomic coordinates in the identified exon skipping events (*/Data/SE.MATS.JunctionCountOnly.txt*).

### LINDA_Hard_Analysis
Here are provided the scripts for the calculation of LINDA networks (splice-aware and splice-unaware) from ENCORE data by using the Hard-Constrained mode (for [HepG2](https://github.com/enio23/LINDA_ENCORE/tree/main/LINDA_Hard_Analysis/HepG2) and [K562](https://github.com/enio23/LINDA_ENCORE/tree/main/LINDA_Hard_Analysis/K562) cell-lines). Additionally here are provided the scripts needed to generate networks in a tabular format which can also be loaded in Cytoscape for visualization (see ```combine_networks.R``` scripts).

### LINDA_Soft_Analysis
Here are provided the scripts for the calculation of LINDA networks (splice-aware and splice-unaware) from ENCORE data by using the Soft-Constrained mode (for [HepG2](https://github.com/enio23/LINDA_ENCORE/tree/main/LINDA_Soft_Analysis/HepG2) and [K562](https://github.com/enio23/LINDA_ENCORE/tree/main/LINDA_Soft_Analysis/K562) cell-lines). Additionally here are provided the scripts needed to generate networks in a tabular format which can also be loaded in Cytoscape for visualization (see ```combine_networks.R``` scripts).

### CARNIVAL
Here are provided the scripts for the calculation of ENCORE networks by using the [CARNIVAL](https://github.com/saezlab/CARNIVAL) method ([Liu A. et.al.](https://www.nature.com/articles/s41540-019-0118-z)). Similar to LINDA, CARNIVAL is able to calculate regulatory networks affecting the downstream regulated TF's, however while CARNIVAL is bale to infer protein activity and interaction signs, it cannot account for splicing effects on protein interactions as LINDA does. CARNIVAL networks were used later for computational benchmarking purposes and comparisons to the LINDA networks.

### PerturbSeq_Benchmarking
In order to evaluate the global relevance of introducing information about splicing effects in protein interaction networks, we have used the pool of LINDA networks generated for the K562  cells from ENCORE and aimed to do a comparison with the other state-of-the-art CARNIVAL method. For such a benchmark, we have additionally relied on  an independent large-scale Perturb-seq experiments using CRISPRi study in which all expressed genes in K562 chronic myeloid leukemia cells were targeted (n=9876 targets repressed) ([Replogle et al. 2022](https://www.sciencedirect.com/science/article/pii/S0092867422005979?via%3Dihub)). For details about the benchmarking steps followed, please refer to the LINDA [manuscript](to_be_added).