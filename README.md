# CBASS_Nursery

Project: Short-term thermal stress assays using CBASS to assess thermal threshold of donor versus year-long in situ nursery-propagated corals.

Data and code for RNASeq analysis and determining ED50 of FvFm over temperature.

### Workflow

1.	Reference-based RNA-Seq pipeline for paired end reads `RNASeq_RACH.sh`  
    Using Trimmomatic to trim reads, STAR to align reads to reference genome  
  	& StringTie to assembly transcripts & create gene count matrix.
  	
2.  Gene-level expression overview `Differential_gene_expression.R`  
    Consisting of Principal Components Analysis (PCA), differential gene expression using DESeq2  
    & venn diagram comparsions.  
    
3.  Inferring functional gene annotation `Functional_annotation.R`  
    Using topGO for Gene Ontology enrichment analysis  
    & EggNOG-based functional inferences to annotate differentially expressed genes of interest.

4.  Dose-response curve modelling `DRC_modelling.R`  
    Using DRC R package to model effective dose 50 of FvFm over temperature.

Publication: 
Alderdice R, Voolstra CR, Lendo CIN, Boote C, Suggett DJ, Edmondson J, et al. Loss of coral thermotolerance following year-long in situ nursery propagation with a consecutively high summer heat-load. Coral Reefs. 2024. doi:10.1007/s00338-024-02505-9 [link](https://link.springer.com/article/10.1007/s00338-024-02505-9)

Raw fastq files available at: [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1011835](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1011835)  
Additional large input files for shell pipeline available at: [https://doi.org/10.5281/zenodo.8431668](https://doi.org/10.5281/zenodo.8431668)

R version 4.3.0 used in scripts
