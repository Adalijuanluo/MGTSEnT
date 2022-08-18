# MGTSEnT
Custom scripts to report the key information based on the output of Multilevel Genome Typing (MGT)  for global _Salmonella_ Enteritidis https://mgtdb.unsw.edu.au/enteritidis/

The key information includes:
1. Population structure of typed _Salmonella_ Enteritidis isolates.
2. Isolates belonging to the multidrug resistance associated MGT STs (Sequence Types). 
3. Summarisation of closely related clusters based on the highest reolution typing level MGT9 GC (Genomic Cluster).
4. Generating Microreact input for visualisation of a cluster of interest (i.e. MGT-GC based dentrogram, metadata table). 
## Author
Ada Lijuan Luo, Ruiting Lan Laboratory, University of New South Wales
## Input and Output
### Input
1. The whole MGT typing dataset from the MGTdb https://mgtdb.unsw.edu.au/enteritidis/ (.csv).
2. Accession number list for the newly sequenced strains or strains of interest (.txt).
### Output
1. CSV file about the population structure (clades & lineages) of the strains, and multidrug resistance associated sequence types (STs). 
2. CSV file of closely related clusters at different resolution levels (by different allele difference cutoffs). 
