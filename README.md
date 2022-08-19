# MGTSEnT
Custom pipeline to report the key information based on the output of Multilevel Genome Typing (MGT)  for global _Salmonella_ Enteritidis https://mgtdb.unsw.edu.au/enteritidis/

The key information includes:
1. Population structure of typed _Salmonella_ Enteritidis isolates.
2. Isolates that belong to the multidrug resistance associated MGT STs (Sequence Types). 
3. Summarisation of closely related clusters based on the highest resolution typing level MGT9 GC (Genomic Cluster).
4. Generating Microreact input for visualisation of a cluster of interest (i.e. MGT-GC based dendrogram, metadata table). 

## Author
Ada Lijuan Luo, Ruiting Lan Laboratory, University of New South Wales
## Input and Output
### Input 
1. The whole MGT typing dataset (.csv) from the MGTdb https://mgtdb.unsw.edu.au/enteritidis/. 
2. Accession number list (.txt) for the newly sequenced strains or strains of interest.

   #### The following files are required for MGTSEnT_Microreactinput.py

   * List of closely related clusters (.txt) for microreact visualisation (e.g. GC10_1813). 
   * Allele profile downloaded from https://mgtdb.unsw.edu.au/enteritidis/.
   
   #### Input file provided by Github repository in the updating_input directory
   * Flag file (.csv) that includes the population structure and multidrug resistance associated MGT types. This file will be regularly updated. Please use the latest version. 
   * Input 1 of the whole MGT typing dataset from the MGTdb is also provided, however, this file may not include the newly sequenced or private genomes. 
   
### Output
1. CSV file about the population structure (clades & lineages) of the strains, and multidrug resistance associated sequence types (STs). 
2. CSV files of closely related clusters at different resolution levels of clustering (using different pair-wise allele difference cutoffs, 0, 1, 2, 5 and 10). The following information are provided: distribution of country and year, classification of 4-week sliding window (if include >= 2 isolates collected within 4 weeks), National/International, Within-State/Inter-State. 
3. Microreact input files for investigating a cluster of interest. 
    * Metadata (.txt) file including geographic and collection time information, and MGT9-based genomic clusters at different clustering threshold from 0 to 10 using the single linkage clustering algorithm. 
    * Dendrogram (.txt) for the hierarchical single linkage clusters. 

## Installation
1. Update or Install Miniconda and add the bioconda channel:
   - conda update -n base -c defaults conda
   
     or 
   
   - Install miniconda3  -> https://docs.conda.io/en/latest/miniconda.html
   - Add bioconda -> http://www.ddocent.com//bioconda/
2. Create miniconda3 environment:
````
conda create -n mgtsent
conda activate mgtsent
pip install geopy
````
3. Download this repository:
````
git clone https://github.com/Adalijuanluo/MGTSEnT.git
cd MGTSEnT
````
## Usage
* Output 1
````
python MGTSEnT.py -h
python MGTSEnT.py -m path/mgtdb.csv -i path/accession_number_list.txt -f updating_input/latest_flag.csv -o outputpath/prefix_
````
* Output 2
````
python MGTSEnT_GC_summary.py -h
python MGTSEnT_GC_summary.py -m path/mgtdb.csv -i path/accession_number_list.txt -f updating_input/latest_flag.csv -o outputpath/prefix_
````
* Output 3
````
python MGTSEnT_Microreactinput.py -h
python MGTSEnT_Microreactinput.py -m path/mgtdb.csv -c path/list_of_clusters.txt -outputpath/prefix_
````


