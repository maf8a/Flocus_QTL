# Scripts

***
***

### `DESeq2` 
***
R script and associated samples file (input) used for differential gene expression analysis between QTL parent clones. In manuscript, see S5_File for raw output data and see Fig 4, and S6_Figure and S7_Figure for plots. Raw RNA-seq reads used for analysis were generated in a [previous study](https://www.nature.com/articles/sdata201630) and are availabe on [GenBank](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA284518). 
### `minimap2`
***
bash script used to map transcripts from each QTL parent clone to the respective genome. Directory also includes fasta files (output from script) containing those transcripts that mapped to the F-locus region. Mapping data (bam files) generated here were used to curate structural gene annotations. For additional details, see S1_Methods and S6_File in manuscript.

### `prank`
***
R script and associated input files for aligning orthologous coding sequences (cds) for genes that are shared between xF and iF (i.e. F-locus region from QTL parent clones Xinb3 and Iinb1, respectively). Subdirectories `xinb3` and `iinb1` contain fasta files of the cds sequences from the respective parent clone. Note that the file names in each subdirectory are identical, which is required for the script to run properly. The alignments created here were used for calculating piN/piS values (see Table 2 in manuscript). 

### `RepeatModeler`
***
bash scripts for running RepeatModeler and RepeatMasker on each QTL parent genome to identify repeat elements in the F-locus region. For output data, see S4_File in manuscript.

### `Rqtl`
***
R scripts for running QTL analysis. For input data, see S1_File, sheet "QTL SNP map" in manuscript.