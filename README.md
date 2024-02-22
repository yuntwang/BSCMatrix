# Dependency Environment Configuration

##  Conda Installation

### download and installation

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
sh Miniconda3-py39_4.12.0-Linux-x86_64.sh
source ~/.bashrc
```
## Conda Environment Configuration

```
conda create -n your_env_name python=3.9
conda activate your_env_name
source activate your_env_name

conda config –add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/
conda config –add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/
conda config –add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free/
conda config –add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main/

conda config –show-sources
```

## Installing R (v4.2) and R Packages

```
conda install -c conda-forge r-base=4.2
conda install -c conda-forge r-seurat=4.3.0
conda install -c conda-forge r-dplyr=1.1.0
conda install -c conda-forge r-tibble=3.1.8
conda install -c conda-forge r-kableextra=1.3.4
conda install -c conda-forge r-knitr=1.42
conda install -c conda-forge r-rmarkdown=2.20
conda install -c conda-forge r-optparse=1.7.3
conda install -c conda-forge r-scales=1.3.0
install.packages('DT')
```

## Installing CellRanger

```
wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1675964753&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWl
jcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NzU5NjQ3NTN9fX1dfQ__&Signature=IFr7ONDqEkZRR6QpU~A6719a9Mc2SD2
tI1z6RrGldFFTCiY6Z7VR0x0Gr90jtvTUmYTJ2S0NyuK6SVmdeIZUCcbjz9elG1ImGx7AprTCRD3m~0se-xha2lFr87bEsbAa-7uoyW14wXRlj17b0oG9WomNvVSNNKJSzSSfCkqX3Ev9B82b~DMD-7-Hlb8lAsorv18R8y41T4UihIRdY-LE-I5Gk3fTod
mBUjvSEuI3VEalsrVsrN5AdBDpwiCPqSiExODVM0RIsUDV158ceAYFiu5Y9wgbwQVOFMZGYI0d-6tO1VPo4RwwWl0X7c2q21im6BNSQrhQzoDv5cj5COesmQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -zxvf cellranger-7.1.0.tar.gz
echo "PATH=PATH/cellranger/7.1.0/bin:$PATH" >>~/.bashrc
```
##  Installing seqkit
```
conda install -c bioconda seqkit
```

# Operating Instructions for BSCMatrix

## Filling in the Configuration File
```
FQ1      /path/to/read_1.fq.gz
FQ2      /path/to/read_2.fq.gz
FA       /path/to/ref/file
GTF      /path/to/gtf/file
OUTDIR   /path/to/result/dir/
PREFIX    outfile-prefix
## other parameters   Other configuration parameters (expected number of cells, number of threads, maximum memory limit of 100Gb, "Nobam" set to 0 represents no output of bam, and setting other characters represents output)
EC          3000
Threads      8
RAM         100
Nobam       1
##ENV  Optional. If not provided, it will be placed in the environment variable
Rscript   /path/to/R/bin/
```

## BSCMatrix Execution

The process is divided into four steps with the following functions:

Step 1: Run fastq2BcUmiSC_v1.1 to identify barcodes and UMIs from fastq data.

Step 2: Call the Cell Ranger program to obtain the gene expression matrix.

Step 3: Run QC to perform UMAP and t-SNE analysis.

Step 4: Run WebReport to generate a web-based report.

## Reference Commands

./BSCMatrix -c config.txt -s 0

./BSCMatrix -c config.txt -s 1,2,3,4

./BSCMatrix -c config.txt -s 1,2

