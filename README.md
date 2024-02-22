## 1 conda 安装

### 下载和安装

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
sh Miniconda3-py39_4.12.0-Linux-x86_64.sh
source ~/.bashrc
```
## 2 conda 环境配置

```
conda create -n (环境名) python=3.9
conda activate (环境名) 
source activate (环境名)
#添加镜像源
conda config –add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/
conda config –add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/
conda config –add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free/
conda config –add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main/
##查看镜像源
conda config –show-sources
```

## 3 安装R（v4.2）以及R包
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


## 4 Cell Ranger安装

```
wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1675964753&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWl
jcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NzU5NjQ3NTN9fX1dfQ__&Signature=IFr7ONDqEkZRR6QpU~A6719a9Mc2SD2
tI1z6RrGldFFTCiY6Z7VR0x0Gr90jtvTUmYTJ2S0NyuK6SVmdeIZUCcbjz9elG1ImGx7AprTCRD3m~0se-xha2lFr87bEsbAa-7uoyW14wXRlj17b0oG9WomNvVSNNKJSzSSfCkqX3Ev9B82b~DMD-7-Hlb8lAsorv18R8y41T4UihIRdY-LE-I5Gk3fTod
mBUjvSEuI3VEalsrVsrN5AdBDpwiCPqSiExODVM0RIsUDV158ceAYFiu5Y9wgbwQVOFMZGYI0d-6tO1VPo4RwwWl0X7c2q21im6BNSQrhQzoDv5cj5COesmQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -zxvf cellranger-7.1.0.tar.gz
echo "PATH=安装路径/cellranger/7.1.0/bin:$PATH" >>~/.bashrc
```
## 5 seqkit安装
```
conda install -c bioconda seqkit
```
