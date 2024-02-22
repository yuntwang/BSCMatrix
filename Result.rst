1 数据准备
1.1.测序数据，双端FASTQ数据
1.2.参考基因组数据，基因组序列文件和gtf文件（配置文件使用的注释文件需要包含gene和exon）

2 配置文件填写
### Globle parameters
## data  测序数据（压缩非压缩均可）
FQ1      /path/to/read_1.fastq
FQ2      /path/to/read_2.fastq

## Ref genome  参考基因组和gtf文件
FA       /path/to/ref/file
GTF      /path/to/gtf/file

## out put    输出路径以及输出结果前缀
OUTDIR   /path/to/result/dir/
PREFIX    outfile-prefix

## other parameters   其它配置参数（期望细胞数、线程数目、内存上限100Gb，Nobam 设置为0代表不输出bam，设置其它字符代表输出）
EC          3000
Threads      8
RAM         100
Nobam       1
##ENV  ###R（可选，如果不提供，则放到环境变量中）
Rscript   /path/to/R/bin/

3 流程运行

3.1.流程说明
流程分为4个步骤，功能如下：
3.1.1步骤1，运行fastq2BcUmiSC_v1.1，识别fastq数据中的barcode、umi
3.1.2步骤2，调用Cell Ranger程序，获取基因表达矩阵
3.1.3步骤3，运行QC，进行umap和tsne分析
3.1.4步骤4，运行WebReport，得到网页版报告
3.2.流程参数
-c  config.txt 配置文件
-s  步骤选择，0为运行1-4所有步骤，也可选择个别步骤单独运行，多个步骤中间使用”,”分割。
3.3.参考命令
./BSCMatrix -c config.txt -s 0
./BSCMatrix -c config.txt -s 1,2,3,4
./BSCMatrix -c config.txt -s 1,2

4 结果说明
01.fastq2BcUmiSC           ##步骤1运行结果目录
├── prefix.bc_stat        ##不同barcode统计
├── prefix.filter          ##barcode 过滤信息
├── prefix.qual.stat       ##Reads 统计
├── prefix.stat           ##barcode检测统计
└── prefix.umi           ##Reads barcode umi信息

02.cellranger/               ##步骤2 运行结果目录
├── bmk_10x_barcode.xls    ## barcode对照表
├── data/                 ##数据存放目录
├── INDEX/                ##参考基因组索引文件夹
├── Log.out                ##索引构建日志文件
├── R1.fq.gz               ## R1端数据
├── R2.fq.gz               ## R2 端数据
├── ReadID.xls             ##ReadID 信息
└── prefix/                 ##调用Cell Ranger 软件分析结果目录

03.cluster/                 ## 步骤3 运行结果目录
├── cluster.csv             ##细胞聚类结果
├── marker_gene.csv        ##marker gene
├── tsne_df.xls             ##tsne聚类结果
├── tsne.pdf               ##tsne pdf格式图片
├── umap_df.xls            ##umap 聚类结果
└── umap.pdf              ##umap pdf格式图片

04.WebReport/             ##步骤4 运行结果目录
├── 10x                  ##调用Cell Ranger程序执行的结果
└── bmk                  ##网页报告
