library("optparse")

option_list <- list(
	make_option(c("-i", "--indir"), help="indir"),
	make_option(c("-o", "--outdir"), help="outdir", default="."),
	make_option(c("--res"), help="resolution, default %default",type="double",default=0.5),
	make_option(c("--MinCell"), help="MinCell, default %default",type="double",default=5),
	make_option(c("--MinFeatures"), help="MinFeatures, default %default",type="double",default=100),
	make_option(c("--point_size"), help="point_size, default %default",type="double",default=1)
)
opt <- parse_args(OptionParser(option_list=option_list))
#######################################
# 首先将需要用到的包load进工作环境
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# 记录运行时间
StartTime=proc.time()


##绘图函数
#输出ggplot对象
SavePlot <- function(od, filename, data, width = 8, height = 8){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height)
}



# 参数初始化
#选择不同的resolution值可以获得不同的cluster数目，值越大cluster数目越多，默认值是0.5
res=opt$res
FilePath=opt$indir
he_file=opt$he_png
#过滤参数：细胞数5，基因数100
MinCell=opt$MinCell
MinFeatures=opt$MinFeatures

#输出目录
outdir=opt$outdir
if(!dir.exists(outdir)){dir.create(path=outdir)}


# 数据导入，给定输入文件路径，路径下需包括matrix.mtx,genes.tsv(或 features.tsv),和barcodes.tsv文件
expr=Read10X(FilePath,cell.column = 1)

#创建Seurat对象，数据集中测到的少于200个基因的位点（min.features = 200）和少于5个位点覆盖的基因（min.cells = 3）被过滤掉
object=CreateSeuratObject(counts = expr,assay = "RNA",min.cells=MinCell,min.features=MinFeatures)
#标准化，鉴定出表达量高变的基因
#object <- SCTransform(object, assay = "Spatial", verbose = FALSE,return.only.var.genes = FALSE)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
object <- ScaleData(object)

object <- RunPCA(object=object,pc.genes = VariableFeatures(object))
object <- FindNeighbors(object=object, reduction = "pca",dims = 1:30,verbose = F)
object <- FindClusters(object=object, verbose = F,resolution = res)
object <- RunUMAP(object=object, reduction = "pca",dims = 1:30)
object <- RunTSNE(object=object, reduction = "pca",dims = 1:30,check_duplicates = FALSE)

#整理聚类结果
Cluster=object@meta.data %>%
  select(.,seurat_clusters) %>%
  rownames_to_column(.,var = "Barcode") %>%
  arrange(.,seurat_clusters)

#输出聚类结果
write.csv(Cluster,file.path(outdir,"cluster.csv"),quote = F,row.names = F)

# 计录聚类时间
EndTime=proc.time()
RunTime <- EndTime - StartTime
print(paste0('ProcessingTime: ',RunTime[3][[1]],'seconds'))


#找到各个类别的marker gene
object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25)

#提取每个类别的marker gene，n=10表示提取每个类别中top10的marker gene
top10 <- object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#将每个类别中的marker gene输出
write.csv(top10,file.path(outdir,"marker_gene.csv"),quote = F, row.names = F)

Cluster=object@meta.data %>%
  select(.,c(seurat_clusters)) %>%
  rownames_to_column(.,var = "Barcode") %>%
  arrange(.,seurat_clusters)

col = c("#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#70e014","#787878","#DB4C6C","#0430e0","#554236","#AF5F3C","#ff7700","#e00417","#DAB370","#fcfc05","#268785","#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f","#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1", "#51f59b")
UMAP_df <- object[['umap']]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = 'Barcode') %>%
  left_join(Cluster, by = 'Barcode') %>%
  dplyr::rename(cluster = seurat_clusters) %>%
  select(Barcode, cluster, everything())
UMAP_df$cluster <- as.factor(UMAP_df$cluster)

###Output umap
write.table(x = UMAP_df,file = paste(outdir,'umap_df.xls',sep="/"),row.names = F,
            col.names = T,quote = F,append = F,sep = "\t")

umap_pic <- ggplot(UMAP_df, aes(x = UMAP_1, y = UMAP_2))+
      geom_point(aes(color = cluster),size=opt$point_size, shape=16)+
      scale_color_manual(values = col)+
      labs(colour = "cluster")+
      theme_bw(base_size = 15)+
      theme(panel.grid = element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=3.5)))

ggsave(filename = paste(outdir,'umap.pdf',sep="/"),width = 12/3.8*2,height = 5,dpi = 600)
ggsave(filename = paste(outdir,'umap.png',sep="/"),width = 12/3.8*2,height = 5,dpi = 600)



TSNE_df <- object[['tsne']]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = 'Barcode') %>%
  left_join(Cluster, by = 'Barcode') %>%
  dplyr::rename(cluster = seurat_clusters) %>%
  select(Barcode, cluster, everything())
TSNE_df$cluster <- as.factor(TSNE_df$cluster)

###Output tsne
write.table(x = TSNE_df,file = paste(outdir,'tsne_df.xls',sep="/"),row.names = F,
            col.names = T,quote = F,append = F,sep = "\t")

tsne_pic <- ggplot(TSNE_df, aes(x = tSNE_1, y = tSNE_2))+
      geom_point(aes(color = cluster),size=opt$point_size, shape=16)+
      scale_color_manual(values = col)+
      labs(colour = "cluster")+
      theme_bw(base_size = 15)+
      theme(panel.grid = element_blank())+
      guides(colour = guide_legend(override.aes = list(size=3.5)))
ggsave(filename = paste(outdir,'tsne.pdf',sep="/"),width = 12/3.8*2,height = 5,dpi = 600)
ggsave(filename = paste(outdir,'tsne.png',sep="/"),width = 12/3.8*2,height = 5,dpi = 600)



