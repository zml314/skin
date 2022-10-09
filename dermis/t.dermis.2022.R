####t.dermis####
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(scater)
library(SingleR)
#memory.limit(102400)
library(monocle)
library(tibble)
library(harmony)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(DOSE)
setwd('/home/zhaolab/R_analasis/dermis.harmony/dermis.res1/t.dermis/')
setwd('E:/R_analysis/dermis.2022.1/')
dir.create('./t')
setwd('./t')
####2.创建存储图形的文件夹####
sam.name <- "multi"
if(!dir.exists(sam.name)){
              dir.create(sam.name)
}
load('../T.dermis.Rdata')
merged=T
rm(T)


#merged@meta.data=transform(merged@meta.data,sample=ifelse(orig.ident%in%c('HC_1','HC_2','HC_3'),'HC',merged$sample))
merged=subset(merged,seurat_clusters!='5')
merged <- NormalizeData(merged, normalization.method = "LogNormalize",scale.factor = 10000)


##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)


#展示标准化之后的整体表达水平
top10 <- head(x = VariableFeatures(merged), 10)
plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()

##而对所有基因进行标准化的方法如下：
#memory.limit(102400)
all.genes <- rownames(merged)
memory.limit(102400)
merged <- ScaleData(merged, features = all.genes)#需要足够大的内存
merged <- ScaleData(merged, vars.to.regress = "percent.mt")
####7.线性降维（PCA）并存储####
#默认用高变基因集,但也可通过features参数自己指定；
merged <- RunPCA(merged, features = VariableFeatures(object = merged))


#PCA结果展示-1
pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 20,height = 20)
VizDimLoadings(merged, dims = 1:2, reduction = "pca")
dev.off()

#PCA结果展示-2
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(merged, reduction = "pca",group.by = 'orig.ident')
dev.off()

#PCA结果展示-3
pdf(paste0("./",sam.name,"/PCA-DimHeatmap.pdf"),width = 20,height = 20)
DimHeatmap(merged, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(merged,ndims = 40)
dev.off()
dim.use=1:20
#TSNE算法
library(harmony)

merged=RunHarmony(merged,group.by.vars ='orig.ident')
library(clustree)
merged.1=FindClusters(merged,resolution = c(seq(0.4,2.8,0.2)))
clustree(merged.1@meta.data,prefix='RNA_snn_res.')
merged<- FindNeighbors(merged, dims = dim.use,reduction = 'harmony',nn.method = "rann")
merged <- FindClusters(merged, resolution =2.2)
merged <- RunTSNE(merged, dims = dim.use,reduction = 'harmony')

pdf(paste0("./",sam.name,"/tsne_Plot.harmony_res2.2_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = T)
dev.off()
save(merged,file = 't.res0.3.20220117.Rdata')
write.table(merged@meta.data,file = paste0("./",sam.name,"/",sam.name,"_cells_details_tsne_res0.3_",max(dim.use),"PC.txt"))

#PCA结果展示-2
pdf(paste0("./",sam.name,"/t_dermis..pdf"),width = 5,height = 4)
DimPlot(merged, reduction = "pca",group.by = 'orig.ident')
dev.off()
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.3_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = T,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.3.nolegend_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res0.3_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="orig.ident", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_sample_res0.3_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="sample", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()

#all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                              #min.pct = 0.25, logfc.threshold = 0.25)
#write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker_genes_tsne_res0.3_",max(dim.use),"PC.txt"),sep="\t",quote = F)
#pdf(file = 'T.marker.pdf',width = 10,height = 8)
FeaturePlot(merged, features = c('CD3D','CD3E','COL1A1','CD3G'),cols = c("gray", "red"),min.cutoff = 0)
#dev.off()
fib.kera.t.merged=subset(merged,seurat_clusters%in%c('10','12'))
save(fib.kera.t.merged,file='fib.kera.t.merged.Rdata')
merged=subset(merged,seurat_clusters!='10')
merged=subset(merged,seurat_clusters!='12')
save(merged,file = 't.removed.merged.Rdata')
library(harmony)

merged=RunHarmony(merged,group.by.vars ='orig.ident')


dim.use=1:30
pdf(paste0("./",sam.name,"/t.removed.dermis.harmony.pdf"),width = 5,height = 4)
DimPlot(merged, reduction = "harmony",group.by = 'orig.ident',cols = cols)
dev.off()

merged<- FindNeighbors(merged, dims = dim.use,reduction = 'harmony')
merged <- FindClusters(merged, resolution =0.3)
merged <- RunTSNE(merged, dims = dim.use,reduction = 'harmony')
pdf(file = 'T.marker.removed.pdf',width = 10,height = 8)
FeaturePlot(merged, features = c('CD3D','CD3E','TRBC2','CD3G'),cols = c("gray", "red"),min.cutoff = 0)
dev.off()

write.table(merged@meta.data,file = paste0("./",sam.name,"/",sam.name,"_cells_details_tsne_res0.3.removed_",max(dim.use),"PC.txt"))

#PCA结果展示-2
pdf(paste0("./",sam.name,"/t_dermis.removed..pdf"),width = 5,height = 4)
DimPlot(merged, reduction = "pca",group.by = 'orig.ident')
dev.off()
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
save(merged,file = 't.removed.tsne.Rdata')
merged$seurat_clusters.1=c(paste0('SC',merged$seurat_clusters))
merged$seurat_clusters.1=factor(merged$seurat_clusters.1,levels = paste0('SC',0:12))
save(merged,file = 't.final.20201203.Rdata')
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.3.removed.nolegend_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = T,cols = cols,group.by = 'seurat_clusters.1')+NoLegend()
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.3.removed_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = T,cols = cols,group.by = 'seurat_clusters.1')
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res0.3.removed.nolegend_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="orig.ident", pt.size=0.1,reduction = "tsne",label = F,cols = cols)+NoLegend()
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res0.3.removed_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="orig.ident", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_sample_res0.3.removed_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="sample", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_sample_res0.3.removed.nolegeng_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="sample", pt.size=0.1,reduction = "tsne",label = F,cols = cols)+NoLegend()
dev.off()

all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker_genes_tsne_res.2.2.renamed.singleR_",max(dim.use),"PC.txt"),sep="\t",quote = F)
all.markers.noRibo=all.markers[!grepl('^RP[SL]',all.markers$gene,ignore.case = F),]
all.markers.noRibo.noMito=all.markers.noRibo[!grepl('^MT-',all.markers.noRibo$gene,ignore.case = F),]
write.table(all.markers.noRibo.noMito,file = 'all.markers.t..res2.2.subcluster.noRI.moMT.txt',sep = '\t')
pdf(file = 'GPR183.pdf',width = 5,height = 4)
FeaturePlot(merged, features = c('GPR183'),cols = c("gray", "red"),min.cutoff = 0)
dev.off()
####12.细胞周期归类 #cc.genes为两个细胞周期的基因集####
merged<- CellCycleScoring(object = merged, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
#head(x = merged@meta.data)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_cellcycle.res0.3removed_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 0.1)
dev.off()
#save(merged,file = 't.dermis.res0.4.dim30.merged.Rdata')
table(Idents(merged))

#merged$orig.ident=factor(merged$orig.ident,levels = c(paste0('HC_',1:3),paste0('DLE_',1:5),paste0('SLE_',1:8)))
#merged$sample=factor(merged$sample,levels = c('HC','DLE','SLE'))
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.3.removed_",max(dim.use),"PC.pdf"),width = 6,height = 5)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols,group.by = 'sample')
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.3.removed_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
save(merged,file = 't.dermis.res2.2.2022.2.20.Rdata')
####15.热图展示Top marker基因####
#筛选top5的marker基因，可以通过参数改为其他数值
marker.sig <- all.markers.noRibo.noMito[all.markers.noRibo.noMito$p_val_adj<=0.05,]

top10 <- marker.sig %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file=paste0("./",sam.name,"/",sam.name,"_top10_marker.res.2.2.nomt.noribo_",max(dim.use),"PC.txt"),sep="\t",quote = F)
#### 11.查看每个cluster的不同来源占比####
a=table(merged$orig.ident,merged$seurat_clusters)
a=as.data.frame(a)
write.table(a,file = 'orig.t.res2.2.subcluster.txt',sep='/t')

p1=ggplot(a,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='T.orig.subcluster.pdf',plot=p3,width=10,height=4)

p1=ggplot(a,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'subcluster',y = 'Frequency ') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                             axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                            panel.background = element_blank(), axis.line =  element_line(colour = "black"))
p2
ggsave(filename = 'orig.subcluster.t.pdf',plot = p2,width = 6,height = 4)
b=table(merged$seurat_clusters,merged$sample)
b=as.data.frame(b)
write.table(b,file = 'sample.t.res2.2.sbucluster.txt',sep = '\t')
p1=ggplot(b,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='T.subcluster.sample.res0.3.pdf',plot=p3,width=5,height=4)

p1=ggplot(b,mapping = aes(Var2,Freq,fill=Var1))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'Frequency ') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                             axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                            panel.background = element_blank(), axis.line =  element_line(colour = "black"))
p2
ggsave(filename = 't.sample.subcluster.pdf',plot = p2,width = 5,height = 4)

c=table(merged$seurat_clusters.1)
c=as.data.frame(c)
p1=ggplot(c, mapping = aes(x = Var1, y = Freq)) + geom_bar(stat = 'identity')+labs(x = 'subcluster',y = 'Number')
p2=p1+scale_fill_manual(values=cols)+theme_bw()+theme(legend.background = element_rect(fill = 'white', colour = 'black'))

p2=p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line.y =  element_line(colour = "black"))
ggsave(filename = 'scT.number.20201203.pdf',plot = p2,width = 6,height = 4)
table(Idents(merged))
Idents(merged)=merged$sample
VlnPlot(merged,features = c('IFI44L','ISG15','IFI27','ISG20','IFI6','IFITM3','IFI44','IFI16'),pt.size = 0,stack = T,)
merged$sample.1=ifelse(merged$sample=='HC','HC','Lupus')
VlnPlot(merged,features = c('CRYAB','TXNRD1'),pt.size = 0)

Idents(merged)=merged$sample.1
all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                                      min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker.sample.1_",max(dim.use),"PC.txt"),sep="\t",quote = F)
genes=all.markers$gene[all.markers$cluster=='Lupus']
genes
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
#BiocManager::install('ReactomePA')
library(ReactomePA)
library(data.table)
####挑出marker基因for GSEA-GO/KEGG分析####
table(merged$seurat_clusters)
a=0:6
for (i in a) {
              gene_list=all.markers[all.markers$cluster==i,]
              genelist_input=gene_list[,c(2,7)]
              
              genename <- as.character(genelist_input[,2])
              genename=mapIds(x = org.Hs.eg.db,keys =genename,
                              keytype = "SYMBOL",column = "ENTREZID")
              # Go_result <- enrichGO(genename, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", pvalueCutoff=0.01)
              geneList = genelist_input[,1]
              names(geneList) = as.character(genename)
              geneList = sort(geneList, decreasing = TRUE)
              Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
              KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
              Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
              go <- setReadable(Go_gseresult, 'org.Hs.eg.db', 'ENTREZID')
              kegg=setReadable(KEGG_gseresult,'org.Hs.eg.db', 'ENTREZID')
              gsea=setReadable(Go_Reactomeresult,'org.Hs.eg.db', 'ENTREZID')
              write.table(go,file = paste0(i,'_t.gsego.txt'),sep='\t')
              write.table(kegg,file =paste0(i,'_t.gsekegg.txt') ,sep='\t')
              write.table(gsea,file =paste0(i,'_t.gsereact.txt') ,sep='\t')
}
####enrichGo/KEGG####
class(all.markers$cluster)
class(a)
for (i in a) {
              gene=all.markers.noRibo.noMito$cluster==i
              gene=all.markers.noRibo.noMito[gene,]
              genename <- as.character(gene$gene)
              genename=mapIds(x = org.Hs.eg.db,keys =genename,
                              keytype = "SYMBOL",column = "ENTREZID")
              enrich.go.BP=enrichGO(gene = genename,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = 'ENTREZID',
                                    ont = 'BP',
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.05,readable = T)
              enrich.KEGG=enrichKEGG(gene = genename,organism = "hsa",
                                     keyType = "kegg",pvalueCutoff = 0.01,
                                     use_internal_data = F,qvalueCutoff = 0.05)
              
              write.table(enrich.go.BP,file =paste0(i,'_T.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_T.kegg.txt') ,sep='\t')
}
####找sample的marker基因####
table(Idents(merged))
Idents(merged)=merged@meta.data$sample
all.markers.sample <- FindAllMarkers(merged, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers.sample,file=paste0("./",sam.name,"/",sam.name,"_sample_marker_genes_t_",30,"PC.txt"),sep="\t",quote = F)
####将marker宽数据改为长数据####
all.markers=read.table('./multi/multi_marker_genes_tsne_res0.4_30PC.txt',sep = '\t',header = T)
top5=read.table('./multi/multi_top5_marker_genes_tsne_30PC.txt',sep = '\t',header = T)
library(pheatmap)
library(reshape2) # 使用的函数 melt & dcast
marker=all.markers[,c(2,6,7)]
marker=marker[marker$gene%in%top5$gene,]
marker<-dcast(marker,marker$gene~marker$cluster,value.var = 'avg_logFC',fill = 0)
row.names(marker)=marker$`marker$gene`
marker=marker[,-1]
write.table(marker,file='T.markertop5.changshuju.txt',sep = '\t')
p1=pheatmap(marker,cluster_cols = T,cluster_rows = T,clustering_distance_rows = "correlation",scale = 'column',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 15,cellheight =15,border_color = 'white',fontsize = 15,legend_labels = 'scale.avg_logFC',legend = T)
p1
ggsave(filename = 'T.subtype.marker.top5.dermis.pdf',plot = p1,width = 20,height = 40)
table(top10$gene)
a=top10$gene
a=a[-61]
table(a)
a=a[-c(57,3,25)]
table(a)
a=a[-2]
Idents(merged)=merged$seurat_clusters
pdf('t.subcluster,dotplot.pdf',width = 15,height = 7)
DotPlot(merged, features = a)+RotatedAxis()+
              scale_x_discrete("")+scale_y_discrete("")
dev.off()
####marker基因展示####
top10$subcluster=c(paste0('SC',top10$cluster))
table(top10$subcluster)
marker=merged@misc$averageExpression$RNA
colnames(marker)=c(paste0('SC',0:6))
a=c(paste0('SC',0:6))
for (i in a) {genelist=top10[top10$subcluster==i,]
genelist=marker[row.names(marker)%in%genelist$gene,]
write.table(genelist,file = paste0('genelist',i,'.txt'),sep = '\t')
}
markers=read.table('./top10.marker.subcluster.txt',sep = '\t',header = T,row.names = 1)

markers=t(markers)
library(pheatmap)
p1=pheatmap(markers,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'column',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')
p1
ggsave(filename = 'T.subtype.top10.marker.dermis.expression.pdf',plot = p1,width = 45,height = 20)

####经典marker基因匹配####

classical.marker=read.table('./classical',sep = '\t',header = T,stringsAsFactors = F,row.names = 1)
naive_T=classical.marker$na.ve
Tcm=classical.marker$Tcm

Tem=classical.marker$Tem
Trm=classical.marker$Trm
MAIT=classical.marker$MAIT
Exhua_T=classical.marker$exhuastion.T
Tfh=classical.marker$Tfh
Th1_like=classical.marker$Th1.like.T
Th17=classical.marker$Th17
Tfr=classical.marker$Tfr
Th2=classical.marker$Th2
Treg=classical.marker$Treg
CTL=classical.marker$CTL




a=c(naive_T,Tcm,Tem,Trm,MAIT,Exhua_T,Tfh,Th1_like,Th17,Tfr,Th2,Treg,CTL,CD4_T,CD8_T)
class(a)
row.names(classical.marker)
annotation_col=data.frame(gene = row.names(classical.marker),subtype = c(rep('ISG_T',8),rep('CD4_T',1),rep('CD8_T',2),
                                      rep('Treg',3),
                                      rep('Exhausted_T',2),
                                      rep('HSP_T',5),
                                      rep('Th2',4),
                                      rep('Th1',2),
                                      rep('TH17',3),
                                      rep('CTL',7),
                                      rep('Naive_T',3),
                                      rep('Memory_T',3),
                                      rep('Resident_Memory_T',1),
                                      rep('Tfh',3)))
all.markers=read.table('./multi/multi_marker_genes_tsne_res0.7.removed_30PC.txt',sep = '\t',header = T,stringsAsFactors = F)
a=read.table('./average.expression.subt.removed.cluster.txt',sep = '\t',header = T,stringsAsFactors = F)
colnames(a)=0:12
select.marker=all.markers[all.markers$gene%in%row.names(annotation_col),]

#annotation_col$gene%in%row.names(select.marker)
#marker.3=top10[top10$cluster=='3',]
row.names(select.marker)=reorder(annotation_col$gene,row.names(select.marker))
select.marker=all.markers[all.markers$gene%in%annotation_col$gene,]
#select.marker=select.marker[,c(3,7,8)]
select.marker<-dcast(select.marker,select.marker$gene~select.marker$cluster,value.var = 'avg_logFC',fill = 0)

row.names(select.marker)=select.marker$`select.marker$gene`
select.marker=select.marker[,-1]
annotation_col=annotation_col[annotation_col$gene%in%row.names(select.marker),]
write.table(annotation_col,file = 'annotation_col.txt',sep = '\t')
#annotation_col=annotation_col[unique(annotation_col$gene),]
annotation_col=read.table('annotation_col.txt',sep = '\t',header = T,stringsAsFactors = F)

row.names(annotation_col)=annotation_col$gene

#annotation_col=annotation_col[,-c(1,2,4)]
#dim(annotation_col)
annotation_col=annotation_col[annotation_col$gene%in%row.names(select.marker),]
annotation_col=as.data.frame(annotation_col$subtype.1,row.names = row.names(annotation_col))
colnames(annotation_col)='subtype'
row.names(select.marker)=reorder(row.names(annotation_col),row.names(select.marker))

#class(annotation_col$subtype)=as.character(annotation_col$subtype)
annotation_col=annotation_col[,-c(1,2)]
annotation_col=as.data.frame(annotation_col,row.names = row.names(select.marker))
colnames(annotation_col)='subtype'
p1=pheatmap(select.marker,cluster_rows = F,cluster_cols = F,annotation_row =annotation_col ,clustering_distance_rows = "correlation",scale = 'column',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 20,cellheight =12,border_color = 'white',fontsize = 10,legend_labels = 'scale.avg_logFC',legend = T)
p1
ggsave(filename = 't.classical.subtype.marke.dermis.pdf',plot = p1,width = 6,height = 15)
####与经典基因匹配####
merged@misc$averageExpression=AverageExpression(merged)
marker=merged@misc$averageExpression$RNA
table(Idents(merged))
classical.gene=read.table('./classical.marker.txt',header = T,row.names = 1)
classical.gene$genes=row.names(classical.gene)
th17=classical.gene$genes[classical.gene$subtype=='Th17']
th17=th17[-2]
th17=marker[row.names(marker)%in%th17,]
ctl=classical.gene$genes[classical.gene$subtype=='CTL']
ctl=marker[row.names(marker)%in%ctl,]
naive_T=classical.gene$genes[classical.gene$subtype=='Naive_T']
naive_T=naive_T[-3]
naive_T=marker[row.names(marker)%in%naive_T,]
memory_t=classical.gene$genes[classical.gene$subtype=='Memory_T']
memory_t=memory_t[-2]
memory_t=marker[row.names(marker)%in%memory_t,]
resident=classical.gene$genes[classical.gene$subtype=='Resident_memory_T']
resident=marker[row.names(marker)%in%resident,]
#tfh=classical.gene$genes[classical.gene$subtype=='Tfh']
#tfh=marker[row.names(marker)%in%tfh,]
treg=classical.gene$genes[classical.gene$subtype=='Treg']
treg=marker[row.names(marker)%in%treg,]
#th2=classical.gene$genes[classical.gene$subtype=='Th2']
#th2=marker[row.names(marker)%in%th2,]
#exhua_T=classical.gene$genes[classical.gene$subtype=='Exhausted_T']
#exhua_T=marker[row.names(marker)%in%exhua_T,]
th1=classical.gene$genes[classical.gene$subtype=='Th1']
th1=marker[row.names(marker)%in%th1,]
cd4=classical.gene$genes[classical.gene$subtype=='CD4']
cd4=marker[row.names(marker)%in%cd4,]
cd8=classical.gene$genes[classical.gene$subtype=='CD8']
cd8=marker[row.names(marker)%in%cd8,]
#hsp=classical.gene$genes[classical.gene$subtype=='HSP_T']
#hsp=marker[row.names(marker)%in%hsp,]
isg=classical.gene$genes[classical.gene$subtype=='ISG_T']
isg=isg[-c(3,5)]

isg=marker[row.names(marker)%in%isg,]
classical=rbind(th17,ctl,naive_T,memory_t,resident,treg,th1,cd4,cd8,isg)
colnames(classical)=c(paste0('SC',0:6))
row.names(classical)[14]='ITGAE'
row.names(classical)[20]='CD4'
library(pheatmap)
classical.gene=data.frame(subtype=classical.gene$subtype,row.names = classical.gene$genes)
annotation.col=data.frame(sub=c('CD4','CD8','CD8','CD4','CD4','CD8','CD8'),row.names = c(paste0('SC',0:6)))
p1=pheatmap(classical,annotation_row = classical.gene,scale = 'row',cluster_rows = F,cluster_cols = T,colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color = 'white',cellwidth = 20,cellheight = 12,clustering_distance_cols = 'correlation',treeheight_col = 20)
ggsave(filename = 'classical.t.subcluster.genes.pdf',plot = p1,width = 20,height = 40)

####Go热图####
####EC.GO富集分析画图####
getwd()
GO=read.table('./T.all.go.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
GO.select=read.table('./T.select.go.txt',sep = '\t',header = T,stringsAsFactors = F)

head(GO)


GO=GO[,c(3,8,10,11)]



GO2=GO[,c(1,3,4)]


GO2=GO2[GO2$Description%in%GO.select$Description,]
library(pheatmap)
library(reshape2) # 使用的函数 melt & dcast

GO3<-dcast(GO2,Description~GO2$subcluster,value.var = 'qvalue',fill = 0)
GO4<-dcast(GO2,Description~GO2$cluster,value.var = 'Count',fill = 0)
#heatmap(GO3)
row.names(GO4)=GO4$Description
GO4=GO4[,-1]

p1=pheatmap(GO4,scale = 'row',clustering_distance_rows = "correlation",colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 50,cellheight =25,border_color = 'white',fontsize = 25,legend_labels = 'scale.avg_logFC',legend = T,angle_col = 45)
p1

ggsave(filename = 't.GO.heatmap.count.pdf',plot = p1,width = 30,height = 30)
####找sample的marker基因####
table(Idents(merged))
Idents(merged)=merged@meta.data$sample
all.markers.sample <- FindAllMarkers(merged, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers.sample,file=paste0("./",sam.name,"/",sam.name,"_sample.marker_genes_tsne_res0.4_",30,"PC.txt"),sep="\t",quote = F)
all.markers.sample=read.table('./multi/multi_sample.marker_genes_tsne_res0.4_30PC.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
a=c('HC','DLE','SLE')
for (i in a) {
              gene_list=all.markers.sample[all.markers.sample$cluster==i,]
              genelist_input=gene_list[,c(2,7)]
              
              genename <- as.character(genelist_input[,2])
              genename=mapIds(x = org.Hs.eg.db,keys =genename,
                              keytype = "SYMBOL",column = "ENTREZID")
              # Go_result <- enrichGO(genename, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", pvalueCutoff=0.01)
              geneList = genelist_input[,1]
              names(geneList) = as.character(genename)
              geneList = sort(geneList, decreasing = TRUE)
              Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
              KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
              Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
              go <- setReadable(Go_gseresult, 'org.Hs.eg.db', 'ENTREZID')
              kegg=setReadable(KEGG_gseresult,'org.Hs.eg.db', 'ENTREZID')
              gsea=setReadable(Go_Reactomeresult,'org.Hs.eg.db', 'ENTREZID')
              write.table(go,file = paste0(i,'_t.sample.gsego.txt'),sep='\t')
              write.table(kegg,file =paste0(i,'_t.sample.gsekegg.txt') ,sep='\t')
              write.table(gsea,file =paste0(i,'_t.sample.gsereact.txt') ,sep='\t')
}
####enrichGo/KEGG####
for (i in a) {
              gene=all.markers.sample$cluster==i
              gene=all.markers.sample[gene,]
              genename <- as.character(gene$gene)
              genename=mapIds(x = org.Hs.eg.db,keys =genename,
                              keytype = "SYMBOL",column = "ENTREZID")
              enrich.go.BP=enrichGO(gene = genename,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = 'ENTREZID',
                                    ont = 'BP',
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.05,readable = T)
              enrich.KEGG=enrichKEGG(gene = genename,organism = "hsa",
                                     keyType = "kegg",pvalueCutoff = 0.01,
                                     use_internal_data = F,qvalueCutoff = 0.05)
              
              write.table(enrich.go.BP,file =paste0(i,'_T.sample.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_T.sample.kegg.txt') ,sep='\t')
}
####sample.select.go展示####
####T.GO富集分析画图####
####气泡图####
GO=read.table('./t.dermis.all.go.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
GO=GO[,c(3,8,10,11)]
GO1=read.table('./t.dermis.selected.go.txt',sep = '\t',header = T,stringsAsFactors = F)
GO7=GO[GO$Description%in%GO1$Description,]
library(forcats)
GO7$Description=fct_infreq(GO7$Description)

#GO7$subcluster=paste0('SC',GO7$subcluster)
#write.table(GO7,file = 'subT.GO7.burble.txt',sep = '\t')
GO7$sample=factor(GO7$sample,levels = c('HC','DLE','SLE'))
p <- ggplot(GO7, aes(sample,Description, size = Count, color=qvalue)) + geom_point() +xlab("sample")
#### 气泡图颜色####
p2=p + scale_color_gradient(low='red',high='blue') +theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.major=element_line(colour=NA))+theme(text = element_text(size = 20))
p2
ggsave(filename = 'burbble.pathway.sample.t.20220217.pdf',plot = p2,width =10,height = 10)

p<-ggplot(GO7,aes(sample,Description)) +
              geom_point(aes(fill=qvalue,size=Count),alpha=0.9,pch=21,colour="gray25") +  #fill对应点的填充色，colour对应点的边框色
              scale_fill_gradient(low='red', high='blue')+ #设定颜色的变化范围
              labs(x='sample',fill='qvalue',size='Count')
p2=p+ theme(legend.background = element_rect(fill = 'white', colour = 'white'))
p2

GO8=read.table('./subT.GO7.burble.txt',sep = '\t',header = T)
GO8$subcluster=c(paste0('SC',GO8$cluster))
GO8$Description.1=paste0('GO:',GO8$Description)
class(GO8$Description)
library(forcats)
GO8$Description=fct_infreq(GO8$Description)
levels(GO8$Description)
GO8$subcluster=factor(GO8$subcluster,levels = c('SC0','SC3','SC4','SC1','SC2','SC5'))
mytheme <- theme(axis.title=element_text(face="bold", size=10,colour = 'gray25'), #坐标轴标题
                 axis.text=element_text(face="bold", size=10,colour = 'gray25'), #坐标轴标签
                 axis.line = element_line(size=0.5, colour = 'black'), #轴线
                 panel.background = element_rect(color='black'), #绘图区边框
                 legend.key = element_blank() #关闭图例边框
)
p<-ggplot(GO8,aes(subcluster,Description)) +
              geom_point(aes(fill=qvalue,size=Count),alpha=0.9,pch=21,colour="gray25") +  #fill对应点的填充色，colour对应点的边框色
              scale_fill_gradient(low='red', high='blue')+ #设定颜色的变化范围
              labs(y='GO_BP',x='subcluster',fill='qvalue)',size='Count')

p2=p+ theme(legend.background = element_rect(fill = 'white', colour = 'white'))+mytheme
p <- ggplot(GO8, aes(subcluster,Description, size = Count, color=qvalue)) + geom_point() +xlab("subcluster")
#### 气泡图颜色####
p2=p + scale_color_gradient(low='red',high='blue') +theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.major=element_line(colour=NA))
ggsave(filename = 'burbble.pathway.subT20201127.pdf',plot = p2,width = 15,height = 20)
####
getwd()
GO=read.table('b.sample.all.go.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
GO.select=read.table('b.sample.select.go.txt',sep = '\t',header = T,stringsAsFactors = F)

head(GO)


GO=GO[,c(3,8,10,11)]



GO2=GO[,c(1,3,4)]


GO2=GO2[GO2$Description%in%GO.select$Description,]
library(pheatmap)
library(reshape2) # 使用的函数 melt & dcast

GO3<-dcast(GO2,Description~GO2$subcluster,value.var = 'qvalue',fill = 0)
GO4<-dcast(GO2,Description~GO2$sample,value.var = 'Count',fill = 0)
#heatmap(GO3)
row.names(GO4)=GO4$Description
GO4=GO4[,-1]

p1=pheatmap(GO4,scale = 'row',clustering_distance_rows = "correlation",colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 50,cellheight =20,border_color = 'white',fontsize = 20,legend_labels = 'scale.avg_logFC',legend = T,angle_col = 45)
p1

ggsave(filename = 'b.sample.GO.heatmap.count.pdf',plot = p1,width = 20,height = 15)
####cluster基因平均表达水平####
table(Idents(merged))
#Idents(merged)=merged$seurat_clusters
merged@misc$averageExpression=AverageExpression(merged)
write.table(merged@misc$averageExpression$RNA,file = 'sample.average.expression.subt.removed.cluster.txt',sep = '\t',row.names = T,quote = F)
b=merged@misc$averageExpression
b=as.data.frame(b)
b=scale(b)
b=b[row.names(b)%in%top10$gene,]
write.table(b,file = 'scale.average.exp.subt.top10.20201203.txt',sep = '\t')
colnames(b)=c(paste0('SC',0:12))
b=t(b)
library(pheatmap)
p1=pheatmap(b,cluster_cols = F,clustering_distance_rows = "correlation",scale = 'row',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 12,cellheight =12,border_color = 'white',fontsize = 10,legend_labels = 'scale.avg_logFC',legend = T,silent = T)

ggsave(filename = 'top10.subt.scale.exp.pdf',plot = p1,width = 15,height = 15)

Idents(merged)=merged$sample
merged@misc$averageExpression=AverageExpression(merged)
write.table(merged@misc$averageExpression$RNA,file = 'average.expression.sample.removed.cluster.txt',sep = '\t',row.names = T,quote = F)
Idents(merged)=merged$seurat_clusters
merged$subt.smaple=paste(Idents(merged),merged$sample,sep = '.')
Idents(merged)=merged$subt.smaple
table(Idents(merged))

merged@misc$averageExpression=AverageExpression(merged)
write.table(merged@misc$averageExpression$RNA,file = 'average.expression.subt.removed.cluster.txt',sep = '\t',row.names = T,quote = F)

pdf(file = 'RIGI.MDA5.1.pdf',width = 15,height = 16)
FeaturePlot(merged, features = c('DDX58','IFIH1','MAVS','TBK1','IRF3','IRF7','NFKB1','TRIM25'),cols = c("gray", "red"),min.cutoff = 0)
dev.off()
pdf(file = 'RIGI.MDA5.2.pdf',width = 15,height = 16)
FeaturePlot(merged, features = c('PRKCA','PRKCB','CSNK2A1','PPP1CA','PPP1CC','HDAC6','TRAF2','TRAF5'),cols = c("gray", "red"),min.cutoff = 0)
dev.off()
pdf(file = 'RIGI.MDA5.3.pdf',width = 15,height = 16)
FeaturePlot(merged, features = c('TRAF6','IKBKB','CHUK','IFNA1','IFNB1','STAT1','STAT2'),cols = c("gray", "red"),min.cutoff = 0)
dev.off()
####RIG1/MDA5gene expresion####
rig.gene.expre=read.table('./RIG1.MDA5.genes.expre.cluster.txt',sep = '\t',header = T,stringsAsFactors = F)
row.names(rig.gene.expre)=rig.gene.expre$X
rig.gene.expre=rig.gene.expre[,-1]
colnames(rig.gene.expre)=0:12
p1=pheatmap(rig.gene.expre,cluster_cols = T,cluster_rows = T,clustering_distance_rows = "correlation",scale = 'row',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 15,cellheight =15,border_color = 'white',fontsize = 15,legend_labels = 'scale.avg_logFC',legend = T)
p1
ggsave(filename = 'rig.mda5.gene.t.cluster.expression.pdf',plot = p1,width = 20,height = 40)
####singleR####
b=merged@assays$RNA@counts
b=SummarizedExperiment(assays=list(counts=b,logcounts=b))

b <- logNormCounts(b)
library(SingleR)
####HPCA参考数据####
hpca.se <- HumanPrimaryCellAtlasData()
pred.hesc <- SingleR(test = b, ref = hpca.se, labels = hpca.se$label.main)
pred.hesc.1 <- SingleR(test = b, ref = hpca.se, labels = hpca.se$label.fine)
load('E:/R_analysis/BD/hpca.se.Rdata')
save(hpca.se,file='hpca.se.Rdata')
#hpca.se=subset(hpca.se,hpca.se$label.main=='T_cells')
hpca.se
pred.hesc
view(hpca.se@colData)
table(pred.hesc$labels)
table(pred.hesc.1$labels)
####将singleR结果匹配到seurat对象####
merged@meta.data$cell_type=pred.hesc$labels
merged@meta.data$subtype=pred.hesc.1$labels
table(merged$orig.ident)
load('E:/R_analysis/BD/mana20200817.Rdata')
mana
table(mana$label.fine)
View(mana@metadata)
####mana参考数据####
pred.mana <- SingleR(test = b, ref = mana, labels = mana$label.main)
pred.mana.1 <- SingleR(test = b, ref = mana, clusters = merged$seurat_clusters,labels = mana$label.fine)
table(pred.mana.1$labels)
merged@meta.data$cell_type=pred.mana$labels
merged@meta.data$subtype=pred.mana.1$labels
table(merged$orig.ident)
a=table(merged$seurat_clusters,merged$cell_type,merged$subtype)
a=data.frame(a)
write.table(a,file = 't.singler.celltype.txt',sep='\t')

####重新画图加配色####
library(RColorBrewer)
Idents(merged)=merged$cell_type
table(merged$cell_type)
Idents(merged)=merged$subtype
Idents(merged)=factor(Idents(merged),levels = c('CD4+ T cells','CD8+ T cells','B cells','T cells','Neutrophils','NK cells','Monocytes','Dendritic cells','Progenitors','Basophils'))
cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
#cols=brewer.pal(12,'Paired')

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res2.2.singleR.mana_",max(dim.use),"PC.pdf"),width = 6,height = 5)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols)
dev.off()
Idents(merged)=merged$subtype
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.7.singleR.mana.1_",max(dim.use),"PC.pdf"),width = 10,height = 8)
DimPlot(object = merged, pt.size=0.1,label = F)
dev.off()
####singleR结果画图####
t.celltype=read.table('./t.singler.celltype.txt',sep = '\t',header = T,stringsAsFactors = F)
t.celltype=t.celltype[order(t.celltype[,1]),]
t.subtype=t.celltype[,c(1,3,4)]
GO4<-dcast(GO2,Description~GO2$sample,value.var = 'Count',fill = 0)
t.subtype=dcast(t.subtype,Var1~t.subtype$Var3,value.var = 'Freq',fill = 0,fun.aggregate = sum)
row.names(t.subtype)=t.subtype$Var1
t.subtype=t.subtype[,-1]
write.table(t.subtype,file = 'subt.singleR.subtype.number.txt',sep = '\t')
p1=pheatmap(t.subtype,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 15,cellheight =15,border_color = 'white',fontsize = 15,legend_labels = 'scale.avg_logFC',legend = T)
p1
ggsave(filename = 't.singleR.subtype.pdf',plot = p1,width = 20,height = 40)
####SCENIC####

library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
rm(list=ls())

####monocle####

saveRDS(merged, file = "./multi/T_tutorial_res0.3.rds")
save(merged,file="./multi/T_res0.3.Robj")
rm(list = ls())
load('./multi/T_res0.3.Robj')
load('./t.downsample.nishixu.Rdata')
####分层抽样分析####
####1.用sample无放回的抽取####
rm(list = ls())
cells=colnames(merged)
cell_id <- sample(x = seq(1, length(cells)), 3000)
cell=cells[cell_id]
cell
merged.1=subset(merged,colnames(merged)%in%cell)
a=colnames(merged)%in%cell
merged$downsample=a
table(merged$downsample)
merged.1=subset(merged,downsample==T)
save(merged.1,file = 't.downsample.nishixu.Rdata')
table(merged$seurat_clusters)
table(merged.1$seurat_clusters)
####2.根据细胞的元信息，分层抽样####
# 获取细胞编号
cell_ids <- row.names(cells)
# 获取细胞来源
cell_source <- factor(cells$donor_organism.human_specific.ethnicity.ontology_label)
# 根据细胞来源换分细胞，然后抽样
cell_ids_list <- lapply(split(cell_ids, cell_source), function(x){ sample(x, 3000) })
# 提取细胞编号
cell_id <- unlist(cell_ids_list)
# 提取细胞
cell_sml <- cells[cell_id, ]
mtx_sml <- mtx[, cell_id]

table(Idents(merged))
markers=read.table('./multi/multi_marker_genes_tsne_res0.3.removed_30PC.txt',sep='\t')
markers <- FindAllMarkers(merged.1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Idents(merged.1)=merged.1$sample
markers.1=markers[markers$pct.1>0.5&markers$pct.2<0.5,]
#top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
rm(merged)
save(merged,file="./multi/merged_f3.Robj") 

data <- as(as.matrix(merged.1@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = merged.1@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
markers=read.table('./multi/nishixu.ordering.genes.txt',sep = '\t',row.names = 1,header = T)
#Construct monocle cds
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size());


cds <- estimateSizeFactors(cds)
rm(merged)
cds <- estimateDispersions(cds)
#expressed_genes <- row.names(subset(fData(cds)))
save(cds,file="./multi/t_cds_normal.Robj")

#ordering_genes<-as.matrix(top20)

#diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~percent.mt")

ordering_genes <- row.names (subset(markers, p_val < 0.01))
cds.1 <- setOrderingFilter(cds, ordering_genes)

pdf('./multi/plot_ordering_genes_subt.pdf')
plot_ordering_genes(cds.1)
dev.off()

#Trajectory step 2: reduce data dimensionality  
HSMM <- reduceDimension(cds.1, max_components = 2,
                        method = 'DDRTree')

#Trajectory step 3: order cells along the trajectory  
HSMM <- orderCells(HSMM)
pdf('./multi/plot_cell_trajectory_seurat_clusters_sample.pdf')
plot_cell_trajectory(HSMM, color_by = "sample",show_branch_points = F)
dev.off()
pdf('./multi/plot_cell_trajectory_state_subt3.pdf')
plot_cell_trajectory(HSMM, color_by = "State",show_branch_points = F)
dev.off()
HSMM=orderCells(HSMM,root_state = 3)
pdf('./multi/plot_cell_trajectory_pseudotime_subt3.pdf')
plot_cell_trajectory(HSMM, color_by = "Pseudotime",show_branch_points = F)
dev.off()
####找差异基因####
HSMM=detectGenes(HSMM,min_expr = 0.1)
expressed_genes=row.names(subset(fData(HSMM),num_cells_expressed>=10)) #在部分基因里面找
pseudotime_de <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
states_de <- differentialGeneTest(HSMM[expressed_genes,],
                                  fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]
table(states_de$use_for_ordering)
saveRDS(HSMM, file = "t_monocle.rds")
write.table(pseudotime_de, file = "t.dermis.pseudotime_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(states_de, file = "t.dermis.states_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
####分支点分析####
BEAM_res=BEAM(HSMM,branch_point = 1,cores = 1)
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")
#cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
#cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree',auto_param_selection = F) # take long time
#cds <- orderCells(cds)
tmp1=plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 4, #这些基因被分成几个group
                                 cores = 1,
                                 branch_labels = c('Cell_fate1','Cell_fate2'),
                                 #hmcols = NULL, #默认值
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #是否返回一些重要信息
)
pdf("branched_heatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()
save(HSMM,file="./sub_t3_orderCells.Robj")
load('./multi/sub_t3_orderCells.Robj')
write.table(pData(cds),file="./multi/subt3_my_pseudotime.txt")
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")

pdf('./multi/plot_cell_trajectory_pseudotime_no_points_subt3.pdf')
plot_cell_trajectory(HSMM, color_by = "Pseudotime",show_branch_points = F)
dev.off()

pdf('./multi/plot_cell_trajectory_subtype_subt3.pdf')
plot_cell_trajectory(HSMM, color_by = "seurat_clusters.1",show_branch_points = F)
dev.off()
#?plot_cell_trajectory
pdf('./multi/plot_cell_trajectory_pseudotime_heatmap_subt3.row.names.pdf')
plot_pseudotime_heatmap(HSMM[row.names(HSMM)%in%ordering_genes,],
                        num_clusters =3,
                        cores = 1,show_rownames = F)
dev.off()
table(markers$cluster)
marker1=markers$gene[markers$cluster=='HC']
marker1=marker1[1:10]
marker2.1=markers$gene[markers$cluster=='DLE']
marker2.1=marker2.1[1:10]
marker2.2=markers$gene[markers$cluster=='SLE']
marker2.2=marker2.2[1:10]
marker2=c(marker2.1,marker2.2)
class(marker2.1)=as.factor(marker2.1)
marker=c(marker1,marker2)
t_monocle
load('./sub_t3_orderCells.Robj')
pdf('./multi/plot_cell_trajectory_pseudotime_heatmap_subt.1.pdf',width = 6,height = 6)
plot_pseudotime_heatmap(HSMM[row.names(HSMM)%in%marker,],show_rownames = T,num_clusters = 3)
dev.off()
write.table(markers,file = './multi/nishixu.ordering.genes.txt',sep = '\t',row.names = T)
#plot_genes_branched_heatmap(HSMM[row.names(subset(HSMM,
pval < 1e-4)),],
branch_point = 1,
num_clusters = 4,
cores = 1,
use_gene_short_name = T,
show_rownames = T)

#plot_multiple_branches_heatmap(HSMM, branches = c(1,2,3),
cluster_rows = TRUE, hclust_method = "ward.D2", num_clusters = 6,
hmcols = NULL, add_annotation_row = NULL, add_annotation_col = NULL,
show_rownames = FALSE, use_gene_short_name = TRUE,
norm_method = c("vstExprs", "log"), scale_max = 3, scale_min = -3,
trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
cores = 1)
#?plot_multiple_branches_heatmap
genes=c('ISG15','IFITM3')
pdf("genes_branched_pseudotime.pdf",width = 9,height = 4)
plot_genes_branched_pseudotime(HSMM[genes,],
                               branch_point = 1,
                               color_by = "sample",
                               cell_size=2,
                               ncol = 2)
dev.off()

## Azimuth预测
sc.test <-merged@assays$RNA@counts   #APP会做SCT标准化
saveRDS(sc.test, "sc.test.rds")      #保存为rds格式上传
#上传http://azimuth.satijalab.org/app/azimuth
predictions <- read.delim('CellType/azimuth_pred.tsv', row.names = 1)
scRNA <- AddMetaData(scRNA, metadata = predictions)
p1 <- DimPlot(scRNA, reduction = "umap", label = T) + NoLegend()
p2 <- DimPlot(scRNA, reduction="umap", group.by="predicted.id", label=T,
              label.size = 4) + ggtitle("Predicted by Azimuth") + NoLegend()
p <- p1 + p2
ggsave("CellType/CellType_Azimuth.png", p, width = 10, height = 4)


p <- DimPlot(scRNA, reduction="umap", group.by="predicted.id", split.by = "predicted.id", ncol = 3) + 
              ggtitle("Predicted by Azimuth") + NoLegend()
ggsave("CellType/CellType_Azimuth_facet.png", p, width = 10, height = 25)
p <- DimPlot(scRNA, reduction="umap", group.by="predicted.id", split.by = "predicted.id", ncol = 3) + 
              ggtitle("Predicted by Azimuth") + NoLegend()
ggsave("CellType/CellType_Azimuth_facet.png", p, width = 10, height = 25)

##===Marker基因鉴定===##
##提取各个Cluster的marker genes
ClusterMarker <- FindAllMarkers(scRNA, assay = "SCT", slot = "data", only.pos = T,
                                logfc.threshold = 0.25, min.pct = 0.1)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'CellType/ClusterMarker.csv', row.names=F)
#ClusterMarker <- read.csv('CellType/ClusterMarker.csv')
#提取差异显著的marker genes
top = 15   #可根据需要调整
TopMarkers1 <- ClusterMarker %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% 
              top_n(n = top, wt = avg_logFC)
TopMarkers2 <- ClusterMarker %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>%
              top_n(n = top, wt = avg_logFC)
TopMarkers <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers,'CellType/TopMarkers.csv', row.names=F)

##提取没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", 
                                             ClusterMarker$gene, ignore.case = F),]
top = 15   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo %>% filter(p_val_adj == 0) %>% 
              group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers2 <- ClusterMarker_noRibo %>% filter(p_val_adj < 0.01) %>% 
              group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers_noRibo <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% 
              arrange(cluster)
write.csv(TopMarkers_noRibo,'CellType/TopMarkers_noRibo.csv', row.names=F)

##导出一个方便对比的表格(每行一个cluster的marker gene)
TopMarkers2Lines(data = TopMarkers_noRibo, output = "CellType")

##辅助查找marker基因
source("Resource/MarkerGeneList.R")
source("Resource/function.R")
#TopMarkers_noRibo <- read.csv('CellType/TopMarkers_noRibo.csv')
MatchMarkerGene(data=TopMarkers_noRibo, output = "CellType")
####cellphoneDB####
merged$subtype=c(paste0('T_',merged$seurat_clusters.1))
merged$subtype=factor(merged$subtype,levels = c(paste0('T_SC',0:12)))
t_merged=merged
load('/home/zhaolab/R_analasis/dermis.harmony/dermis.res1/fib.dermis/fib.20201115/multi/fib.dermis.res0.3.20201115.Rdata')
merged$seurat_clusters.1=c(paste0('SC',merged$seurat_clusters))
merged$seurat_clusters.1=factor(merged$seurat_clusters.1,levels=c(paste0('SC',0:13)))
merged$subtype=c(paste0('F_',merged$seurat_clusters.1))
merged$subtype=factor(merged$subtype,levels = c(paste0('F_SC',0:13)))
fib_merged=subset(merged,seurat_clusters.1%in%c('SC1','SC3','SC5','SC6','SC11'))
cellphoneDB=merge(fib_merged,t_merged)
Idents(cellphoneDB)=cellphoneDB$subtype
Idents(cellphoneDB)=factor(Idents(cellphoneDB),levels = c('F_SC1','F_SC3','F_SC5','F_SC6','F_SC11',paste0('T_SC',0:12)))
count_raw=as.matrix(cellphoneDB@assays$RNA@counts)
count_norm=apply(count_raw,2,function(x)(x/sum(x))*10000)
count_norm=data.frame(Gene=rownames(count_norm),count_norm,check.names = F)
write.table(count_norm,'cellphonedb.count.txt',sep = '\t',quote = F,row.names = F)


meta.data=data.frame(cell=rownames(cellphoneDB@meta.data),cell_type=cellphoneDB$cell_type,subtype=cellphoneDB$subtype)
write.table(meta.data,'cellphonedb.metadata.txt',sep = '\t',quote = F,row.names = F)

merged$renamed=ifelse(merged$seurat_clusters%in%c(5,7,8,23,24),'Activated CD4 T1',
                      ifelse(merged$seurat_clusters==1,'Activated CD4 T2',
                             ifelse(merged$seurat_clusters%in%c(0,4,11,14,21,18,28),'Naive CD4 T',
                                    ifelse(merged$seurat_clusters%in%c(2,25,27),'CTL2',
                                           ifelse(merged$seurat_clusters%in%c(3,10,16,17),'CTL1',
                                                  ifelse(merged$seurat_clusters%in%c(6,13,20),'Tfh',
                                                         ifelse(merged$seurat_clusters%in%c(12,15,22),'Treg',
                                                                ifelse(merged$seurat_clusters%in%c(19,26),'Memory CD4','Memory CD8'))))))))
merged$renamed=fct_infreq(merged$renamed)
table(merged$renamed)
table(merged$seurat_clusters)
Idents(merged)=merged$renamed
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res2.2.renamed_",max(dim.use),"PC.pdf"),width = 7,height = 5)
DimPlot(object = merged, pt.size=0.1,reduction = "tsne",label = T,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res2.2_",max(dim.use),"PC.pdf"),width = 7,height = 5)
DimPlot(object = merged, group.by="orig.ident", pt.size=0.1,reduction = "tsne",label = F)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_sample_res2.2_",max(dim.use),"PC.pdf"),width = 7,height = 5)
DimPlot(object = merged, group.by="sample", pt.size=0.1,reduction = "tsne",label = F)
dev.off()

b=table(merged$sample,merged$renamed)
b=as.data.frame(b)
write.table(b,file = 'sample.renamed.res2.2.txt',row.names = F,sep = '\t')
c=table(merged$orig.ident,merged$renamed)

c=as.data.frame(c)
write.table(c,file = 'orig.renamed.res2.2.txt',row.names = F,sep = '\t')


#d$Var1=factor(d$Var1,levels = c('HC_1','HC_2','HC_3','DLE_1','DLE_2','SLE_1','SLE_2','SLE_3'))
p1=ggplot(b,mapping = aes(Var2,Freq,fill=Var1))+geom_bar(stat='identity',position='fill') +
              labs(x = 'subclusters',y = 'Frequency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                             axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='renamed.t.res2.2.sample.fill.pdf',plot=p3,width=8,height=5)
p1=ggplot(c,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'subclusters',y = 'Frequency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                             axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='sample.renamed.t.res2.2.pdf',plot=p3,width=8,height=5)
all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_renamed.res2.2.t.marker_genes_tsne_",max(dim.use),"PC.txt"),sep="\t",quote = F)
all.markers.noRibo=all.markers[!grepl('^RP[SL]',all.markers$gene,ignore.case = F),]
all.markers.noRibo.noMito=all.markers.noRibo[!grepl('^MT-',all.markers.noRibo$gene,ignore.case = F),]
write.table(all.markers.noRibo.noMito,file = 'all.markers.renamed.res2.2.t.noRI.moMT.txt',sep = '\t')
table(merged$renamed)
DotPlot(merged, features = c('IL7R','TCF7','CD27','CD69','PRF1','HSPA1B','GZMB','ICOS','PDCD1','FOXP3','IL2RA','IFNG','GZMA','NFKBIA','TNFAIP3','SELL','CCR7','CD8A'))+RotatedAxis()+
              scale_x_discrete("")+scale_y_discrete("")
####真表皮T与PBMC中的T比较####
####表皮的T与PBMC中的T比较####
getwd()
table(all.markers$cluster)
table(Idents(merged))
all.markers.pbmc=read.table('./multi/marker.genes.t.pbmc.txt',header = T)
all.markers.dermis=read.table('./multi/multi_marker_genes_tsne_res0.3.removed_20PC.txt',sep = '\t',header = T)
table(all.markers.dermis$cluster)
table(merged$seurat_clusters)
all.markers.epi=read.table('./multi/epi.t.multi_marker_genes_tsne_res0.15.removed_20PC.txt',header = T,row.names = 1)
table(all.markers.epi$cluster)
table(all.markers.pbmc$cluster)
a=0:5
table(all.markers.dermis$cluster)
b=0:6
commongene.1=c()
for (i in a) {
              all.markers.1=all.markers.pbmc[all.markers.pbmc$cluster==i,]
              for(d in b) {
                            all.markers.2=all.markers.epi[all.markers.epi$cluster==d,]
                            all.number=length(all.markers.2$gene)+length(all.markers.1$gene)
                            commongene=length(intersect(all.markers.1$gene,all.markers.2$gene))/all.number
                            commongene.1=c(commongene.1,commongene)
              }
}

commongene.2=matrix(commongene.1,nrow = 6,ncol = 7,byrow = T)
commongene.2
colnames(commongene.2)=c(paste0('SC',0:6,'.E'))
rownames(commongene.2)=c(paste0('SC',0:5,'.P'))
commongene.2=t(commongene.2)
library(pheatmap)
p1=pheatmap(commongene.2,cluster_cols = T,cluster_rows = T,clustering_distance_rows = "correlation",scale = 'column',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 45,cellheight =25,border_color = 'white',fontsize = 25,legend_labels = 'scale.avg_logFC',legend = T)
p1
library(ggplot2)
ggsave(filename = 'T.subcluster.epi.pbmc.gene.similarity.score.scaled.scores.pdf',plot = p1,width = 40,height = 20)
write.table(commongene.2,file = 'T.epi.pbmc.commongene.similarity.score.txt',row.names = T,col.names = T,sep = '\t')
####真皮t与PBMC比较####
a=0:5
table(all.markers.dermis$cluster)
b=0:6
commongene.1=c()
for (i in a) {
              all.markers.1=all.markers.pbmc[all.markers.pbmc$cluster==i,]
              for(d in b) {
                            all.markers.2=all.markers.dermis[all.markers.dermis$cluster==d,]
                            all.number=length(all.markers.2$gene)+length(all.markers.1$gene)
                            commongene=length(intersect(all.markers.1$gene,all.markers.2$gene))/all.number
                            commongene.1=c(commongene.1,commongene)
              }
}

commongene.2=matrix(commongene.1,nrow = 6,ncol = 7,byrow = T)
commongene.2
colnames(commongene.2)=c(paste0('SC',0:6,'.D'))
rownames(commongene.2)=c(paste0('SC',0:5,'.P'))
commongene.2=t(commongene.2)
library(pheatmap)
p1=pheatmap(commongene.2,cluster_cols = T,cluster_rows = T,clustering_distance_rows = "correlation",scale = 'column',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 45,cellheight =25,border_color = 'white',fontsize = 25,legend_labels = 'scale.avg_logFC',legend = T)
p1
library(ggplot2)
ggsave(filename = 'T.subcluster.der.pbmc.gene.similarity.score.scaled.scores.pdf',plot = p1,width = 40,height = 20)
write.table(commongene.2,file = 'T.der.pbmc.commongene.similarity.score.txt',row.names = T,col.names = T,sep = '\t')
####真表皮T细胞比较####
a=0:6
table(all.markers.dermis$cluster)
b=0:6
commongene.1=c()
for (i in a) {
              all.markers.1=all.markers.epi[all.markers.epi$cluster==i,]
              for(d in b) {
                            all.markers.2=all.markers.dermis[all.markers.dermis$cluster==d,]
                            all.number=length(all.markers.2$gene)+length(all.markers.1$gene)
                            commongene=length(intersect(all.markers.1$gene,all.markers.2$gene))/all.number
                            commongene.1=c(commongene.1,commongene)
              }
}

commongene.2=matrix(commongene.1,nrow = 7,ncol = 7,byrow = T)
commongene.2
colnames(commongene.2)=c(paste0('SC',0:6,'.D'))
rownames(commongene.2)=c(paste0('SC',0:6,'.E'))
commongene.2=t(commongene.2)
library(pheatmap)
p1=pheatmap(commongene.2,cluster_cols = T,cluster_rows = T,clustering_distance_rows = "correlation",scale = 'row',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 45,cellheight =25,border_color = 'white',fontsize = 25,legend_labels = 'scale.avg_logFC',legend = T)
p1
library(ggplot2)
ggsave(filename = 'T.subcluster.der.epi.gene.similarity.score.scaled.scores.pdf',plot = p1,width = 40,height = 20)
write.table(commongene.2,file = 'T.der.epi.commongene.similarity.score.txt',row.names = T,col.names = T,sep = '\t')
p=VlnPlot(merged,features = c('COL1A1','COL3A1','CD3D','CD3E','VWF','CDH5','LYZ','AIF1','CD79A','MS4A1','TPSAB1','TPSB2','NKG7','XCL1','CDH19','MPZ'),stack = T,pt.size = 0)
p
ggsave(filename = 't.dermis.celltypemarkers.pdf',plot = p,width = 10,height = 6)
