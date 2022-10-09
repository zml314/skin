####mac.dermis.2022####
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(scater)
library(SingleR)
library(monocle)
library(tibble)
library(harmony)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)
library(DOSE)
library(devtools)
setwd('E:/R_analysis/dermis.2022.1')
dir.create('./mac')
setwd('./mac/')
sam.name <- "multi"
if(!dir.exists(sam.name)){
              dir.create(sam.name)
}
merged=macro
rm(macro)
merged=subset(merged,seurat_clusters!='4')
merged=subset(merged,seurat_clusters!='6')
merged=subset(merged,seurat_clusters!='7')
merged=subset(merged,seurat_clusters!='9')
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
all.genes <- rownames(merged)
merged <- ScaleData(merged, features = row.names(merged))#需要足够大的内存
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
library(harmony)

merged=RunHarmony(merged,group.by.vars ='orig.ident')
pdf(paste0("./",sam.name,"/PCA-DimPlot.harmony.pdf"),width = 5,height = 4)
DimPlot(merged, reduction = "harmony",group.by = 'orig.ident')
dev.off()
#TSNE算法
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))

merged<- FindNeighbors(merged, dims = dim.use,reduction = 'harmony')
merged <- FindClusters(merged, resolution = 0.1)
table(Idents(merged))
Idents(merged)=merged$seurat_clusters
merged <- RunTSNE(merged, dims = dim.use,check_duplicates = FALSE,reduction = 'harmony')
FeaturePlot(merged,features = c('LYZ','AIF1','CD3D','COL1A1','KRT1'))
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.1_",max(dim.use),"PC.pdf"),width = 6,height = 5)
DimPlot(object = merged, pt.size=0.1,label = T,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.1.nolabel_",max(dim.use),"PC.pdf"),width = 6,height = 5)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols)
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res0.1_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="orig.ident", pt.size=0.1,reduction = "tsne",label = F)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_sample_res0.1_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="sample", pt.size=0.1,reduction = "tsne",label = F)
dev.off()
####12.细胞周期归类 #cc.genes为两个细胞周期的基因集####
merged<- CellCycleScoring(object = merged, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
#head(x = merged@meta.data)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_cellcycle_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(merged,reduction = "tsne",label = F,group.by="Phase",pt.size = 0.1)
dev.off()
#### 13. 计算marker基因 ####
all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker_genes_tsne_res0.3_",max(dim.use),"PC.txt"),sep="\t",quote = F)
all.markers.noRibo=all.markers[!grepl('^RP[SL]',all.markers$gene,ignore.case = F),]
all.markers.noRibo.noMito=all.markers.noRibo[!grepl('^MT-',all.markers.noRibo$gene,ignore.case = F),]
write.table(all.markers.noRibo.noMito,file = 'all.markers.b.res0.1.subcluster.noRI.moMT.txt',sep = '\t')
Idents(merged)=merged$sample
all.markers.sample <- FindAllMarkers(merged, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers.sample,file=paste0("./",sam.name,"/",sam.name,"_marker_genes.sample_",max(dim.use),"PC.txt"),sep="\t",quote = F)

pdf(file = 'macro.dc.marker.1.pdf',width = 15,height = 12)
FeaturePlot(merged, features = c('LYZ','TCF4','FCGR3A','CD14','CD1C','IL1B','AIF1','ITGAE'),cols = c("gray", "purple"),min.cutoff = 0)
dev.off()
pdf(file = 'CXCR3,CCR6.pdf',width = 10,height = 4)
FeaturePlot(merged, features = c('CXCR3','CCR6'),cols = c("gray", "red"),min.cutoff = 0)
dev.off()
View(all.markers)
b=table(merged$sample,merged$seurat_clusters)
b=as.data.frame(b)
write.table(b,file = 'sample.subcluster.txt',row.names = F,sep = '\t')
c=table(merged$orig.ident,merged$seurat_clusters)

c=as.data.frame(c)
write.table(c,file = 'orig.subcluster.txt',row.names = F,sep = '\t')


#d$Var1=factor(d$Var1,levels = c('HC_1','HC_2','HC_3','DLE_1','DLE_2','SLE_1','SLE_2','SLE_3'))
p1=ggplot(b,mapping = aes(Var2,Freq,fill=Var1))+geom_bar(stat='identity',position='fill') +
              labs(x = 'subclusters',y = 'Frequency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                        axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='subcluster.sample.fill.pdf',plot=p3,width=6,height=5)
p1=ggplot(b,mapping = aes(Var2,Freq,fill=Var1))+geom_bar(stat='identity',position='fill') +
              labs(x = 'Subcluster',y = 'Frequency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                            axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='subcluster.sample.fill.pdf',plot=p3,width=6,height=5)
p1=ggplot(c,mapping = aes(Var2,Freq,fill=Var1))+geom_bar(stat='identity',position='fill') +
              labs(x = 'Subcluster',y = 'Frequency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                            axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='subcluster.orig.fill.pdf',plot=p3,width=10,height=5)
p1=ggplot(c,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'Sample',y = 'Frequency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                        axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='orig.subcluster.fill.pdf',plot=p3,width=6,height=5)

marker.sig <- all.markers.noRibo.noMito[all.markers.noRibo.noMito$p_val_adj<=0.05,]

top10 <- marker.sig %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file=paste0("./",sam.name,"/",sam.name,"_top10_marker.res0.2.nomt.noribo_",max(dim.use),"PC.txt"),sep="\t",quote = F)
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
              
              write.table(enrich.go.BP,file =paste0(i,'_mac.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_mac.kegg.txt') ,sep='\t')
}
#### 13. 计算sample.marker基因 ####
all.markers.sample <- FindAllMarkers(merged, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers.sample,file=paste0("./",sam.name,"/",sam.name,"_sample.marker_genes_tsne_res0.3_",max(dim.use),"PC.txt"),sep="\t",quote = F)


a=c('HC','DLE','SLE')
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
              
              write.table(enrich.go.BP,file =paste0(i,'_M.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_M.kegg.txt') ,sep='\t')
}
####marker基因展示####
Idents(merged)=merged$seurat_clusters
merged@misc$averageExpression=AverageExpression(merged)
write.table(merged@misc$averageExpression$RNA,file = 'average.expression.submac.removed.cluster.txt',sep = '\t',row.names = T,quote = F)

top10$subcluster=c(paste0('SC',top10$cluster))
table(top10$subcluster)
marker=merged@misc$averageExpression$RNA
colnames(marker)=c(paste0('SC',0:6))
a=c(paste0('SC',0:6))
for (i in a) {genelist=top10[top10$subcluster==i,]
genelist=marker[row.names(marker)%in%genelist$gene,]
write.table(genelist,file = paste0('genelist',i,'.txt'),sep = '\t')
}
markers=read.table('./mac.subcluster.top10.gene.expression.txt',sep = '\t',header = T,row.names = 1)

markers=t(markers)
library(pheatmap)
p1=pheatmap(markers,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'column',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')
p1
ggsave(filename = 'mac.subtype.top10.marker.dermis.expression.pdf',plot = p1,width = 45,height = 20)
DotPlot(merged, features = top10$gene)+RotatedAxis()+
              scale_x_discrete("")+scale_y_discrete("")+theme(panel.grid.major=element_line(colour=NA))+
              theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
marker=read.table('./average.expression.submac.removed.cluster.txt',sep = '\t',header = T,row.names = 1)
IL1b_Macro=c('IL1B','CD68')
IL1b_Macro=marker[row.names(marker)%in%IL1b_Macro,]
HSP=c('HSPH1','HSPA1A','HSP90AA1','HSPA6','HSPE1')
HSP=marker[row.names(marker)%in%HSP,]
ISG=c('STAT1','IFI44','IFIH1','ISG15','IFI6','DDX58','IFIT3','IFIT1')
ISG=marker[row.names(marker)%in%ISG,]
macro=c('CD14','CD68','FCGR3A','RNASE1')
CD1C=c('CD1C','CD1A','PKIB','FCGR2B','CLEC10A')
pDC=c('GZMB','JCHAIN','IRF4','IRF7','TCF4','IL3RA','MZB1')
CLEC9A=c('CLEC9A')
langerhans=c('CD207')
macro=marker[row.names(marker)%in%macro,]
cd1c=marker[row.names(marker)%in%CD1C,]
pdc=marker[row.names(marker)%in%pDC,]
clec9a=marker[row.names(marker)%in%CLEC9A,]
langerhans=marker[row.names(marker)%in%langerhans,]
classical=rbind(macro,cd1c,pdc,clec9a,IL1b_Macro,HSP,ISG,langerhans)
row.names(classical)
row.names(classical)[17]='CLEC9A'
anno.col=data.frame(cell_type=c('Macrophage','Macrophage','DC','DC','DC','DC','Macrophage'),row.names = c(paste0('SC',0:6)))
anno.row=data.frame(subtype=c('CD16_Macrophage','CD14_Macrophage','cDC','pDC','pDC','pDC','pDC','pDC','pDC','pDC',
                              'CLEC9A_DC','HSP','HSP','HSP','HSP','HSP','ISG','ISG','ISG','ISG','ISG','ISG','ISG','ISG','Langerhans'),
                    row.names = c('FCGR3A','CD14','CD1C','JCHAIN','MZB1','IRF4','IRF7','GZMB','TCF4','IL3RA','CLEC9A','HSPA6',
                                  'HSPE1','HSPA1A','HSPH1','HSP90AA1','ISG15','IFI6','IFI44','IFIH1','STAT1','DDX58','IFIT3','IFIT1','CD207'))
library(pheatmap)
colnames(classical)=c(paste0('SC',0:5))
p1=pheatmap(classical,annotation_row = anno.row,annotation_col = anno.col,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 30,cellheight =25,fontsize = 22,legend = T,border_color = 'white')
p1
ggsave(filename = 'classicla.mac.subtype.top10.genes.epidermis.expression.pdf',plot = p1,width = 45,height = 20)
save(merged,file = 'mac.merged.res0.1.Rdata')
####气泡图####
GO=read.table('./mac.sample.all.go.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
GO=GO[,c(3,8,10,11)]
GO1=read.table('./mac.sample.selected.go.txt',sep = '\t',header = T,stringsAsFactors = F)
GO7=GO[GO$Description%in%GO1$Description,]
GO7$subcluster=paste0('SC',GO7$subcluster)
library(forcats)
library(ggplot2)
GO7$Description=fct_infreq(GO7$Description)
GO7$sample=factor(GO7$sample,levels = c('HC','DLE','SLE'))
p <- ggplot(GO7, aes(sample,Description, size = Count, color=qvalue)) + geom_point() +xlab("sample")
#### 气泡图颜色####
p2=p + scale_color_gradient(low='red',high='blue') +theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.major=element_line(colour=NA))+theme(text = element_text(size = 20))
p2
ggsave(filename = 'burbble.pathway.sample.mac.20220213.pdf',plot = p2,width =10,height = 10)
####与经典marker比较####
IL1b_Macro=c('IL1B','CD68')
IL1b_Macro=marker[row.names(marker)%in%IL1b_Macro,]
HSP=c('HSPH1','HSPA1A','HSP90AA1','HSPA6','HSPE1')
HSP=marker[row.names(marker)%in%HSP,]
ISG=c('STAT1','IFI44','IFIH1','ISG15','IFI6','DDX58','IFIT3','IFIT1')
ISG=marker[row.names(marker)%in%ISG,]
macro=c('CD14','CD68','FCGR3A','RNASE1')
CD1C=c('CD1C','CD1A','PKIB','FCGR2B','CLEC10A')
pDC=c('GZMB','JCHAIN','IRF4','IRF7','TCF4','IL3RA','MZB1')
langer=c('CD207')
CLEC9A=c('CLEC9A')

macro=marker[row.names(marker)%in%macro,]
cd1c=marker[row.names(marker)%in%CD1C,]
pdc=marker[row.names(marker)%in%pDC,]
langer=marker[row.names(marker)%in%langer,]
clec9a=marker[row.names(marker)%in%CLEC9A,]
classical=rbind(macro,cd1c,pdc,langer,clec9a,IL1b_Macro,HSP,ISG)
classical=t(classical)
row.names(classical)=c(paste0('SC',0:5))
anno.col=data.frame(cell_type=c('Macrophage','DC','DC','Mix','DC','DC'),row.names = colnames(classical))
p1=pheatmap(classical,annotation_col = anno.col,cluster_cols = T,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white',angle_col = 45)
p1
ggsave(filename = 'classicla.subtype.top10.genes.epidermis.expression.pdf',plot = p1,width = 45,height = 20)
p=VlnPlot(merged,features = c('COL1A1','COL3A1','CD3D','CD3E','VWF','CDH5','LYZ','AIF1','CD79A','MS4A1','TPSAB1','TPSB2','NKG7','XCL1','CDH19','MPZ'),stack = T,pt.size = 0)
p
ggsave(filename = 'mac.dermis.celltype.genes.pdf',plot = p,width = 10,height = 6)
####macro/dc打分####
####基因集打分####
####巨噬细胞得分箱线图####
kegggmt2 <- read.gmt("c8.all.v7.5.1.symbols.gmt")
genelist=kegggmt2$gene[kegggmt2$term=='TRAVAGLINI_LUNG_MACROPHAGE_CELL']
genelist=kegggmt2$gene[kegggmt2$term=='HAY_BONE_MARROW_DENDRITIC_CELL']
genelist=c(genelist,'CLEC9A','JCHAIN','MZB1','IRF4','IRF7','GZMB','TCF4','IL3RA','CD207')
a=merged@assays$RNA@var.features
b=genelist%in%a
table(b)
b=genelist[b]
b
b=list(b)
data <- FetchData(merged,vars = c("subtype","sample","chemokines1"))
merged <- Seurat::AddModuleScore(object = merged,features = b,ctrl = 200,name = "dc")
colnames(merged@meta.data)

VlnPlot(merged,features = "dc1",pt.size = 0)
data <- FetchData(merged,vars = c("seurat_clusters","macrophage1"))
#colnames(data)[3]='IRGs'
data=data[data$subtype%in%c('Langerhans'),]
p <- ggplot(data = data,aes(seurat_clusters,chemokines,fill = sample))
p1=p + geom_boxplot() + RotatedAxis() + labs(title = "Chemokines and Receptors  Score",y = "Score") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme_bw()+theme(panel.grid=element_blank()) 
p2=p1+theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1))

p3=p2+stat_compare_means(label = "p.format")###添加统计学分析P值
p3
ggsave(filename = 'IRGs score.fib.subtype.sample.xiangxiantu.pdf',plot = p2,width = 8,height = 5)
####与STM文章对比marker####
cDC1=c('CLEC9A','WDFY4','CCND1','DNASE1L3','CPNE3','TCEA3','CPVL')
cDC1=marker[row.names(marker)%in%cDC1,]
cDC2A=c('LAMP3','EBI3','CCL19','FSCN1','CCL17','CCL22','BIRC3')
cDC2A=marker[row.names(marker)%in%cDC2A,]
cDC2B=c('CD1C','CLEC10A','FCER1A','CXCL8','IL1B','G0S2','CXCL3')
cDC2B=marker[row.names(marker)%in%cDC2B,]
pDClike=c('PPP1R14A','TRPM4','AXL','KLF4','LILRA4','AREG','PTGDS')
pDClike=marker[row.names(marker)%in%pDClike,]
pDC=c('GZMB','IGKC','JCHAIN','PLAC8','TCF4','IL3RA','CLEC4C')
pDC=marker[row.names(marker)%in%pDC,]
CD16DC=c('CXCL11','S100A9','APOBEC3A','CXCL10','CCL3','S100A8','CXCL9','FCGR3A')
CD16DC=marker[row.names(marker)%in%CD16DC,]
#LC=c('CD207','FCGBP','CD1A','HLA-DQB2','TACSTD2','S100B','ACOT7')
#LC=marker[row.names(marker)%in%LC,]
LAM=c('SPP1','APOE','APOC1','CTSL','CCL2','CTSD','FABP4')
LAM=marker[row.names(marker)%in%LAM,]
PVM=c('SELENOP','RNASE1','STAB1','MAF','C1QA','F13A1','CCL18')
PVM=marker[row.names(marker)%in%PVM,]
MAC=c('MRC1','CD163','CD68','CCR1','MACRO','TREM2')
MAC=marker[row.names(marker)%in%MAC,]
ISG=c('STAT1','IFI44','IFIH1','ISG15','IFI6','DDX58','IFIT3','IFIT1')
ISG=marker[row.names(marker)%in%ISG,]
HSP=c('HSPA6','HSPA1A','HSPA1B','HSPB1','HSPA5','HSPH1','HSP90AA1')
HSP=marker[row.names(marker)%in%HSP,]
genes=rbind(cDC1,cDC2A,cDC2B,pDClike,pDC,CD16DC,LAM,PVM,MAC,ISG,HSP)
anno.row=data.frame(subtype=c(rep('cDC1',7),rep('cDC2A',7),rep('cDC2B',7),
                              rep('pDClike',7),rep('pDC',7),rep('CD16DC',8),
                              rep('LAM',7),rep('PVM',7),rep('MAC',5),
                              rep('ISG',8),rep('HSP',7)),
                    row.names = row.names(genes))
library(pheatmap)
#genes=t(genes)
anno.col=data.frame(cell_type=c('Macro','DC','Mix','DC','DC','DC','DC'),row.names = colnames(genes))
p1=pheatmap(genes,annotation_col = anno.col,annotation_row = anno.row,cluster_cols = T,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 30,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')
p1
b=subset(merged,seurat_clusters=='2')
save(b,file = 'dermis.mix.mac.dc.Rdata')
####mix再聚类####
dim(b)
merged=b
merged <- NormalizeData(merged, normalization.method = "LogNormalize",scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged)
merged <- ScaleData(merged, features = row.names(merged))#需要足够大的内存
merged <- ScaleData(merged, vars.to.regress = "percent.mt")
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
ElbowPlot(merged,ndims = 40)
dim.use=1:20
library(harmony)
merged=RunHarmony(merged,group.by.vars ='orig.ident')
#TSNE算法
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
merged<- FindNeighbors(merged, dims = dim.use,reduction = 'harmony')
merged <- FindClusters(merged, resolution = 0.1)
table(Idents(merged))
Idents(merged)=merged$seurat_clusters
merged <- RunTSNE(merged, dims = dim.use,check_duplicates = FALSE,reduction = 'harmony')

CD16DC=c('FCGR3A','CXCL11','S100A9','APOBEC3A','CXCL10','CCL3','S100A8','CXCL9')
LAM=c('SPP1','APOE','APOC1','CTSL','CCL2','CTSD','FABP4')
LAM=marker[row.names(marker)%in%LAM,]
PVM=c('SELENOP','RNASE1','STAB1','MAF','C1QA','F13A1','CCL18')
PVM=marker[row.names(marker)%in%PVM,]
MAC=c('MRC1','CD163','CD68','CCR1','MACRO','TREM2')
MAC=marker[row.names(marker)%in%MAC,]
genes.1=c(CD16DC,LAM,PVM,MAC)
VlnPlot(merged,features = genes.1,stack = T,flip = T)+NoLegend()
merged@misc$averageExpression=AverageExpression(merged)
marker=merged@misc$averageExpression$RNA
CD16DC=marker[row.names(marker)%in%CD16DC,]
classical=rbind(CD16DC,LAM,PVM,MAC)
anno.row=data.frame(classical=c(rep('CD16 DC',8),rep('LAM',7),rep('PVM',7),rep('MAC',5)),row.names = row.names(classical))
p1=pheatmap(classical,annotation_row = anno.row,cluster_cols = T,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white',angle_col = 45)
p1

all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
marker.sig <- all.markers[all.markers$p_val<=0.05,]

top10 <- marker.sig %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DotPlot(merged,features = top10$gene)+RotatedAxis()+
              scale_x_discrete("")+scale_y_discrete("")
pdf(paste0("./",sam.name,"/Mix.CellCluster-TSNEPlot_res0.1_",max(dim.use),"PC.pdf"),width = 6,height = 5)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols)
dev.off()
save(merged,file = 'mix.res0.1.Rdata')
