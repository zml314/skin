####b.dermis.2022.1####
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(scater)
library(SingleR)
memory.limit(102400)
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
dir.create('./b')
setwd('./b')
####2.创建存储图形的文件夹####
sam.name <- "multi"
if(!dir.exists(sam.name)){
              dir.create(sam.name)
}

merged=B
rm(B)
merged=subset(merged,seurat_clusters!='4')
merged=subset(merged,seurat_clusters!='8')
merged=subset(merged,seurat_clusters!='5')
merged <- NormalizeData(merged, normalization.method = "LogNormalize",scale.factor = 10000)


##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
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

merged<- FindNeighbors(merged, dims = dim.use,reduction = 'harmony',nn.method = "rann")
merged <- FindClusters(merged, resolution =0.2)
merged <- RunTSNE(merged, dims = dim.use,reduction = 'harmony')
pdf(paste0("./",sam.name,"/tsne_Plot.harmony_res0.2_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = T)
dev.off()
FeaturePlot(merged,c('CD79A','MS4A1','CD3D','COL1A1'))
save(merged,file = 'b.res0.2.20220117.Rdata')
write.table(merged@meta.data,file = paste0("./",sam.name,"/",sam.name,"_cells_details_tsne_res0.2_",max(dim.use),"PC.txt"))
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
Idents(merged)=merged$RNA_snn_res.0.1
table(Idents(merged))
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.2.nolable_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res0.2_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="orig.ident", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_sample_res0.2_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="sample", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker_genes_tsne_res0.2.removed_",max(dim.use),"PC.txt"),sep="\t",quote = F)
all.markers.noRibo=all.markers[!grepl('^RP[SL]',all.markers$gene,ignore.case = F),]
all.markers.noRibo.noMito=all.markers.noRibo[!grepl('^MT-',all.markers.noRibo$gene,ignore.case = F),]
write.table(all.markers.noRibo.noMito,file = 'all.markers.b.res0.2.subcluster.noRI.moMT.txt',sep = '\t')
####12.细胞周期归类 #cc.genes为两个细胞周期的基因集####
merged<- CellCycleScoring(object = merged, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
#head(x = merged@meta.data)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_cellcycle.res0.2.removed_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 0.1)
dev.off()
####15.热图展示Top marker基因####
#筛选top5的marker基因，可以通过参数改为其他数值
marker.sig <- all.markers.noRibo.noMito[all.markers.noRibo.noMito$p_val_adj<=0.05,]

top10 <- marker.sig %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file=paste0("./",sam.name,"/",sam.name,"_top10_marker.res0.2.nomt.noribo_",max(dim.use),"PC.txt"),sep="\t",quote = F)
a=table(merged$orig.ident,merged$seurat_clusters)
a=as.data.frame(a)
write.table(a,file = 'orig.res0.2.b.subcluster.txt',sep='/t')

p1=ggplot(a,mapping = aes(Var2,Freq,fill=Var1))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='b.subcluster.orig.pdf',plot=p3,width=10,height=4)

b=table(merged$seurat_clusters,merged$sample)
b=as.data.frame(b)
write.table(b,file = 'sample.b.res0.2.sbucluster.txt',sep = '\t')
p1=ggplot(b,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='b.res0.2.subcluster.sample.pdf',plot=p3,width=6,height=4)
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
              
              write.table(enrich.go.BP,file =paste0(i,'_T.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_T.kegg.txt') ,sep='\t')
}
####找sample的marker基因####
table(Idents(merged))
Idents(merged)=merged@meta.data$sample
all.markers.sample <- FindAllMarkers(merged, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers.sample,file=paste0("./",sam.name,"/",sam.name,"_sample_marker_genes_b_",30,"PC.txt"),sep="\t",quote = F)
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
              
              write.table(enrich.go.BP,file =paste0(i,'_T.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_T.kegg.txt') ,sep='\t')
}
####marker基因展示####
Idents(merged)=merged$seurat_clusters
merged@misc$averageExpression=AverageExpression(merged)
write.table(merged@misc$averageExpression$RNA,file = 'average.expression.subb.removed.cluster.txt',sep = '\t',row.names = T,quote = F)

top10$subcluster=c(paste0('SC',top10$cluster))
table(top10$subcluster)
marker=merged@misc$averageExpression$RNA
colnames(marker)=c(paste0('SC',0:6))
a=c(paste0('SC',0:6))
for (i in a) {genelist=top10[top10$subcluster==i,]
genelist=marker[row.names(marker)%in%genelist$gene,]
write.table(genelist,file = paste0('genelist',i,'.txt'),sep = '\t')
}
markers=read.table('./top10.subcluser.b.marker.genes.txt',sep = '\t',header = T,row.names = 1)

markers=t(markers)
library(pheatmap)
p1=pheatmap(markers,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'column',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')
p1
ggsave(filename = 'b.subtype.top10.marker.dermis.expression.pdf',plot = p1,width = 45,height = 20)
####经典marker基因匹配####
Naive_b=c('IGHM','IGHD','IL4R','TCL1A')
Naive_b=marker[row.names(marker)%in%Naive_b,]
Mempry_b=c('CD27')
memory_b=marker[row.names(marker)%in%Mempry_b,]
DN2=c('TBX21','FGR', 'FCRL2', 'FCRL3', 'FCRL5','IL10RA')
dn2=marker[row.names(marker)%in%DN2,]
Activate_B=c('CD69','CD86')
Activate_B=marker[row.names(marker)%in%Activate_B,]
ABC=c('CD86','SLAMF7', 'IL2RA', 'FAS','TBK1', 'EBI3', 'ZBTB32','TRAF5')
abc=marker[row.names(marker)%in%ABC,]
Plasma=c('CD38', 'JCHAIN','MZB1','IGHG4', 'IGHA1','IGHG1','IGHG3')
Plasma=marker[row.names(marker)%in%Plasma,]
hsp=c('HSPA6','HSPA1A','HSPA1B','HSPB1','HSPH1','HSP90AA1')
hsp=marker[row.names(marker)%in%hsp,]
isg=c('ISG15','IFIT1','IFITM1','ISG20')
isg=marker[row.names(marker)%in%isg,]
classical=rbind(memory_b,abc,Plasma,hsp,isg)
colnames(classical)=c(paste0('SC',0:6))
a=c(Naive_b,Mempry_b,DN2,Activate_B,ABC,Plasma)
annotation_col=data.frame(subtype = c(rep('Memory_B',1),rep('ABC',8),rep('Plasma',7),rep('HSP_B',6),rep('ISG_B',4)),row.names = row.names(classical))
length(row.names(classical))
##辅助查找marker基因
source("./MarkerGeneList.R")
source("./function.R")
#TopMarkers_noRibo <- read.csv('CellType/TopMarkers_noRibo.csv')
MatchMarkerGene(data=all.markers, output = "./")

select.marker=select.marker[,c(2,6,7)]
select.marker<-dcast(select.marker,select.marker$gene~select.marker$cluster,value.var = 'avg_logFC',fill = 0)
row.names(select.marker)=select.marker$`select.marker$gene`
select.marker=select.marker[,-1]
annotation_col=annotation_col[annotation_col$gene%in%row.names(select.marker),]
row.names(annotation_col)=annotation_col$gene
dim(annotation_col)
annotation_col=annotation_col[,-1]
annotation_col=as.data.frame(annotation_col$subtype,row.names = annotation_col$gene)
colnames(annotation_col)='subtype'
#class(annotation_col$subtype)=as.character(annotation_col$subtype)
library(pheatmap)
row.names(classical)[1]='CD27'
row.names(annotation_col)[1]='CD27'
annotation_col$subtype
p1=pheatmap(classical,annotation_row =annotation_col ,clustering_distance_cols  = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 15,cellheight =12,border_color = 'white',fontsize = 10,cluster_rows = F)
p1
ggsave(filename = 'b.classical.subtype.marker.dermis.pdf',plot = p1,width = 6,height =10)
####气泡图####
GO=read.table('./sample.b.go.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
GO=GO[,c(3,8,10,11)]
GO1=read.table('./sample.b.go.selected.txt',sep = '\t',header = T,stringsAsFactors = F)
GO7=GO[GO$Description%in%GO1$Description,]
#GO7$subcluster=paste0('SC',GO7$subcluster)
GO7$sample=factor(GO7$sample,levels = c('HC','DLE','SLE'))
library(forcats)
GO7$Description=fct_infreq(GO7$Description)
p <- ggplot(GO7, aes(sample,Description, size = Count, color=qvalue)) + geom_point() +xlab("sample")
#### 气泡图颜色####
p2=p + scale_color_gradient(low='red',high='blue') +theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.major=element_line(colour=NA))+theme(text = element_text(size = 20))
p2
ggsave(filename = 'burbble.pathway.subb.dermis.20220213.pdf',plot = p2,width =10,height = 10)
####与经典marker匹配####
classical.marker=read.table('./classical.gene.b.txt',header = T,sep = '\t',row.names = 1)
table(classical.marker$subtype)
marker=read.table('./average.expression.subb.removed.cluster.txt',sep = '\t',header = T,row.names = 1)
colnames(marker)=c(paste0('SC',0:6))
#abs=classical.marker[classical.marker$subtype=='ABC',]
#abs=marker[row.names(marker)%in%row.names(abs),]
ABC=c('CD86','SLAMF7', 'IL2RA', 'FAS','TBK1', 'EBI3', 'ZBTB32','TRAF5')
abc=marker[row.names(marker)%in%ABC,]

DN2=classical.marker[classical.marker$subtype=='DN2',]
dn2=marker[row.names(marker)%in%row.names(DN2),]
HSP=classical.marker[classical.marker$subtype=='HSP_B',]
hsp=marker[row.names(marker)%in%row.names(HSP),]
ISG=classical.marker[classical.marker$subtype=='ISG_B',]
isg=marker[row.names(marker)%in%row.names(ISG),]
Memory=classical.marker[classical.marker$subtype=='Memory_B',]
mem=marker[row.names(marker)%in%row.names(Memory),]
Naive=classical.marker[classical.marker$subtype=='Naive_B',]
Naive=marker[row.names(marker)%in%row.names(Naive),]
plasma=classical.marker[classical.marker$subtype=='Plasma',]
plasma=marker[row.names(marker)%in%row.names(plasma),]
gene=rbind(abc,dn2,hsp,isg,mem,Naive,plasma)
classical.marker=data.frame(subtype=classical.marker$subtype,row.names = row.names(classical.marker))
write.table(gene,file = 'classical.gene.subcluster.b.dermis.txt',sep = '\t')
row.names(classical.marker)
classsical.marker.1=classical.marker[!classical.marker$subtype=='ABC',]
genes=row.names(classical.marker)[classical.marker$subtype!='ABC']
classical.marker=data.frame(subtype=classsical.marker.1,row.names = genes)
abc=data.frame(subtype=rep('ABC',length(ABC)),row.names = row.names(abc))
classical.marker=rbind(classical.marker,abc)
gene=read.table('classical.gene.subcluster.b.epidermis.txt',sep = '\t',header = T,row.names = 1)
library(pheatmap)
p1=pheatmap(gene,cluster_cols = F,annotation_row =classical.marker ,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')
p1
ggsave(filename = 'B.classical.gene.dermis.expression.2022.pdf',plot = p1,width = 45,height = 20)
gene=t(gene)
p1=pheatmap(gene,cluster_cols = F,annotation_col =classical.marker ,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'column',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white',angle_col = 45)
p1
ggsave(filename = 'B.classical.gene.dermis.expression.pdf',plot = p1,width = 45,height = 20)
####真表皮B与PBMC中的B比较####
####表皮的B与PBMC中的B比较####
getwd()
table(all.markers$cluster)
table(Idents(merged))
all.markers.pbmc=read.table('./sc.markergenes.pbmc.txt',sep = '\t',header = T)
all.markers.dermis=read.table('./multi/multi_marker_genes_tsne_res0.2.removed_20PC.txt',sep = '\t',header = T)
table(all.markers.dermis$cluster)
table(merged$seurat_clusters)
all.markers.epi=read.table('./multi/epi.b.multi_marker_genes_tsne_res0.2.removed_15PC.txt',header = T,row.names = 1)
table(all.markers.epi$cluster)
colnames(all.markers.pbmc)=c('genes','cluster')
table(all.markers.pbmc$cluster)
a=c(paste0('B_SC',0:6),paste0('P_SC',0:1))

b=0:4
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

commongene.2=matrix(commongene.1,nrow = 9,ncol = 5,byrow = T)
colnames(commongene.2)=c(paste0('SC',0:4))
rownames(commongene.2)=c(paste0('B_SC',0:6),paste0('P_SC',0:1))
anno.row=data.frame(subtype.2=c(rep('subclusters of B cell in PBMC',7),rep('subclusters of plasma in PBMC',2)),row.names = c(paste0('B_SC',0:6),paste0('P_SC',0:1)))
library(pheatmap)
p1=pheatmap(commongene.2,annotation_row = anno.row,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 45,cellheight =25,border_color = 'white',fontsize = 25,legend = T,angle_col = 45)
p1
library(ggplot2)
ggsave(filename = 'b.subcluster.epider.capbmc.gene.similarity.score.scaled.scores.pdf',plot = p1,width = 40,height = 20)
write.table(commongene.2,file = 'b.epi.pbmc.commongene.similarity.score.txt',row.names = T,col.names = T,sep = '\t')
####真皮t与PBMC比较####
a=c(paste0('B_SC',0:6),paste0('P_SC',0:1))

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

commongene.2=matrix(commongene.1,nrow = 9,ncol = 7,byrow = T)
commongene.2
colnames(commongene.2)=c(paste0('SC',0:6))
rownames(commongene.2)=c(paste0('B_SC',0:6),paste0('P_SC',0:1))

commongene.2=t(commongene.2)
library(pheatmap)
p1=pheatmap(commongene.2,cluster_cols = T,cluster_rows = T,clustering_distance_rows = "correlation",scale = 'row',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 45,cellheight =25,border_color = 'white',fontsize = 25,legend_labels = 'scale.avg_logFC',legend = T)
p1
library(ggplot2)
ggsave(filename = 'b.subcluster.der.pbmc.gene.similarity.score.scaled.scores.pdf',plot = p1,width = 40,height = 20)
write.table(commongene.2,file = 'b.der.pbmc.commongene.similarity.score.txt',row.names = T,col.names = T,sep = '\t')
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
ggsave(filename = 'b.dermis.celltypemarkers.pdf',plot = p,width = 10,height = 6)
####
classical.marker=read.table('./classical.gene.b.txt',header = T,sep = '\t',row.names = 1)
colnames(classical.marker)
classical.marker=data.frame(subtype=classical.marker$subtype,row.names = row.names(classical.marker))
colnames(marker)=c(paste0('SC',0:6))
ABC=c('ITGAX','TBX21','FCR','FCRL2','EBI3','CD86')
abc=marker[row.names(marker)%in%ABC,]
HSP=row.names(classical.marker)[classical.marker$subtype=='HSP_B']
hsp=marker[row.names(marker)%in%HSP,]
ISG=row.names(classical.marker)[classical.marker$subtype=='ISG_B']
isg=marker[row.names(marker)%in%ISG,]
Memory=row.names(classical.marker)[classical.marker$subtype=='Memory_B']
mem=marker[row.names(marker)%in%Memory,]
Naive=row.names(classical.marker)[classical.marker$subtype=='Naive_B']
Naive=marker[row.names(marker)%in%Naive,]
plasma=row.names(classical.marker)[classical.marker$subtype=='Plasma']
plasma.1=marker[row.names(marker)%in%plasma,]
genes=rbind(abc,hsp,isg,mem,Naive,plasma.1)
anno.row=data.frame(subtype=c(rep('ABC',5),rep('HSP',7),
                              rep('ISG',4),rep('Memory',2),
                              rep('Naive',5),rep('Plasma',8)),
                    row.names = row.names(genes))
colnames(genes.1)=c(paste0('B_',colnames(genes.1)))
library(pheatmap)
p1=pheatmap(genes.1,cluster_cols = F,annotation_row =anno.row,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')+NoLegend()
p1
genes.2
row.names(genes.2)
genes.1=genes[c(1:5),]
genes.2=genes[-c(1:5),]
genes.2=genes.2[-c(13),]
