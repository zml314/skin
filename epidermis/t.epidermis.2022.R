####t.epidermis.2022####
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
setwd('E:/R_analysis/epidermis.2022.1/')
dir.create('./t')
setwd('./t')
####2.创建存储图形的文件夹####
sam.name <- "multi"
if(!dir.exists(sam.name)){
              dir.create(sam.name)
}
merged=T
rm(T)
merged=subset(merged,seurat_clusters!='4')
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
table(merged$seurat_clusters)
merged=subset(merged,seurat_clusters!=6)
merged<- FindNeighbors(merged, dims = dim.use,reduction = 'harmony',nn.method = "rann")
merged <- FindClusters(merged, resolution =0.15)
merged <- RunTSNE(merged, dims = dim.use,reduction = 'harmony')

pdf(paste0("./",sam.name,"/tsne_Plot.harmony_res0.2_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = T)
dev.off()
FeaturePlot(merged,features = c('CD3D','CD3E','CD4','CD8A','KRT1'))
save(merged,file = 't.res0.15.epidermis.20220207.Rdata')
write.table(merged@meta.data,file = paste0("./",sam.name,"/",sam.name,"_cells_details_tsne_res0.15_",max(dim.use),"PC.txt"))

#PCA结果展示-2
pdf(paste0("./",sam.name,"/t_epidermis.res0.15.pdf"),width = 5,height = 4)
DimPlot(merged, reduction = "pca",group.by = 'orig.ident')
dev.off()
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.15_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = T,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.15.nolegend_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res0.15_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="orig.ident", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_sample_res0.3_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="sample", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker_genes_tsne_res0.15.removed_",max(dim.use),"PC.txt"),sep="\t",quote = F)
all.markers.noRibo=all.markers[!grepl('^RP[SL]',all.markers$gene,ignore.case = F),]
all.markers.noRibo.noMito=all.markers.noRibo[!grepl('^MT-',all.markers.noRibo$gene,ignore.case = F),]
write.table(all.markers.noRibo.noMito,file = 'all.markers.t.subcluster.noRI.moMT.txt',sep = '\t')
####12.细胞周期归类 #cc.genes为两个细胞周期的基因集####
merged<- CellCycleScoring(object = merged, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
#head(x = merged@meta.data)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_cellcycle.res0.3removed_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 0.1)
dev.off()
####15.热图展示Top marker基因####
#筛选top5的marker基因，可以通过参数改为其他数值
marker.sig <- all.markers[all.markers$p_val_adj<=0.05,]

top10 <- marker.sig %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file=paste0("./",sam.name,"/",sam.name,"_top10_marker.res0.15_",max(dim.use),"PC.txt"),sep="\t",quote = F)
a=table(merged$orig.ident,merged$seurat_clusters)
a=as.data.frame(a)
write.table(a,file = 'orig.res0.2.t.subcluster.txt',sep='/t')

p1=ggplot(a,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='t.orig.res0.2.subcluster.pdf',plot=p3,width=10,height=4)

b=table(merged$seurat_clusters,merged$sample)
b=as.data.frame(b)
write.table(b,file = 'sample.t.res0.2.sbucluster.txt',sep = '\t')
p1=ggplot(b,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='t.res0.2.subcluster.sample.pdf',plot=p3,width=5,height=4)
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
####找sample的marker基因####
table(Idents(merged))
Idents(merged)=merged@meta.data$sample
all.markers.sample <- FindAllMarkers(merged, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers.sample,file=paste0("./",sam.name,"/",sam.name,"_sample_marker_genes_b_",30,"PC.txt"),sep="\t",quote = F)
####marker基因展示####
table(Idents(merged))
Idents(merged)=merged$seurat_clusters
merged@misc$averageExpression=AverageExpression(merged)
write.table(merged@misc$averageExpression$RNA,file = 'average.expression.subb.removed.cluster.txt',sep = '\t',row.names = T,quote = F)
Idents(merged)=merged$sample
merged@misc$averageExpression=AverageExpression(merged)
write.table(merged@misc$averageExpression$RNA,file = 'average.expression.subt.removed.sample.txt',sep = '\t',row.names = T,quote = F)

top10$subcluster=c(paste0('SC',top10$cluster))
table(top10$subcluster)
marker=merged@misc$averageExpression$RNA
colnames(marker)=c(paste0('SC',0:6))
a=c(paste0('SC',0:6))
for (i in a) {genelist=top10[top10$subcluster==i,]
genelist=marker[row.names(marker)%in%genelist$gene,]
write.table(genelist,file = paste0('genelist',i,'.txt'),sep = '\t')
}
markers=read.table('./top10.t.epidermis.average.expression.txt',sep = '\t',header = T,row.names = 1)
#row.names(markers)=markers$X
markers=t(markers)
library(pheatmap)
p1=pheatmap(markers,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'column',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')
p1
ggsave(filename = 't.subtype.top10.marker.epidermis.expression.pdf',plot = p1,width = 45,height = 20)
####与经典基因匹配####
marker=merged@misc$averageExpression$RNA

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
####气泡图####
GO=read.table('./t.sample.go.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
GO=GO[,c(3,8,10,11)]
GO1=read.table('./t.sample.go.selected.txt',sep = '\t',header = T,stringsAsFactors = F)
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
ggsave(filename = 'burbble.pathway.sample.t.20220214.pdf',plot = p2,width =10,height = 10)
VlnPlot(merged,features = c('KRT14','KRT1','CD3D','CD3E','LYZ','AIF1','PMEL','MLANA','NKG7','XCL2','CD79A','MS4A1'),stack = T)
merged$subcluster=c(paste0('T_SC',merged$seurat_clusters))
table(merged$subcluster)
Idents(merged)=merged$subcluster
FeaturePlot(merged,features =c('KRT14','KRT1','CD3D','CD3E','LYZ','AIF1','PMEL','MLANA','NKG7','XCL2','CD79A','MS4A1') )
p=VlnPlot(object = merged,features = c('KRT14','KRT1','CD3D','CD3G','LYZ','AIF1','PMEL','MLANA','NKG7','XCL2','CD79A','MS4A1'),stack = T,pt.size = 0)
p
ggsave(filename = 't.epidermis.celltype.genes.pdf',plot = p,width = 8,height = 5)
pdf(paste0("./",sam.name,"基因线粒体比例核糖体基因比例红细胞比例.pdf"),width = 15,height = 5)
VlnPlot(merged, features = c("nFeature_RNA", "percent.mt","percent.redcell"), ncol = 4,group.by = "orig.ident")
dev.off()