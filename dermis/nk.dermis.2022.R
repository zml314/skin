####nk.dermis.2022####
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
setwd('E:/R_analysis/dermis.2022.1/')
dir.create('./nk')
setwd('./nk')
####2.创建存储图形的文件夹####
sam.name <- "multi"
if(!dir.exists(sam.name)){
              dir.create(sam.name)
}

merged=nk
rm(nk)
merged=subset(merged,seurat_clusters!='2')
merged=subset(merged,seurat_clusters!='6')
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
FeaturePlot(merged,c('GNLY','XCL1','NKG7','XCL2','CD3D','COL1A1'))
save(merged,file = 'b.res0.2.20220117.Rdata')
write.table(merged@meta.data,file = paste0("./",sam.name,"/",sam.name,"_cells_details_tsne_res0.2_",max(dim.use),"PC.txt"))
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.2_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols)
dev.off()

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res0.2_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="orig.ident", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_SamGroup_sample_res0.2_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, group.by="sample", pt.size=0.1,reduction = "tsne",label = F,cols = cols)
dev.off()
Idents(merged)=merged$sample
all.markers <- FindAllMarkers(merged, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
all.markers$gene[all.markers$cluster=='SLE']
VlnPlot(merged,features = c('IFI27','IFITM3','IFI6','ISG15','ISG20'),pt.size = 0,stack = T)
write.table(all.markers,file=paste0("./",sam.name,"/",sam.name,"_marker_genes_tsne_res0.2.removed_",max(dim.use),"PC.txt"),sep="\t",quote = F)
all.markers.noRibo=all.markers[!grepl('^RP[SL]',all.markers$gene,ignore.case = F),]
all.markers.noRibo.noMito=all.markers.noRibo[!grepl('^MT-',all.markers.noRibo$gene,ignore.case = F),]
write.table(all.markers.noRibo.noMito,file = 'all.markers.NK.res0.2.subcluster.noRI.moMT.txt',sep = '\t')
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
write.table(a,file = 'orig.res0.2.NK.subcluster.txt',sep='/t')
p1=ggplot(a,mapping = aes(Var2,Freq,fill=Var1))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='NK.subcluster.orig.pdf',plot=p3,width=8,height=4)

b=table(merged$seurat_clusters,merged$sample)
b=as.data.frame(b)
write.table(b,file = 'sample.NK.res0.2.sbucluster.txt',sep = '\t')
p1=ggplot(b,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'subcluster',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='nk.res0.2.subcluster.sample.pdf',plot=p3,width=5,height=4)
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
a=0:4
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
              
              write.table(enrich.go.BP,file =paste0(i,'_nk.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_nk.kegg.txt') ,sep='\t')
}
####找sample的marker基因####
table(Idents(merged))
Idents(merged)=merged@meta.data$sample
all.markers.sample <- FindAllMarkers(merged, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers.sample,file=paste0("./",sam.name,"/",sam.name,"_sample_marker_genes_nk_",30,"PC.txt"),sep="\t",quote = F)
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
              
              write.table(enrich.go.BP,file =paste0(i,'_nk.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_nk.kegg.txt') ,sep='\t')
}
####marker基因展示####
Idents(merged)=merged$seurat_clusters
merged@misc$averageExpression=AverageExpression(merged)
write.table(merged@misc$averageExpression$RNA,file = 'average.expression.subb.removed.cluster.txt',sep = '\t',row.names = T,quote = F)

top10$subcluster=c(paste0('SC',top10$cluster))
table(top10$subcluster)
marker=merged@misc$averageExpression$RNA
colnames(marker)=c(paste0('SC',0:4))
a=c(paste0('SC',0:4))
for (i in a) {genelist=top10[top10$subcluster==i,]
genelist=marker[row.names(marker)%in%genelist$gene,]
write.table(genelist,file = paste0('genelist',i,'.txt'),sep = '\t')
}
markers=read.table('./nk.subcluster.top10.txt',sep = '\t',header = T,row.names = 1)

markers=t(markers)
library(pheatmap)
p1=pheatmap(markers,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'column',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')
p1
ggsave(filename = 'nk.subtype.top10.marker.dermis.expression.pdf',plot = p1,width = 45,height = 20)
####气泡图####
GO=read.table('./nk.sample.go.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
GO=GO[,c(3,8,10,11)]
GO1=read.table('./nk.sample.go.selected.txt',sep = '\t',header = T,stringsAsFactors = F)
GO7=GO[GO$Description%in%GO1$Description,]
GO7$subcluster=paste0('SC',GO7$subcluster)
library(forcats)
GO7$Description=fct_infreq(GO7$Description)
GO7$sample=factor(GO7$sample,levels = c('HC','DLE','SLE'))
p <- ggplot(GO7, aes(sample,Description, size = Count, color=qvalue)) + geom_point() +xlab("sample")
#### 气泡图颜色####
p2=p + scale_color_gradient(low='red',high='blue') +theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.major=element_line(colour=NA))+theme(text = element_text(size = 20))
p2
ggsave(filename = 'burbble.pathway.sample.nk.dermis.20220213.pdf',plot = p2,width =10,height = 10)
p=VlnPlot(merged,features = c('COL1A1','COL3A1','CD3D','CD3E','VWF','CDH5','LYZ','AIF1','CD79A','MS4A1','TPSAB1','TPSB2','NKG7','XCL1','CDH19','MPZ'),stack = T,pt.size = 0)
p
ggsave(filename = 'nk.dermis.celltypemarkers.pdf',plot = p,width = 10,height = 6)
