####B.epidermis.2022####
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
setwd('E:/R_analysis/epidermis.2022.1')
dir.create('./B')
setwd('./B/')
sam.name <- "multi"
if(!dir.exists(sam.name)){
              dir.create(sam.name)
}

merged=B
rm(B)
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
dim.use=1:15
library(harmony)

merged=RunHarmony(merged,group.by.vars ='orig.ident')

merged<- FindNeighbors(merged, dims = dim.use,reduction = 'harmony',nn.method = "rann")
merged <- FindClusters(merged, resolution =0.1)
merged <- RunTSNE(merged, dims = dim.use,reduction = 'harmony')
merged1=RunUMAP(merged,dims = dim.use,reduction = 'harmony')
pdf(paste0("./",sam.name,"/tsne_Plot.harmony_res0.2_",max(dim.use),"PC.pdf"),width = 5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = T)
dev.off()
FeaturePlot(merged,c('CD79A','MS4A1','CD3D','KRT5'))
save(merged,file = 'b.epidermis.res0.1.20220117.Rdata')
write.table(merged@meta.data,file = paste0("./",sam.name,"/",sam.name,"_cells_details_tsne_res0.2_",max(dim.use),"PC.txt"))
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))

pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.1.nolabel_",max(dim.use),"PC.pdf"),width =5,height = 4)
DimPlot(object = merged, pt.size=0.1,label = F,cols = cols)
dev.off()
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_orig_res0.1_",max(dim.use),"PC.pdf"),width = 5,height = 4)
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
marker.sig <- all.markers.noRibo.noMito[all.markers.noRibo.noMito$p_val<=0.05,]

top10 <- marker.sig %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file=paste0("./",sam.name,"/",sam.name,"_top10_marker.res0.2.nomt.noribo_",max(dim.use),"PC.txt"),sep="\t",quote = F)
a=table(merged$orig.ident,merged$seurat_clusters)
a=as.data.frame(a)
write.table(a,file = 'orig.res0.1.b.subcluster.txt',sep='/t')

p1=ggplot(a,mapping = aes(Var2,Freq,fill=Var1))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='b.subcluster.orig.pdf',plot=p3,width=10,height=4)

b=table(merged$seurat_clusters,merged$sample)
b=as.data.frame(b)
write.table(b,file = 'sample.b.res0.1.sbucluster.txt',sep = '\t')
p1=ggplot(b,mapping = aes(Var1,Freq,fill=Var2))+geom_bar(stat='identity',position='fill') +
              labs(x = 'sample',y = 'frenquency') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'),
                                                         axis.text.x = element_text(angle = 45, hjust = 1))
p2=p1+scale_fill_manual(values=cols)

p3=p2 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3
ggsave(filename='b.res0.1.subcluster.sample.pdf',plot=p3,width=5,height=4)
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
              gene=all.markers.noRibo.noMito$cluster==i0001-4842|2380-8195|1614-6832|0935-9648|1943-8206|0001-8732|1552-5260|1073-449X|0003-4819|0923-7534|0066-4146|0066-4154|0732-0582|1553-4006|1543-5008|0066-4308|0163-7525|0006-4971|0959-8138|1535-6108|2159-8274|0161-4940|0092-8674|1931-3128|1550-4131|1001-0602|1934-5909|2451-9294|0009-2665|0306-0012|0009-7322|0893-8512|0010-8545|1754-5692|0195-668X|0302-2838|1560-2745|0016-5085|0017-5749|1074-7613|1750-984X|2168-6106|2374-2437|2168-622X|0098-7484|2542-4351|0732-183X|2001-3078|0168-8278|0735-1097|0140-6736|2213-8587|2589-7500|2214-109X|1473-3099|1474-4422|1470-2045|2215-0366|2468-2667|2213-2600|1433-8351|0927-796X|1369-7021|1057-5987|1546-0738|1476-4598|1748-0132|0028-0836|2157-846X|1087-0156|2520-1158|1465-7392|1755-4330|1758-678X|2520-1131|2058-7546|1061-4036|1529-2908|1476-1122|1078-8956|1548-7091|1748-3387|1097-6256|1749-4885|1745-2473|1474-175X|1759-5002|2397-3358|1759-4774|2056-676X|1474-1776|1759-5029|1759-5045|1471-0056|1474-1733|2058-8437|1740-1526|1471-0072|1759-5061|1759-4758|1471-003X|2522-5820|1759-4790|0028-4793|0031-6997|0370-1573|0031-9333|0360-1285|0079-6425|0079-6700|1350-9462|1529-1006|8755-1209|0034-6861|0036-8075|2470-9476|0962-8924|2589-7209|1364-6613|1759-0876|1723-8617
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

              write.table(enrich.go.BP,file =paste0(i,'_B.go.txt') ,sep='\t')
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
              
              write.table(enrich.go.BP,file =paste0(i,'_b.go.txt') ,sep='\t')
              write.table(enrich.KEGG,file =paste0(i,'_b.kegg.txt') ,sep='\t')
}
####找sample的marker基因####
table(Idents(merged))
Idents(merged)=merged@meta.data$sample
all.markers.sample <- FindAllMarkers(merged, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers.sample,file=paste0("./",sam.name,"/",sam.name,"_sample_marker_genes_b_",30,"PC.txt"),sep="\t",quote = F)
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
markers=read.table('./top10.subb.average.expression.txt',sep = '\t',header = T,row.names = 1)

markers=t(markers)
library(pheatmap)
p1=pheatmap(markers,cluster_cols = F,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'column',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')
p1
ggsave(filename = 'b.subtype.top10.marker.dermis.expression.pdf',plot = p1,width = 45,height = 20)
####经典marker基因匹配####
Mempry_b=c('CD27')
memory_b=marker[row.names(marker)%in%Mempry_b,]
ABC=c('TBX21','ITGAX', 'FCRL2', 'IL10RA')
abc=marker[row.names(marker)%in%ABC,]
Plasma=c('CD38', 'JCHAIN','MZB1','IGHG4', 'IGHA1','IGHG1','IGHG3')
Plasma=marker[row.names(marker)%in%Plasma,]
isg=c('ISG15','IFIT1','IFITM1','ISG20')
isg=marker[row.names(marker)%in%isg,]
hsp=c('HSPA6','HSPA1A','HSPA1B','HSPB1','HSPA5','HSPH1','HSP90AA1')
hsp=marker[row.names(marker)%in%hsp,]
classical=rbind(memory_b,abc,Plasma,isg,hsp)
anno.row=data.frame(subtype=c(rep('Memory',1),rep('ABC',4),
                              rep('Plasma',7),rep('ISG',4),rep('HSP',7)),
                    row.names = row.names(classical))
p1=pheatmap(classical,annotation_row =anno.row ,clustering_distance_cols  = "correlation",scale = 'row',colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 15,cellheight =12,border_color = 'white',fontsize = 10,cluster_rows = F)
p1
colnames(classical)=c(paste0('SC',0:4))
#a=c(Naive_b,Mempry_b,DN2,Activate_B,ABC,Plasma)
annotation_col=data.frame(subtype = c(rep('Naive_B',4),rep('Memory_B',2),rep('DN2',6),rep('ABC',8),rep('Plasma',7)),row.names = row.names(classical))
length(row.names(classical))
p1=pheatmap(classical,annotation_row =annotation_col ,clustering_distance_cols  = "correlation",scale = 'row',colorRampPalette(colors = c("blue","white","red"))(100),treeheight_row = 20,treeheight_col = 20,cellwidth = 15,cellheight =12,border_color = 'white',fontsize = 10,cluster_rows = F)
p1
ggsave(filename = 'b.classical.subtype.marker.epidermis.pdf',plot = p1,width = 6,height =10)
####气泡图####
GO=read.table('./b.epidermis.all.go.txt',sep = '\t',header = T,stringsAsFactors = F,quote = '')
GO=GO[,c(3,8,10,11)]
GO1=read.table('./b.epidermis.selected.go.txt',sep = '\t',header = T,stringsAsFactors = F)
GO7=GO[GO$Description%in%GO1$Description,]
#GO7$subcluster=paste0('SC',GO7$subcluster)
library(forcats)
library(ggplot2)
GO7$Description=fct_infreq(GO7$Description)
GO7$sample=factor(GO7$sample,levels = c('HC','DLE','SLE'))
p <- ggplot(GO7, aes(sample,Description, size = Count, color=qvalue)) + geom_point() +xlab("sample")
#### 气泡图颜色####
p2=p + scale_color_gradient(low='red',high='blue') +theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.major=element_line(colour=NA))+theme(text = element_text(size = 20))
p2
ggsave(filename = 'burbble.pathway.sample.b.20220214.pdf',plot = p2,width =10,height = 10)
####与经典marker匹配####
rm(list = ls())
classical.marker=read.table('./classical.gene.b.txt',header = T,sep = '\t',row.names = 1)
colnames(classical.marker)
classical.marker=data.frame(subtype=classical.marker$subtype,row.names = row.names(classical.marker))
marker=read.table('./average.expression.subb.removed.cluster.txt')
colnames(marker)=c(paste0('SC',0:4))
abc=row.names(classical.marker)[classical.marker$subtype=='ABC']
abc
abc=marker[row.names(marker)%in%row.names(abc),]
#ABC=c('CD86','SLAMF7', 'IL2RA', 'FAS','TBK1', 'EBI3', 'ZBTB32','TRAF5','ITGAX','TBX21','FCRL2')
#abc=marker[row.names(marker)%in%ABC,]

DN2=row.names(classical.marker)[classical.marker$subtype=='DN2']
DN2
AD=c(abc,DN2)
AD=c()
dn2=marker[row.names(marker)%in%row.names(DN2),]
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
gene=c(ABC,DN2,HSP,ISG,Memory,Naive,plasma)
genes=rbind(abc,dn2,hsp,isg,mem,Naive,plasma.1)
VlnPlot(merged,features = gene,stack = T,flip = T)+NoLegend()
classical.marker=data.frame(subtype=classical.marker$subtype,row.names = row.names(classical.marker))
write.table(gene,file = 'classical.gene.subcluster.b.epidermis.txt',sep = '\t')
row.names(classical.marker)
classsical.marker.1=classical.marker[!classical.marker$subtype=='ABC',]
genes=row.names(classical.marker)[classical.marker$subtype!='ABC']
classical.marker=data.frame(subtype=classsical.marker.1,row.names = genes)
abc=data.frame(subtype=rep('ABC',length(ABC)),row.names = row.names(abc))
classical.marker=rbind(classical.marker,abc)
classical.marker=classical.marker[classical.marker$subtype!='DN2',]
gene=read.table('classical.gene.subcluster.b.epidermis.txt',sep = '\t',header = T,row.names = 1)
library(pheatmap)
row.names(genes)
genes=genes[-13,]
genes=genes[-c(19),]
cd69=marker[row.names(marker)%in%c('CD69'),]
genes=rbind(genes,cd69)
row.names(genes)
genes=genes[-28,]
classical.marker=ifelse(classical.marker$subtype%in%c('DN2','ABC'),'ABC/DN2',classical.marker$subtype)
classical.marker=data.frame(subtype=classical.marker,row.names = )
p1=pheatmap(genes,cluster_cols = F,annotation_row =classical.marker ,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')+NoLegend()
p1
ggsave(filename = 'B.classical.gene.epidermis.expression.2022.pdf',plot = p1,width = 45,height = 20)
gene=t(gene)
p1=pheatmap(genes,cluster_cols = T,
            cluster_rows = F,
            clustering_distance_rows = "correlation",
            scale = 'row',
            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
            treeheight_row = 20,treeheight_col = 20,
            cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',
            legend = T,border_color = 'white',angle_col = 45)
p1
ggsave(filename = 'B.classical.gene.dermis.expression.pdf',plot = p1,width = 45,height = 20)
p=VlnPlot(object = merged,features = c('KRT14','KRT1','CD3D','CD3G','LYZ','AIF1','PMEL','MLANA','NKG7','XCL2','CD79A','MS4A1'),stack = T,pt.size = 0)
ggsave(filename = 'b.epidermis.celltype.genes.pdf',plot = p,width = 8,height = 5)
####
rm(list = ls())
classical.marker=read.table('./classical.gene.b.txt',header = T,sep = '\t',row.names = 1)
colnames(classical.marker)
classical.marker=data.frame(subtype=classical.marker$subtype,row.names = row.names(classical.marker))
marker=read.table('./average.expression.subb.removed.cluster.txt')
colnames(marker)=c(paste0('SC',0:4))
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
colnames(genes.2)=c(paste0('B_',colnames(genes.2)))
p1=pheatmap(genes.2,cluster_cols = F,annotation_row =anno.row,cluster_rows = F,clustering_distance_rows = "correlation",scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_row = 20,treeheight_col = 20,cellwidth = 40,cellheight =25,fontsize = 22,legend_labels = 'scale.avg_logFC',legend = T,border_color = 'white')+NoLegend()
p1
genes.2
row.names(genes.2)
genes.1=genes[c(1:5),]
genes.2=genes[-c(1:5),]
genes.2=genes.2[-c(12,5),]
