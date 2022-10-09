####cellphonedb可视化####
setwd('D:/SCLE/cellphonedb')
library(tidyverse)
library(RColorBrewer)
library(scales)
pvalues=read.table("./scle/out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
pvalues=pvalues[,12:dim(pvalues)[2]] #此时不关注前11列
statdf=as.data.frame(colSums(pvalues < 0.05)) #统计在某一种细胞pair的情况之下，显著的受配体pair的数目；阈值可以自己选
colnames(statdf)=c("number")

#排在前面的分子定义为indexa；排在后面的分子定义为indexb
statdf$indexb=str_replace(rownames(statdf),"^.*\\.","")
statdf$indexa=str_replace(rownames(statdf),"\\..*$","")
statdf$indexb=ifelse(statdf$indexb%in%c('Macro','DC'),'Macro/DC',statdf$indexb)
statdf$indexa=ifelse(statdf$indexa%in%c('Macro','DC'),'Macro/DC',statdf$indexa)

#设置合适的细胞类型的顺序
rankname=c('T','Kera','Fib','Macro/DC','vEC','LEC','B','Plasma','Mast','NK','Melano','Schwann') 
#转成因子类型，画图时，图形将按照预先设置的顺序排列
table(statdf$indexb)
statdf$indexa=factor(statdf$indexa,levels = rankname)
statdf$indexb=factor(statdf$indexb,levels = rankname)

statdf%>%ggplot(aes(x=indexa,y=indexb,fill=number))+geom_tile(color="white")+
              scale_fill_gradientn(colours = c("#4393C3","#ffdbba","#B2182B"),limits=c(0,20))+
              scale_x_discrete("cluster 1 produces molecule 1")+
              scale_y_discrete("cluster 2 produces molecule 2")+
              theme_minimal()+
              theme(
                            axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45),
                            panel.grid = element_blank()
              )
ggsave(filename = "interaction.num.hc.pdf",device = "pdf",width = 10,height = 8,units = c("cm"))
write.table(statdf,file = 'interaction.number.scle.txt',sep = '\t',row.names = T)
write.table(statdf,file = 'interaction.number.hc.txt',sep = '\t',row.names = T)
write.table(statdf,file = 'interaction.number.dle.txt',sep = '\t',row.names = T)

he.interaction=read.table('./hc/hc.dermis.cell.communication.number.txt',sep = '\t',header = T,row.names = 1)
pheatmap::pheatmap(he.interaction,scale = 'none',cluster_rows = F,cluster_cols = F)
dle.interaction=read.table('./dle/dle.dermis.cell.communication.number.txt',sep = '\t',header = T,row.names = 1)
sle.interaction=read.table('./sle/sle.dermis.cell.communication.number.txt',sep = '\t',header = T,row.names = 1)
he.interaction                   
row.names(he.interaction)=c(paste0(row.names(he.interaction),'.HC'))
colnames(he.interaction)=c(paste0(colnames(he.interaction),'.HC'))
row.names(dle.interaction)=c(paste0(row.names(dle.interaction),'.DLE'))
colnames(dle.interaction)=c(paste0(colnames(dle.interaction),'.DLE'))
row.names(sle.interaction)=c(paste0(row.names(sle.interaction),'.SLE'))
colnames(sle.interaction)=c(paste0(colnames(sle.interaction),'.SLE'))
write.table(he.interaction,file = 'hc.number.txt',sep = '\t',row.names = T)
write.table(dle.interaction,file = 'dle.number.txt',sep = '\t',row.names = T)
write.table(sle.interaction,file = 'sle.number.txt',sep = '\t',row.names = T)
all.interaction=read.table('./all.dermis.number.txt',sep = '\t',header = T,row.names = 1)
pheatmap::pheatmap(all.interaction,scale = 'none',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',annotation_row = row.annotation,annotation_col = col.annotation,angle_col = 45)
col.annotation=data.frame(sample=c(rep('HC',8),rep('DLE',8),rep('SLE',8)),row.names = colnames(all.interaction))
row.annotation=data.frame(sample=c(rep('HC',8),rep('DLE',8),rep('SLE',8)),row.names = row.names(all.interaction))
rm(list = ls())


all.interaction=read.table('./all.dermis.number.1.txt',sep = '\t',header = T,row.names = 1)
pheatmap::pheatmap(all.interaction,scale = 'none',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation',display_numbers = T,number_color = 'white')
#col.annotation=data.frame(sample=c(rep('HC',8),rep('DLE',8),rep('SLE',8)),row.names = colnames(all.interaction))
#row.annotation=data.frame(sample=c(rep('HC',8),rep('DLE',8),rep('SLE',8)),row.names = row.names(all.interaction))
all.interaction=read.table('./fib.all.number.dermis.txt',sep = '\t',header = T,row.names = 1)
pheatmap::pheatmap(all.interaction,scale = 'none',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation',display_numbers = T,number_color = 'white')

all.interaction=read.table('./fib.EC.all.number.dermis.txt',sep = '\t',header = T,row.names = 1)
pheatmap::pheatmap(all.interaction,scale = 'none',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation',display_numbers = T,number_color = 'black')

####点图####
setwd('../')
new_path <- './scle/out/'
mypvals <- read.delim(file.path(new_path,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(file.path(new_path, "means.txt"), check.names = FALSE)
mysigmeans <- read.delim(file.path(new_path, "significant_means.txt"), check.names = FALSE)
dim(mymeans)
mymeans[1:4, 1:4]
# 调整受配体位置，配体在前，受体在后Rearrange data column sequence
library(dplyr)
order_sequence <- function(df){
              da <- data.frame()
              for(i in 1:length(df$gene_a)){
                            sub_data <- df[i, ]
                            if(sub_data$receptor_b=='False'){
                                          if(sub_data$receptor_a=='True'){
                                                        old_names <- colnames(sub_data)
                                                        my_list <- strsplit(old_names[-c(1:11)], split="\\|")
                                                        my_character <- paste(sapply(my_list, '[[', 2L), sapply(my_list, '[[', 1L), sep='|')
                                                        new_names <- c(names(sub_data)[1:4], 'gene_b', 'gene_a', 'secreted', 'receptor_b', 'receptor_a', "annotation_strategy", "is_integrin", my_character)
                                                        sub_data = dplyr::select(sub_data, new_names)
                                                        # print('Change sequence!!!')
                                                        names(sub_data) <- old_names
                                                        da = rbind(da, sub_data) 
                                          }
                            }else{
                                          da = rbind(da, sub_data)
                            }
              }
              return(da)
}

# 受体配体对，要求一个为受体，一个为配体
df <- subset(mymeans, receptor_a == 'True' & receptor_b == 'False' | receptor_a == 'False' & receptor_b == 'True')

# 筛选单基因的受体配体对，待商榷！
df <- df %>% dplyr::mutate(na_count = rowSums(is.na(df) | df == "")) %>% subset(na_count == 0) %>% dplyr::select(-na_count)
dim(df)
## [1] 72 27

# 运行函数，重排受体配体顺序，并生成新列
means_order <- order_sequence(df) %>% tidyr::unite(Pairs, gene_a, gene_b)
pvals_order <- order_sequence(mypvals) %>% tidyr::unite(Pairs, gene_a, gene_b)

# 保存文件
write.table(means_order %>% dplyr::select(-interacting_pair), file.path(new_path,"means_order.txt"), sep = '\t', quote=F, row.names=F)
write.table(pvals_order %>% dplyr::select(-interacting_pair), file.path(new_path,"pvalues_order.txt"), sep = '\t', quote=F, row.names=F)
#合并表达量文件和pvalue文件#
setwd('./scle/')
means_sub <- means_order[, c('Pairs', colnames(mymeans)[-c(1:11)])]
pvals_sub <- pvals_order[, c('Pairs', colnames(mymeans)[-c(1:11)])]
means_gather <- tidyr::gather(means_sub, celltype, mean_expression, names(means_sub)[-1])
pvals_gather <- tidyr::gather(pvals_sub, celltype, pval, names(pvals_sub)[-1])
mean_pval <- dplyr::left_join(means_gather, pvals_gather, by = c('Pairs', 'celltype'))
library(DT)

create_dt <- function(x){
              DT::datatable(x,
                            extensions = 'Buttons',
                            options = list(dom = 'Blfrtip',
                                           buttons = c('csv', 'excel', 'pdf'),
                                           lengthMenu = list(c(10,25,50,-1), c(10,25,50,"All"))))
}
create_dt(mean_pval)
write.table(mean_pval, "mean_pval_order.txt", sep = '\t', quote=F, row.names=F)
# 提取显著性表达的受体配体对 #
# 至少在一组细胞类型两两组合中，pvalue显著的受体配体对，认为是显著性表达的受体配体对
?unique.data.frame
mean_pval.1 =mean_pval%>% distinct(Pairs, celltype, .keep_all = TRUE)
dim(mean_pval.1)
dim(mean_pval)
a <- mean_pval.1 %>% dplyr::select(Pairs, celltype, pval) %>% tidyr::spread(key=celltype, value=pval)
sig_pairs <- a[which(rowSums(a<=0.05)!=0), ]
dim(sig_pairs)
?select
?spread
# 保存显著性表达的受体配体对
mean_pval_sub <- subset(mean_pval, Pairs %in% sig_pairs$Pairs)
write.table(mean_pval_sub, "mean_pval_sig.txt", sep = '\t', quote=F, row.names=F)
#可视化受配体对#
library(RColorBrewer)
library(scales)
library(ggplot2)
library(cowplot)

# 比较均值和P值数据变换前后的分布
p1 <- mean_pval_sub %>% ggplot(aes(x=mean_expression)) + geom_density(alpha=0.2, color='red')
p2 <- mean_pval_sub %>% ggplot(aes(x=log2(mean_expression))) + geom_density(alpha=0.2, color='red')
mean_distribution <- plot_grid(p1, p2, labels = "AUTO", nrow=1)

p1 <- mean_pval_sub %>% ggplot(aes(x=pval)) + geom_density(alpha=0.2, color='red')
p2 <- mean_pval_sub %>% ggplot(aes(x=-log10(pval+1*10^-3))) + geom_density(alpha=0.2, color='red')
pval_distribution <- plot_grid(p1, p2, labels = "AUTO", nrow=1)
plot_grid(mean_distribution, pval_distribution, labels = "AUTO", ncol=1)

# 展示部分（10对）显著性表达的受体配体对
p <- mean_pval_sub %>% dplyr::arrange(Pairs) %>% head(16*40) %>% 
              ggplot(aes(x=Pairs, y=celltype)) +
              # geom_point(aes(color=log2(mean_expression), size=pval)) +
              # scale_size(trans = 'reverse') +
              geom_point(aes(color=log2(mean_expression), size=-log10(pval+1*10^-3)) ) +
              guides(colour = guide_colourbar(order = 1),size = guide_legend(order = 2)) +
              labs(x='', y='') +
              scale_color_gradientn(name='Expression level \n(log2 mean expression \nmolecule1, molecule2)', colours = terrain.colors(100)) +
              # scale_color_gradient2('Expression level \n(log2 mean expression \nmolecule1, molecule2)', low = 'blue', mid = 'yellow', high = 'red') +
              theme(axis.text.x= element_text(angle=45, hjust=1)) +
              # coord_flip() +
              theme(
                            panel.border = element_rect(color = 'black', fill = NA),
                            panel.grid.major.x = element_blank(),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.background   = element_blank(),
                            axis.title.x       = element_blank(),
                            axis.title.y       = element_blank(),
                            axis.ticks         = element_blank()
                            # plot.title         = element_text(hjust = 0.5),
                            # legend.position = 'bottom' # guides(fill = guide_legend(label.position = "bottom"))
                            # legend.position    = "bottom"
                            # axis.text.y.right  = element_text(angle=270, hjust=0.5)
              ) +
              theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
                    # legend.key.height = unit(1, 'cm'), #change legend key height
                    # legend.key.width = unit(1, 'cm'), #change legend key width
                    legend.title = element_text(size=9), #change legend title font size
                    legend.text = element_text(size=8)) #change legend text font size

p
save_plot(paste0(new_path,"/Ligands_Receptors_dotplot.png"), p, base_height = 10, base_aspect_ratio = 1, base_width = 6, dpi=600)
save_plot(paste0(new_path,"/Ligands_Receptors_dotplot.pdf"), p, base_height = 10, base_aspect_ratio = 1, base_width = 6)
####成纤维点图####
fib=read.table('./fib.mean.pval.sig.txt',sep = '\t',header = T)
#fib1=separate(fib,col = celltype,into = c('cell1','cell2'),sep = '|')
?separate
fib1=fib[fib$celltype=='Fibroblasts'|fib$X=='Fibroblasts',]
fib.hc=fib1[fib1$sample=='HC',]
fib.dle=fib1[fib1$sample=='DLE',]
fib.sle=fib1[fib1$sample=='SLE',]
fib.hc=fib.hc[fib.hc$celltype!='Mast_cells',]
fib.hc=fib.hc[fib.hc$celltype!=fib.hc$X,]

fib.dle=fib.dle[fib.dle$celltype!='Endothelial_cells',]
fib.dle=fib.dle[fib.dle$X!='Endothelial_cells',]
fib.dle=fib.dle[fib.dle$X!=fib.dle$celltype,]
fib.dle=fib.dle[fib.dle$X!='Schwann_cells',]
fib.dle=fib.dle[fib.dle$celltype!='Schwann_cells',]
fib.dle=fib.dle[fib.dle$X!='Mast_cells',]
fib.dle=fib.dle[fib.dle$celltype!='Mast_cells',]

fib.sle=fib.sle[fib.sle$celltype!='Endothelial_cells',]
fib.sle=fib.sle[fib.sle$X!='Endothelial_cells',]
fib.sle=fib.sle[fib.sle$X!=fib.sle$celltype,]
fib.sle=fib.sle[fib.sle$X!='Schwann_cells',]
fib.sle=fib.sle[fib.sle$celltype!='Schwann_cells',]
fib.sle=fib.sle[fib.sle$X!='Mast_cells',]
fib.sle=fib.sle[fib.sle$celltype!='Mast_cells',]
a=fib.hc$Pairs
b=fib.dle$Pairs
c=fib.sle$Pairs
d=intersect(intersect(a,b),c)
fib.hc.1=fib.hc[fib.hc$Pairs%in%d,]
fib.dle.1=fib.dle[fib.dle$Pairs%in%d,]
fib.sle.1=fib.sle[fib.sle$Pairs%in%d,]
fib.dle.1$mean_expression.1
a=fib.hc.1$Pairs
b=fib.dle.1$Pairs
c=fib.sle.1$Pairs
names(fib.hc.1)
names(fib.dle.1)
fib.hc.1$cell_cell=c(paste(fib.hc.1$celltype,fib.hc.1$X,sep = '|'))
fib.dle.1$cell_cell=c(paste(fib.dle.1$celltype,fib.dle.1$X,sep = '|'))
fib.sle.1$cell_cell=c(paste(fib.sle.1$celltype,fib.sle.1$X,sep = '|'))
fib.hc.1$Pairs.1=c(paste(fib.hc.1$Pairs,fib.hc.1$cell_cell,sep = '_'))
fib.dle.1$Pairs.1=c(paste(fib.dle.1$Pairs,fib.dle.1$cell_cell,sep = '_'))
fib.sle.1$Pairs.1=c(paste(fib.sle.1$Pairs,fib.sle.1$cell_cell,sep = '_'))


fib.immune=merge(x = fib.hc.1,y=fib.dle.1, by = 'Pairs.1')

fib.immune=merge(x=fib.immune,y=fib.sle.1,by='Pairs.1')
write.table(fib.immune,file = 'fib.immune.txt',sep = '\t')
 write.table(fib.hc,file = 'fib.hc.mean.pval.txt',sep = '\t')
write.table(fib.dle,file = 'fib.dle.mean.pval.txt',sep = '\t')
write.table(fib.sle,file = 'fib.sle.mean.pval.txt',sep = '\t')
fib.dle$hc.mean=ifelse(fib.dle$Pairs%in%fib.hc$Pairs,fib.dle$mean_expression/fib.hc$mean_expression,0)
fib.sle$hc.mean=ifelse(fib.sle$Pairs%in%fib.hc$Pairs,fib.sle$mean_expression/fib.hc$mean_expression,0)
a=fib.dle$Pairs
b=fib.hc$Pairs
c=intersect(intersect(a,b),d)
fib2=fib1[fib1$Pairs%in%c,]
fib2.hc=fib2[fib2$sample=='HC',]
fib2.dle=fib2[fib2$sample=='DLE',]
fib2.sle=fib2[fib2$sample=='SLE',]
d=fib.sle$Pairs
e=fib.dle$Pairs[fib.dle$hc.mean>=10]
e
f=fib.sle$Pairs[fib.sle$hc.mean>=10]
f
g=intersect(e,f)
g
g1=fib1[fib1$Pairs%in%g,]
g1$cell_cell=c(paste(g1$celltype,g1$X,sep = '|'))
g1$sample=factor(g1$sample,levels = c('HC','DLE','SLE'))
g1$cell_cell.1=c(paste(g1$cell_cell,g1$sample,sep = '.'))


table(fib.dle$hc.mean==0)
a=fib.dle$Pairs[fib.dle$hc.mean==0]
b=fib.sle$Pairs[fib.sle$hc.mean==0]
d=a==b
table(d)
d=a[d]
d
fib.gene=fib1[fib1$Pairs%in%d,]

fib.gene$cell_cell=c(paste(fib.gene$celltype,fib.gene$X,sep = '|'))

fib.gene$cell_cell.1=c(paste(fib.gene$cell_cell,fib.gene$sample,sep = '.'))
p <- ggplot(g1, aes(cell_cell.1, Pairs,size = log2(mean_expression), color=pval)) + geom_point() +xlab("cell_cell")
#### 气泡图颜色####
p2=p + scale_color_gradient(low='red',high='blue') +theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.major=element_line(colour=NA))
p3=p2+theme(axis.text.x = element_text(angle = 90,hjust = 0.5, vjust = 0.5))

p3
ggsave(filename = 'burbble.pathway.subT20201127.pdf',plot = p2,width = 15,height = 20)
table(abs(fib.dle$hc.mean)>4)
fib.dle$Pairs[abs(fib.dle$hc.mean)>0.5]
a
b
d=intersect(a,b)
d
#dle及sle中有，HC中没有的
fib.gene=fib1[fib1$Pairs%in%d,]

write.table(fib.dle,file = 'fib.dle.txt',sep = '\t')
write.table(fib.sle,file = 'fib.sle.txt',sep = '\t')
hc=read.table('./hc/hc.receptor.ligand.txt',header = T,row.names = 1)
library(pheatmap)
p=pheatmap(hc,scale = 'row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation')
ggsave(filename = './hc/hc.dermis.receptor.ligand.communication.score.pdf',plot = p,width = 8,height = 6)
dle=read.table('./dle/dle.receptor.ligand.txt',header = T,row.names = 1)
p=pheatmap(dle,scale = 'row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation')
ggsave(filename = './dle/dle.dermis.receptor.ligand.communication.score.pdf',plot = p,width = 8,height = 6)
sle=read.table('./sle/sle.receptor.ligand.txt',header = T,row.names = 1)
p=pheatmap(sle,scale = 'row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation')
ggsave(filename = './sle/sle.dermis.receptor.ligand.communication.score.pdf',plot = p,width = 8,height = 6)
colnames(hc)=c(paste(colnames(hc),'HC',sep = '.'))
colnames(dle)=c(paste(colnames(hc),'DLE',sep = '.'))
colnames(sle)=c(paste(colnames(hc),'SLE',sep = '.'))
all=cbind(hc,dle,sle)
p=pheatmap(all,scale = 'row',cluster_rows = F,cluster_cols = F,annotation_col = anno.col,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation',display_numbers = T,number_color = 'black')
colnames(all)=c(rep(c('B','Endothelial_cells','Fibroblasts','Macro/DC','Mast_cells','NK','Schwann_cells','T'),3))
anno.col=data.frame(
                    sample=c(rep('HC',8),rep('DLE',8),rep('SLE',8)),
                    row.names = colnames(all))
anno.row=data.frame(type=c(rep('ligand',8)),row.names = row.names(all))
ggsave(filename = './all.dermis.receptor.ligand.communication.score.number.pdf',plot = p,width = 14,height = 6)
p=pheatmap(all,scale = 'row',cluster_rows = F,cluster_cols = F,annotation_col = anno.col,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation')
p
ggsave(filename = './all.dermis.receptor.ligand.communication.score.nonumber.pdf',plot = p,width = 14,height = 6)
####表皮####
setwd('E:/R_analysis/epidermis.2022.1/cellphonedb/')
hc=read.table('./hc/hc.epidermis.ligand.receptor.txt',header = T,row.names = 1)
library(pheatmap)
p=pheatmap(hc,scale = 'row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation')
ggsave(filename = './hc/hc.epidermis.receptor.ligand.communication.score.pdf',plot = p,width = 8,height = 6)
dle=read.table('./dle/dle.epidermis.ligand.receptor.txt',header = T,row.names = 1)
p=pheatmap(dle,scale = 'row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation')
ggsave(filename = './dle/dle.epidermis.receptor.ligand.communication.score.pdf',plot = p,width = 8,height = 6)
sle=read.table('./sle/sle.epidermis.ligand.receptor.1.txt',header = T,row.names = 1)
p=pheatmap(sle,scale = 'row',cluster_rows = F,cluster_cols = F,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation')
ggsave(filename = './sle/sle.epidermis.receptor.ligand.communication.score.pdf',plot = p,width = 8,height = 6)
colnames(hc)=c(paste(colnames(hc),'HC',sep = '.'))
colnames(dle)=c(paste(colnames(hc),'DLE',sep = '.'))
colnames(sle)=c(paste(colnames(hc),'SLE',sep = '.'))
all=cbind(hc,dle,sle)
p=pheatmap(all,scale = 'row',cluster_rows = F,cluster_cols = F,annotation_col = anno.col,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation',display_numbers = T,number_color = 'black')
colnames(all)=c(rep(c('B','Endothelial_cells','Fibroblasts','Macro/DC','Mast_cells','NK','Schwann_cells','T'),3))
anno.col=data.frame(
              sample=c(rep('HC',6),rep('DLE',6),rep('SLE',6)),
              row.names = colnames(all))
anno.row=data.frame(type=c(rep('ligand',8)),row.names = row.names(all))
ggsave(filename = './all.epidermis.receptor.ligand.communication.score.number.pdf',plot = p,width = 12,height = 6)
p=pheatmap(all,scale = 'row',cluster_rows = F,cluster_cols = F,annotation_col = anno.col,color = colorRampPalette(c('navy','white', 'red'))(100),border_color = 'white',angle_col = 45,clustering_distance_cols = 'correlation')
p
ggsave(filename = './all.epidermis.receptor.ligand.communication.score.nonumber.pdf',plot = p,width = 12,height = 6)
####可视化点图####
setwd('E:/R_analysis/dermis.2022.1/cellphonedb/')
setwd('../dle')
scle=read.table('./mean_pval_sig.txt',sep = '\t',header = T)
scle.sig=scle[scle$pval<=0.05,]

hc=read.table('../hc/mean_pval_sig.txt',sep = '\t',header = T)
hc.sig=hc[hc$pval<=0.05,]

?str_split
?strsplit
dim(hc.sig)
#hc.sig=hc.sig[,-c(5,6)]
####字符串分割####
a=strsplit(hc.sig$celltype,'|',fixed = T)
cell1 = sapply(a,function(x){x[1]})
cell2=sapply(a,function(x){x[2]})
hc.sig$cell1=cell1
hc.sig$cell2=cell2
hc.sig$sample='HC'
scle.sig$sample='SCLE'

table(scle.sig$cell1)
sle.sig.kera.macro=sle.sig[sle.sig$cell1%in%c('Keratinocytes','Macro/DC')|sle.sig$cell2%in%c('Keratinocytes','Macro/DC'),]
sle.sig.kera.macro=sle.sig.kera.macro[sle.sig.kera.macro$cell1!=sle.sig.kera.macro$cell2,]
write.table(sle.sig.kera.macro,file = 'sle.sig.kera.macro.txt',sep = '\t',row.names = T)
table(hc.sig.EC.Fib.macro$celltype)
hc.sig.kera.macro$sample='HC'
dle.sig.kera.macro$sample='DLE'
sle.sig.kera.macro$sample='SLE'
all=merge(hc.sig,scle.sig,by = c('Pairs','celltype'))
all=rbind(hc.sig.kera.macro,dle.sig.kera.macro,sle.sig.kera.macro)
top10 <- all %>% group_by(sample) %>% top_n(n = 10, wt = mean_expression)
hc.sig.kera.macro$PC=c(paste(hc.sig.kera.macro$Pairs,hc.sig.kera.macro$celltype,sep = '.'))
dle.sig.kera.macro$PC=c(paste(dle.sig.kera.macro$Pairs,dle.sig.kera.macro$celltype,sep = '.'))
sle.sig.kera.macro$PC=c(paste(sle.sig.kera.macro$Pairs,sle.sig.kera.macro$celltype,sep = '.'))
hc.dle=merge(hc.sig.kera.macro,dle.sig.kera.macro,by = 'PC')
hc.dle.sle=merge(hc.dle,sle.sig.kera.macro,by='PC')
all$ratio=all$mean_expression.y/all$mean_expression.x
write.table(all,file = 'hc.scle.interaction.pairs.20220609.txt',sep = '\t',row.names = T)

hc.dle.sle=read.table('../hc.dle.sle.epidermis.KERA.MACRO.txt',header = T)
hc.dle.sle$Pairs[1:10]
a=c('CCL19_CCR7','FAM3C_CLEC2D','TNFSF9_HLA-DPA1','MIF_CD74','TNFSF13B_CD40')
a=c('CCL20_CXCR3','TNFSF4_TNFRSF4','SPN_ICAM1','TNFSF9_TNFRSF9','FAM3C_HLA-C','FAM3C_CLEC2D')
hc.sig$sample='HC'
dle.sig$sample='DLE'
sle.sig$sample='SLE'
hc$sample='HC'
dle$sample='DLE'
sle$sample='SLE'
all=rbind(hc.sig,scle.sig)
all=rbind(hc.sig,dle.sig,sle.sig)
selected.1=all[all$Pairs%in%a,]
#### 气泡图颜色####
selected.1$sample=factor(selected.1$sample,levels = c('HC','SCLE'))
selected.1$cell_cell.1=c(paste(selected.1$celltype,selected.1$sample,sep = '.'))
selected.1$Pairs=fct_infreq(selected.1$Pairs)
table(selected.1$Pairs)
write.table(selected.1,file = '../selected.ggplot.1.txt',sep = '\t')
selected=read.table('../scle/cell.communication.pairs.cellphonedb.txt',sep = '\t',header = T)
selected.1$cell_cell.1=forcats::fct_inorder(selected$cell_cell.1)
selected.1=selected.1[selected.1$cell1%in%c('B','Fib','Kera','Macro/DC','Mast','NK','Plasma','T'),]
table(selected.1$cell1)
selected.1=selected.1[selected.1$cell1!=selected.1$cell2,]
selected.1=selected.1[selected.1$cell2%in%c('B','Fib','Kera','Macro/DC','Mast','NK','Plasma','T'),]

p <- ggplot(selected, aes(Pairs, cell_cell.1 ,size = mean_expression, color=pval)) + geom_point() +xlab("")

p2=p + scale_color_gradient(low='red',high='blue') +theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.major=element_line(colour=NA))
p3=p2+theme(axis.text.x = element_text(angle = 90,hjust = 0.5, vjust = 0.5))

p3
ggsave(filename = '../burbble.scle.cell_cell.communication.20220609.pdf',plot = p2,width = 12,height = 15)
write.table(selected.1,file = 'cell.communication.pairs.cellphonedb.txt',sep = '\t',row.names = T)
getwd()
