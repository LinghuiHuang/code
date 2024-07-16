top10=read.table(file = 'analysis/4_PPI/top10.txt',sep = "\t",header = T)
top10=top10$name

ppi.hub.genes=c(intersect(top10,mcode.cluster1),
                intersect(top10,mcode.cluster2))
length(ppi.hub.genes) ## 7

library(pROC)
setdiff(ppi.hub.genes,rownames(gse57065.exp))
rownames(gse57065.exp)[grep("H2A",rownames(gse57065.exp))]
rownames(gse57065.exp)[which(rownames(gse57065.exp)=='HIST2H2AA3')]='H2AC18'

#########################
gse57065_data=cbind(gse57065.cli,t(gse57065.exp[,rownames(gse57065.cli)]))
gse65682_data=cbind(gse65682.cli,t(gse65682.exp[,rownames(gse65682.cli)]))
gse95233_data=cbind(gse95233.cli,t(gse95233.exp[,rownames(gse95233.cli)]))
gse185263_data=cbind(gse185263.cli,t(gse185263.exp[,rownames(gse185263.cli)]))

get_ROC=function(dt=NULL){
  colnames(dt) <- c('disease', 'Gene')
  roc.res <- roc(disease ~ Gene,
                 data=dt,
                 aur=TRUE,
                 ci=TRUE, 
                 smooth=F)
  auc=round(roc.res$auc,3)
}

gse57065.hub.genes.auc=c()
for(gene in ppi.hub.genes){
  dt <- gse57065_data[, c('disease', gene)]
  auc=get_ROC(dt)
  gse57065.hub.genes.auc=c(gse57065.hub.genes.auc,auc)
}

gse65682.hub.genes.auc=c()
for(gene in ppi.hub.genes){
  dt <- gse65682_data[, c('disease', gene)]
  auc=get_ROC(dt)
  gse65682.hub.genes.auc=c(gse65682.hub.genes.auc,auc)
}

gse95233.hub.genes.auc=c()
for(gene in ppi.hub.genes){
  dt <- gse95233_data[, c('disease', gene)]
  auc=get_ROC(dt)
  gse95233.hub.genes.auc=c(gse95233.hub.genes.auc,auc)
}

gse185263.hub.genes.auc=c()
for(gene in ppi.hub.genes){
  dt <- gse185263_data[, c('disease', gene)]
  auc=get_ROC(dt)
  gse185263.hub.genes.auc=c(gse185263.hub.genes.auc,auc)
}

bulk.hub.genes.auc=rbind(gse57065.hub.genes.auc,gse65682.hub.genes.auc,
                         gse95233.hub.genes.auc,gse185263.hub.genes.auc)
rownames(bulk.hub.genes.auc)=c('GSE57065','GSE65682','GSE95233','GSE185263')
colnames(bulk.hub.genes.auc)=c(ppi.hub.genes)
bulk.hub.genes.auc=data.frame(bulk.hub.genes.auc)
bulk.hub.genes.auc
#           LTF CEACAM8 PGLYRP1  MMP8 MAPK14 S100A12  LCN2
# GSE57065  0.896   0.885   0.993 0.983  1.000   0.996 0.933
# GSE65682  0.906   0.871   0.977 0.971  0.982   0.997 0.930
# GSE95233  0.885   0.892   0.986 0.991  0.981   0.999 0.930
# GSE185263 0.691   0.621   0.826 0.841  0.786   0.956 0.793
########################
p.all <- list()
for (gene in ppi.hub.genes) {
  dt <- gse57065_data[, c('disease', gene)]
  colnames(dt) <- c('disease', 'Gene')
  roc.res <- roc(disease ~ Gene,
                 data=dt,
                 aur=TRUE,
                 ci=TRUE, 
                 smooth=F)
  # roc 
  p <- ggroc(roc.res, legacy.axes = TRUE )+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color="darkgrey", linetype=8)+
    theme_gray() +
    annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(roc.res$auc,3))) +
    ggtitle(paste0(gene, ' ROC'))
  p.all[[gene]] <- p
}
length(p.all)

ppi.hub.genes.roc <- cowplot::plot_grid(plotlist = p.all,
                                             ncol = 4)
ppi.hub.genes.roc

# ggsave(plot = ppi.hub.genes.roc,
#        filename = 'analysis/5_Model/ppi.hub.genes.roc.pdf',
#        width = 16, height = 8)

library(pROC)
dt <- gse57065_data[, c('disease', ppi.hub.genes)]
table(dt$disease)
# healthy  sepsis 
# 25      82
pdf("analysis/5_Hub_Genes/ppi.hub.genes.roc.pdf",width=5,height=5)
plots=list()
AUCs=c()
for (gene in ppi.hub.genes) {
  pt=mg_surv_pROC_smooth_trycatch(status = dt$disease,score = dt[,gene],col = mycolors[length(plots)+1],add = length(plots)>0)
  AUCs=c(AUCs,round(pt$auc,2))
  plots=c(plots,list(pt))
}
legend('bottomright',legend = paste0(ppi.hub.genes,', AUC=',AUCs),col = mycolors[1:length(AUCs)],
       # lty = c(1,1),border = T,title = 'AUC',bty="n",
       lwd=2)
dev.off()

################
clin.color=mycolors[1:2]
dt <- gse57065_data[, c('disease', ppi.hub.genes)]
p.all=list()
for(gene in ppi.hub.genes){
  p=mg_violin_use(data = dt[,c('disease',gene)]
                ,melt = T
                ,group.col = clin.color
                ,xlab = ''
                ,ylab = "Expression"
                ,jitter=F
                ,test_method = 'kruskal.test'
                ,cmp_test_method = 'wilcox.test'
                ,legend.pos = NULL
                ,show_compare = T)
  p.all=c(p.all,list(p))
}
length(p.all)
figure=mg_merge_plot(p.all,nrow = 2,ncol = 4)
figure
savePDF('analysis/5_Hub_Genes/ppi.hub.genes.tissue.boxplot.pdf',figure,height = 10,width = 20)

################# PPI 基因与 PCD的关系
cor.res=psych::corr.test(x=data.frame(t(gse57065.exp[ppi.hub.genes,])),
                         y = data.frame(gse57065.pcd.ssgsea,check.names = F))
df_cor=cor.res$r
df_pval=cor.res$p.adj
########
library(tidyverse)
g = pivot_longer(data=rownames_to_column(data.frame(df_cor,check.names = F),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(data.frame(df_pval,check.names = F)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
all(g$from==gp$rowname & g$to==gp$gene)
g$p.adj = gp$p
head(g)
################## 
df=data.frame(g,check.names = F,stringsAsFactors = F)
cor.R <- pivot_wider(data=g[,c(1:2,3)],names_from ='to',values_from ='cor')
cor.R=data.frame(cor.R)
rownames(cor.R)=cor.R$from
cor.R=cor.R[,-1]
cor.R=as.matrix(cor.R)

cor.Pval <- pivot_wider(data=g[,c(1:2,4)],names_from ='to',values_from ='p.adj')
cor.Pval=data.frame(cor.Pval)
rownames(cor.Pval)=cor.Pval$from
cor.Pval=cor.Pval[,-1]
cor.Pval=as.matrix(cor.Pval)
cor.Pval <- ifelse(cor.Pval < 0.05, '*', '')
cor.Pval[is.na(cor.Pval)] <- ''

library(ComplexHeatmap)
p=Heatmap(cor.R,
          name = "corr",
          col = circlize::colorRamp2(c(-1, 0, 1), c('#3B4992FF', 'white', '#EE0000FF')),
          border = F,
          show_column_names = T,
          show_column_dend = F,
          show_row_dend = F,
          cluster_columns = T,
          cluster_rows = T,
          column_names_side = 'bottom',
          row_names_side = 'left',
          cell_fun = function(j, i, x, y, w, h, col) {
            # add text to each grid
            grid.text(cor.Pval[i, j], x, y)
          })

pcd.sorted=colnames(cor.R)[column_order(p)]
pcd.sorted

####################
df=g
df <- df %>%
  mutate(type = cut(cor, breaks = c(-1, 0, 1),
                    labels = c("negative", "positive")),
         p.signif = cut(p.adj, breaks = c(0,0.0001, 0.001, 0.01, 0.05,1),
                        labels = c("****", "**", "**","*","ns"),
                        right = FALSE, include.lowest = TRUE))
head(df)
df$to=factor(df$to,levels = rev(pcd.sorted))
head(df)
write.table(df,file = 'analysis/5_Hub_Genes/ppi.hub.genes.PCD.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)

hub.genes.pcd.corplot=ggplot(df,aes(y = to,x = from))+
  geom_point(aes(size = abs(cor), color = (cor))) +
  geom_text(aes(y = to, x = from, label = p.signif)) +
  theme_bw() + 
  # scale_color_continuous(values = c('#0072B5','#EE4C97')) +
  # scale_colour_gradient2(low="#0072B5", mid = "white",high="#EE4C97",midpoint = 0)+
  scale_colour_gradient2(high="#F97B72", mid = "white",low="#0072B5",midpoint = 0)+
  ylab('') + xlab('') + labs(color = "pearson's cor",size="pearson's cor")
hub.genes.pcd.corplot

############### 
get.IOBR.immu.format=function(tcga.t.exp.cibersort){
  tcga.t.exp.cibersort = data.frame(tcga.t.exp.cibersort)
  rownames(tcga.t.exp.cibersort) = tcga.t.exp.cibersort$ID
  tcga.t.exp.cibersort = tcga.t.exp.cibersort[, -1]
  colnames(tcga.t.exp.cibersort) = gsub('(.*)_.*', "\\1", colnames(tcga.t.exp.cibersort))
  return(tcga.t.exp.cibersort)
}

load("raw_datas/Pathways_ssGSEA/gse57065.exp.estimate.RData")
load("raw_datas/Pathways_ssGSEA/gse57065.exp.cibersort.RData")

gse57065.exp.estimate=get.IOBR.immu.format(gse57065.exp.estimate)
gse57065.exp.cibersort=get.IOBR.immu.format(gse57065.exp.cibersort)
gse57065.exp.cibersort=gse57065.exp.cibersort[,1:22]
dim(gse57065.exp.cibersort)

gse57065.tme=cbind(gse57065.exp.estimate,gse57065.exp.cibersort)
dim(gse57065.tme)

tme.type=rbind(data.frame(Cell=colnames(gse57065.exp.estimate),type="ESTIMATE"),
               data.frame(Cell=colnames(gse57065.exp.cibersort),type="CIBERSORT"))
head(tme.type)
#############
cor.res=psych::corr.test(x=data.frame(t(gse57065.exp[ppi.hub.genes,])),
                         y = data.frame(gse57065.tme,check.names = F))
df_cor=cor.res$r
df_pval=cor.res$p.adj
########
library(tidyverse)
g = pivot_longer(data=rownames_to_column(data.frame(df_cor,check.names = F),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(data.frame(df_pval,check.names = F)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
all(g$from==gp$rowname & g$to==gp$gene)
g$p.adj = gp$p
g$TME=tme.type$type[match(g$to,tme.type$Cell)]
head(g)
################## 
df=data.frame(g,check.names = F,stringsAsFactors = F)
cor.R <- pivot_wider(data=g[,c(1:2,3)],names_from ='to',values_from ='cor')
cor.R=data.frame(cor.R)
rownames(cor.R)=cor.R$from
cor.R=cor.R[,-1]
cor.R=as.matrix(cor.R)

cor.Pval <- pivot_wider(data=g[,c(1:2,4)],names_from ='to',values_from ='p.adj')
cor.Pval=data.frame(cor.Pval)
rownames(cor.Pval)=cor.Pval$from
cor.Pval=cor.Pval[,-1]
cor.Pval=as.matrix(cor.Pval)
cor.Pval <- ifelse(cor.Pval < 0.05, '*', '')
cor.Pval[is.na(cor.Pval)] <- ''
####################
df=g
df <- df %>%
  mutate(type = cut(cor, breaks = c(-1, 0, 1),
                    labels = c("negative", "positive")),
         p.signif = cut(p.adj, breaks = c(0,0.0001, 0.001, 0.01, 0.05,1),
                        labels = c("****", "**", "**","*","ns"),
                        right = FALSE, include.lowest = TRUE))
head(df)
write.table(df,file = 'analysis/5_Hub_Genes/ppi.hub.genes.TME.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)

p.all=list()
for(i in c('ESTIMATE','CIBERSORT')){
  dt=df[which(df$TME==i),]
  p=ggplot(dt,aes(y = to,x = from))+
    geom_point(aes(size = abs(cor), color = (cor))) +
    geom_text(aes(y = to, x = from, label = p.signif)) +
    theme_bw() + 
    # scale_color_continuous(values = c('#0072B5','#EE4C97')) +
    # scale_colour_gradient2(low="#0072B5", mid = "white",high="#EE4C97",midpoint = 0)+
    scale_colour_gradient2(high="#F97B72", mid = "white",low="#0072B5",midpoint = 0)+
    ylab('') + xlab('') + labs(color = "pearson's cor",size="pearson's cor")
  p.all=c(p.all,list(p))
}
length(p.all)
hub.genes.TME.corplot=mg_merge_plot(p.all,nrow = 2,ncol = 1,heights = c(8,22),align = 'hv')
hub.genes.TME.corplot

plot1=mg_merge_plot(hub.genes.pcd.corplot,p.all[[1]],
                    nrow = 2,ncol = 1,labels = LETTERS[1:2],
                    heights = c(12,4),align = 'hv')

# figure=mg_merge_plot(hub.genes.pcd.corplot,hub.genes.TME.corplot,nrow = 1,ncol = 2)
figure=mg_merge_plot(plot1,p.all[[2]],nrow = 1,ncol = 2)
figure
savePDF('analysis/5_Hub_Genes/ppi.hub.genes.TME.corplot.pdf',figure,height = 10,width = 16)
###########################
load("raw_datas/Pathways_ssGSEA/gse57065.h.all.ssgsea.RData")
gse57065.h.all.ssgsea=t(gse57065.h.all.ssgsea)
colnames(gse57065.h.all.ssgsea)=gsub("HALLMARK_","",colnames(gse57065.h.all.ssgsea))
dim(gse57065.h.all.ssgsea)

#############
cor.res=psych::corr.test(x=data.frame(t(gse57065.exp[ppi.hub.genes,])),
                         y = data.frame(gse57065.h.all.ssgsea,check.names = F))
df_cor=cor.res$r
df_pval=cor.res$p.adj
########
library(tidyverse)
g = pivot_longer(data=rownames_to_column(data.frame(df_cor,check.names = F),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(data.frame(df_pval,check.names = F)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
all(g$from==gp$rowname & g$to==gp$gene)
g$p.adj = gp$p
head(g)
################## 
df=data.frame(g,check.names = F,stringsAsFactors = F)
cor.R <- pivot_wider(data=g[,c(1:2,3)],names_from ='to',values_from ='cor')
cor.R=data.frame(cor.R)
rownames(cor.R)=cor.R$from
cor.R=cor.R[,-1]
cor.R=as.matrix(cor.R)

cor.Pval <- pivot_wider(data=g[,c(1:2,4)],names_from ='to',values_from ='p.adj')
cor.Pval=data.frame(cor.Pval)
rownames(cor.Pval)=cor.Pval$from
cor.Pval=cor.Pval[,-1]
cor.Pval=as.matrix(cor.Pval)
cor.Pval <- ifelse(cor.Pval < 0.05, '*', '')
cor.Pval[is.na(cor.Pval)] <- ''
####################
df=g
df <- df %>%
  mutate(type = cut(cor, breaks = c(-1, 0, 1),
                    labels = c("negative", "positive")),
         p.signif = cut(p.adj, breaks = c(0,0.0001, 0.001, 0.01, 0.05,1),
                        labels = c("****", "**", "**","*","ns"),
                        right = FALSE, include.lowest = TRUE))
head(df)
write.table(df,file = 'analysis/5_Hub_Genes/ppi.hub.genes.hallmark.cor.res.txt',sep = "\t",col.names = T,row.names = F,quote = F)

p=Heatmap(t(cor.R), 
        name = "Correlation", 
        col = circlize::colorRamp2(c(-0.5, 0, 0.5),
                                   c(pal_d3(alpha =1)(7)[1], 'white', pal_d3(alpha =1)(7)[2])),
        border = T,
        show_column_names = F,
        show_column_dend = F,
        show_row_names = T,
        show_row_dend = F,
        cluster_columns=T,
        cluster_rows=T,
        column_names_side = 'bottom',
        column_names_rot = 60,
        row_names_side = 'right',
        rect_gp = gpar(col = NA, lwd = 1),
        row_names_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(t(cor.Pval)[i, j], x, y)
        })


column.sorted=rownames(cor.R)[column_order(p)]
row.sorted=colnames(cor.R)[row_order(p)]

cor.R_2=cor.R[column.sorted,row.sorted]
cor.Pval_2=cor.Pval[column.sorted,row.sorted]
range(cor.R_2)
pdf('analysis/5_Hub_Genes/ppi.hub.genes.hallmark.cor.heatmap.pdf',height = 8,width = 12)
Heatmap(t(cor.R_2),
        name = "Cor", 
        col = circlize::colorRamp2(c(-1, 0, 1),
                                   c(pal_d3(alpha =1)(7)[1], 'white', pal_d3(alpha =1)(7)[2])),
        border = T,
        show_column_names = T,
        show_column_dend = F,
        show_row_names = T,
        show_row_dend = F,
        cluster_columns=T,
        cluster_rows=T,
        column_names_side = 'bottom',
        column_names_rot = 0,
        column_names_centered=T,
        row_names_side = 'right',
        rect_gp = gpar(col = 'white', lwd = 2),
        row_names_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(t(cor.Pval_2)[i, j], x, y)
        })
dev.off()

