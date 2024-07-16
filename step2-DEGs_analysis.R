dim(gse57065.cli) ## 107*8
table(gse57065.cli$disease)
# healthy  sepsis 
# 25      82
gse57065.degs=mg_limma_DEG(exp = gse57065.exp,group = gse57065.cli$disease,ulab = 'sepsis',dlab = 'healthy')
gse57065.degs$Summary
# 1.2-fold    1.3-fold    1.5-fold    2-fold   
# p<0.05   "2717|3628" "1882|2446" "1106|1252" "416|343"
# p<0.01   "2607|3566" "1849|2428" "1103|1251" "415|343"
# FDR<0.05 "2680|3605" "1869|2439" "1106|1252" "416|343"
# FDR<0.01 "2533|3493" "1822|2408" "1097|1247" "413|343"
gse57065.degs.res=gse57065.degs$DEG
write.table(gse57065.degs.res,file = 'analysis/2_DEGs/gse57065.degs.res.txt',sep = "\t",quote = F)
write.table(gse57065.degs.res,file = 'results/Files/gse57065.degs.res.txt',sep = "\t",quote = F)

library(dplyr)
gse57065.degs.up=gse57065.degs.res %>% filter(adj.P.Val<0.05 & logFC>log2(2)) %>% rownames()
gse57065.degs.dwn=gse57065.degs.res %>% filter(adj.P.Val<0.05 & logFC<(-log2(2))) %>% rownames()
length(gse57065.degs.up) ## 416
length(gse57065.degs.dwn) ## 343

################## NETosis geneset
pcd.genesets=readxl::read_excel('raw_datas/PCD.geneSets.PMID36341760.xlsx')
pcd.genesets=data.frame(pcd.genesets)
pcd.genesets=pcd.genesets[,-1]
head(pcd.genesets)

pcd.list=list()
pcd.genesets.df=c()
for(i in colnames(pcd.genesets)){
  pcd.list[[i]]=as.character(na.omit(pcd.genesets[,i]))
  pcd.genesets.df=rbind(pcd.genesets.df,data.frame(PCD=i,Symbol=pcd.genesets[,i],check.names = F,stringsAsFactors = F))
}

NETosis.genes=pcd.list[['Netotic.cell.death']]
intersect(NETosis.genes,gse57065.degs.up) ## "ELANE" "MPO"   "CAMP"  "PADI4"
intersect(NETosis.genes,gse57065.degs.dwn) ## NULL

degs.volcano=mg_volcano(logfc = gse57065.degs.res$logFC,pvalue = gse57065.degs.res$adj.P.Val,
              colors=c(mycolors[5],'grey',mycolors[4]),
              cutFC = 1,
              cutPvalue = 0.05,
              symbol = rownames(gse57065.degs.res),
              showText = intersect(NETosis.genes,gse57065.degs.up),
              leg='GSE57065',
              legend.pos='tl',
              ylab='-log10(FDR)',
              xlab='log2(FC)')
degs.volcano

############# function enrichment analysis
length(gse57065.degs.up) ## 416
length(gse57065.degs.dwn) ## 343

gse57065.degs.up.enrich.res=mg_clusterProfiler(gse57065.degs.up)
gse57065.degs.dwn.enrich.res=mg_clusterProfiler(gse57065.degs.dwn)

write.table(gse57065.degs.up.enrich.res$GO_BP,file = 'analysis/2_DEGs/gse57065.degs.up.enrich.res.txt',sep = "\t",quote = F)
write.table(gse57065.degs.dwn.enrich.res$GO_BP,file = 'analysis/2_DEGs/gse57065.degs.dwn.enrich.res.txt',sep = "\t",quote = F)

write.table(gse57065.degs.up.enrich.res$GO_BP,file = 'results/Files/gse57065.degs.up.enrich.res.txt',sep = "\t",quote = F)
write.table(gse57065.degs.dwn.enrich.res$GO_BP,file = 'results/Files/gse57065.degs.dwn.enrich.res.txt',sep = "\t",quote = F)

table(gse57065.degs.up.enrich.res$Enrich_tab$DB,gse57065.degs.up.enrich.res$Enrich_tab$FDR<0.05)
# TRUE
# geneontology_Biological_Process  212
# geneontology_Cellular_Component   47
# geneontology_Molecular_Function    9
# pathway_KEGG                       3
table(gse57065.degs.dwn.enrich.res$Enrich_tab$DB,gse57065.degs.dwn.enrich.res$Enrich_tab$FDR<0.05)
# TRUE
# geneontology_Biological_Process  235
# geneontology_Cellular_Component   32
# geneontology_Molecular_Function   16
# pathway_KEGG                      31

theme_custom_3=theme(axis.text.x=element_text(size=10,face="plain", family="Times", colour="black"),
                     axis.text.y=element_text(size=10,face="plain", family="Times", colour="black"))

library(enrichplot)
p.all=list()
for(i in c('KEGG','GO_BP','GO_CC','GO_MF')){
  p=barplot(gse57065.degs.up.enrich.res[[i]],width=0.5,showCategory=20)
  p=p+theme_custom_3+scale_fill_distiller(palette="OrRd")+
    labs(title = 'Up-regulated')+xlab(i)
  p.all=c(p.all,list(p))
}
length(p.all)

up.enrich.boxplot=mg_merge_plot(p.all,nrow = 2,ncol = 2,
                                align = 'hv', labels = 'D')
up.enrich.boxplot=p.all[[2]]

p.all=list()
for(i in c('KEGG','GO_BP','GO_CC','GO_MF')){
  p=barplot(gse57065.degs.dwn.enrich.res[[i]],width=0.5,showCategory=20)
  p=p+theme_custom_3+scale_fill_distiller(palette="OrRd")+
    labs(title = 'Down-regulated')+xlab(i)
  p.all=c(p.all,list(p))
}
length(p.all)

dwn.enrich.boxplot=mg_merge_plot(p.all,nrow = 2,ncol = 2,
                     align = 'hv',labels = 'D')
dwn.enrich.boxplot=p.all[[2]]

plot1=mg_merge_plot(degs.volcano,dwn.enrich.boxplot,
                    nrow = 1, ncol = 2,
                    labels = c('A', 'C'))
plot2=mg_merge_plot(up.enrich.boxplot,mg_getplot_bank('X'),
                    nrow = 1, ncol = 2,widths = c(3,1),
                    labels = 'B')

figure=mg_merge_plot(plot1,plot2,
                     nrow = 2, ncol = 1,
                     labels = c('','B'))
savePDF('analysis/2_DEGs/degs.enrich.boxplot.pdf',figure,height = 14,width = 14)
