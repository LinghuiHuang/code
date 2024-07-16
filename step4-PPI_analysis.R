############# PPI analysis
string.node=read.csv('analysis/4_PPI/string_interactions_short.tsv default node.csv',check.names = F,stringsAsFactors = F)
nodes=string.node$name
length(nodes) ## 163
########
table(string.node$`MCODE::Clusters (1)`)
#             Cluster 1 Cluster 2 Cluster 3 Cluster 4 Cluster 5 Cluster 6 Cluster 7 Cluster 8 
# 107        16        19         4         4         4         3         3         3

mcode.cluster1=string.node$name[which(string.node$`MCODE::Clusters (1)`=='Cluster 1')]
mcode.cluster2=string.node$name[which(string.node$`MCODE::Clusters (1)`=='Cluster 2')]
length(mcode.cluster1) ## 16
length(mcode.cluster2) ## 19

ppi.nodes.enrich.res=mg_clusterProfiler(nodes)
mcode.cluster1.enrich.res=mg_clusterProfiler(mcode.cluster1)
mcode.cluster2.enrich.res=mg_clusterProfiler(mcode.cluster2)

write.table(ppi.nodes.enrich.res$Enrich_tab,file = 'analysis/4_PPI/ppi.nodes.enrich.res.txt',sep = "\t",quote = F)
write.table(ppi.nodes.enrich.res$Enrich_tab,file = 'results/Files/ppi.nodes.enrich.res.txt',sep = "\t",quote = F)

library(enrichplot)
p.all=list()
for(i in c('KEGG','GO_BP','GO_CC','GO_MF')){
  p=barplot(ppi.nodes.enrich.res[[i]],width=0.5,showCategory=20)
  p=p+theme_custom_3+scale_fill_distiller(palette="PuRd")+
    labs(title = i)
  p.all=c(p.all,list(p))
}
length(p.all)

ppi.nodes.enrich.boxplot=mg_merge_plot(p.all,nrow = 2,ncol = 2,
                                align = 'hv', labels = 'D')
ppi.nodes.enrich.boxplot
savePDF('analysis/4_PPI/ppi.nodes.enrich.boxplot.pdf',ppi.nodes.enrich.boxplot,height = 8,width = 14)

########### top10
dt=string.node %>% dplyr::arrange(desc(Degree),desc(BetweennessCentrality),desc(ClosenessCentrality))
top10=dt[1:10,] %>% dplyr::select(name,Degree,BetweennessCentrality,ClosenessCentrality)
top10
dim(top10)
top10[,c(3,4)]=apply(top10[,c(3,4)], 2, function(x){round(x,digits =3)})

write.table(top10,file = 'analysis/4_PPI/top10.txt',sep = "\t",quote = F,row.names = T,col.names = T)

dev.off()
pdf("analysis/4_PPI/MCODE.venn.pdf",height = 6,width = 6)
mg_venn_plot(list(Top10=top10$name,
                  MCODE.Clusters1=mcode.cluster1,
                  MCODE.Clusters2=mcode.cluster2),fill = mycolors[1:3])
dev.off()  
