############################## WGCNA
##### cluster gene
library(WGCNA)
################# GSE57065
load("raw_datas/GEO/gse57065.exp.RData")
load("raw_datas/GEO/gse57065.cli.RData")
all(colnames(gse57065.exp)==rownames(gse57065.cli))
gse57065.cli_dt=data.frame(gse57065.cli,NETosis=gse57065.pcd.ssgsea[rownames(gse57065.cli),'Netotic.cell.death'])
table(gse57065.cli_dt$disease)
# healthy  sepsis 
# 25      82
dim(gse57065.cli_dt) ## 107*9
####################
dim(gse57065.exp) ## 20549*107
mads=apply(gse57065.exp[,rownames(gse57065.cli_dt)],1,mad)
pcg.selected=names(mads[order(mads,decreasing = T)[1:(length(mads)*0.7)]])
length(pcg.selected) ## 14384
###############
edata.wgcna=gse57065.exp[pcg.selected,rownames(gse57065.cli_dt)]
edata.wgcna=t(edata.wgcna)
dim(edata.wgcna) ## 107 * 14384
range(edata.wgcna)
## 1.334591 14.234079

########
pdf('analysis/3_WGCNA/Fig3ABC.pdf',width = 8,height = 8)
power=mg_wgcna_get_power(edata.wgcna,RsquaredCut=0.9)
dev.off()

power$cutPower ## 7
net=mg_WGCNA_getModule(edata.wgcna,power = power$cutPower
                       , deepSplit = 1, mergeCutHeight = 0.25
                       , minModuleSize = 60)
length(table(net$Modules[,2])) ## 19

pdf('analysis/3_WGCNA/Fig3D.pdf',height = 6,width = 8)
plotDendroAndColors(net$Tree, net$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

fig3e=mg_barplot_point(labels = names(table(net$Modules[,2]))
                       ,values = as.numeric(table(net$Modules[,2]))
                       ,point_sizes = 2
                       ,point_cols = names(table(net$Modules[,2]))
                       ,xlab = 'Number of Genes',legend.pos = NULL)
fig3e=fig3e+theme(axis.text.y = element_text(family = "Times", face = "plain"),
                  axis.title.y = element_text(family = "Times", face = "plain"),
                  panel.background = element_rect(fill = "white",
                                                  colour = "black"),
                  plot.background = element_rect(fill = "white",
                                                 colour = "black")
)
fig3e 
savePDF('analysis/3_WGCNA/Fig3E.pdf',fig3e,height = 5,width = 4)

#### 
# Calculate eigengenes
MEs = net$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('analysis/3_WGCNA/Fig3F.pdf',height = 6,width = 10,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()

############### GSE57065
dim(gse57065.cli_dt)
datTraits0 = data.frame(type=gse57065.cli_dt$disease,NETosis=gse57065.cli_dt$NETosis,gse57065.cli_dt[,c("age",'gender')],stringsAsFactors = F)
datTraits0$type=factor(datTraits0$type,levels = c('healthy','sepsis'))
datTraits0$gender=factor(datTraits0$gender,levels = c("Female","Male"))
str(datTraits0)
######
datTraits=datTraits0
datTraits[,c('type',"gender")]=sapply(datTraits[, c('type',"gender")], function(x)as.numeric(as.factor(x)))
head(datTraits)
head(datTraits0)
dim(datTraits)
####### Calculate module eigengenes
MEs<-net$MEs
dim(MEs)
## Define numbers of genes and samples
nGenes = ncol(edata.wgcna)
nSamples = nrow(edata.wgcna)
## 
modTraitCor = WGCNA::cor(MEs[,rownames(MEDiss)[METree$order]]
                         , datTraits
                         , use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, nSamples)
## 
textMatrix = paste(signif(modTraitCor, 2), "(", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

dev.off()
pdf('analysis/3_WGCNA/Fig3G.pdf',width = 6,height = 6)
labeledHeatmap(Matrix = data.frame(modTraitCor), 
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = rownames(modTraitCor), 
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix), setStdMargins = FALSE,
               cex.text = 0.5, zlim = c(-1,1),xLabelsAngle =45,xLabelsAdj =1,
               main = paste("Module-trait relationships"))
dev.off()

#
#
## 
dim(edata.wgcna) ## 107*14384
dim(net$MEs) ## 107 * 19
geneModuleMembership <- as.data.frame(signedKME(edata.wgcna, data.frame(net$MEs), outputColumnName = ""))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

#
all(rownames(edata.wgcna)==rownames(datTraits))
geneTraitSignificance <- as.data.frame(cor(edata.wgcna, datTraits, use = 'pairwise.complete.obs'))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

###
### Intramodular analysis: identifying genes with high GS and MM
############
table(net$Modules[,2])

rownames(modTraitCor)
#### module-GENE
table(net$Modules[,2])
blue.module.genes=rownames(net$Modules)[net$Modules[,2]=='blue']
length(blue.module.genes) ## 1828

df=net$Modules[,2][net$Modules[,2] %in% c('blue')]
df=data.frame(df)
colnames(df)='Module'
df=df %>% tibble::rownames_to_column('Gene')
dim(df) ## 1828*2
head(df)
write.table(df,file = 'results/Files/gse65682.wgcna.blue.module.genes.txt',sep = '\t',quote = F,row.names = F,col.names = T)
write.table(df,file = 'analysis/3_WGCNA/gse65682.wgcna.blue.module.genes.txt',sep = '\t',quote = F,row.names = F,col.names = T)

########### module-blue
blue.hub.genes=get_module_hub_genes(net = net,module = "blue",trait = "NETosis"
                                         ,output = "analysis/3_WGCNA/wgcna.blue.scatterplot.pdf"
                                         ,MM=0.4,GS=0.4,pval=0.05)
length(blue.hub.genes) ## 
#########################
modNames = substring(names(MEs), 3)
moduleColors = unname(net$Modules[,2])
modNames
moduleColors

module = "blue"
column = match(module, modNames)
column
moduleGenes = moduleColors==module
trait="NETosis"
table(moduleColors)[module]
colnames(geneModuleMembership)[column]

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, trait]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for ",trait),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
table(moduleColors)[module]

######################### blue module
dev.off()
pdf('analysis/3_WGCNA/gse57065.module.degs.venn.pdf',height = 7,width = 7)
mg_venn_plot(list(Up=gse57065.degs.up,Blue=blue.hub.genes,Down=gse57065.degs.dwn),fill=c(mycolors[4],'blue',mycolors[5]))
dev.off()

###################
wgcna.degs.genes=Reduce(intersect,list(Blue=blue.hub.genes,DEGs=c(gse57065.degs.up,gse57065.degs.dwn)))
length(wgcna.degs.genes) ## 231
write.table(data.frame(wgcna.degs.genes),'analysis/4_PPI/wgcna.degs.genes.txt',sep = "\t",row.names = F,col.names = F,quote = F)
