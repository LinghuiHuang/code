get_PlotMutiBoxplot=function(exp.immune.score,group,group.val=NULL,test_method='kruskal.test',ylab='Estimated Proportion',legend.pos='tr',group_cols=NULL,xangle=45){
  if(is.null(group_cols)){
    p=mg_PlotMutiBoxplot(exp.immune.score[rownames(group),]
                         , group = group[,group.val]
                         , legend.pos = legend.pos
                         , add = 'boxplot'
                         , xangle=xangle
                         , ylab = ylab
                         , group_cols = pal_nejm()(8)[c(1,3,2,4)]
                         , test_method = test_method
    )
  }else{
    p=mg_PlotMutiBoxplot(exp.immune.score[rownames(group),]
                         , group = group[,group.val]
                         , legend.pos = legend.pos
                         , add = 'boxplot'
                         , xangle=xangle
                         , binwidth=0.02
                         , ylab = ylab
                         , group_cols = group_cols
                         , test_method = test_method
    )
  }
  
  return(p)
}

mg_violin_use=function(data,group.col=NULL,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL){
  library(ggplot2)
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }else if(legend.pos=='none'){
    pos='none'
  }
  
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  if(!is.null(ylim)){
    data_m=data_m[data_m[,2]<=max(ylim),]
    data_m=data_m[which(!is.na(data_m[,1])),]
    ylim[2]=1.2*ylim[2]
  }
  ######## violin
  p1<-ggplot(data_m,aes(x=Group,y=value,color=Group,fill=Group))+
    geom_violin(width=0.5, alpha=0.7)
  
  if(ct<=2){
    if(!is.null(group.col)){
      p1=p1+scale_color_manual(values = group.col)
      p1=p1+scale_fill_manual(values = group.col)
    }else{
      p1=p1+scale_color_manual(values=ggsci::pal_d3("category20", alpha = 0.6)(9))
      p1=p1+scale_fill_manual(values=ggsci::pal_d3("category20", alpha = 0.6)(9))
    }
  }else if(ct<=10){
    # p1=p1+ggsci::scale_fill_npg(name=leg.title)
    if(!is.null(group.col)){
      p1=p1+scale_color_manual(values = group.col)
      p1=p1+scale_fill_manual(values = group.col)
    }else{
      p1=p1+scale_color_manual(values = c(pal_lancet('lanonc',alpha =0.8)(9)[c(2,1,4,3,5:6)]))
      p1=p1+scale_fill_manual(values = c(pal_lancet('lanonc',alpha =0.8)(9)[c(2,1,4,3,5:6)]))
    }
    
  }else if(ct<=20){
    p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  }else if(ct<=30){
    cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }else if(ct<=38){
    cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10)
                ,ggsci::pal_d3("category20", alpha = 0.6)(20)
                ,ggsci::pal_nejm("default", alpha = 0.6)(8))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }
  
  ####### boxplot
  p1=p1+geom_boxplot(width=0.2,fill="white",outlier.shape = NA)
  #######
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.6,shape=21, size=1,show.legend=FALSE,width = 0.15)
    }else{
      p1<-p1+geom_jitter(alpha=0.6,shape=21, size=1,show.legend=FALSE,width = 0.15,size=point_size)
    }
  }
  
  p1=p1+theme_bw()
  p1=p1+theme(axis.text.x=tx, #
              axis.text.y=element_text(family="Times",face="plain"), #
              axis.title.y=element_text(family="Times",face="plain"), #
              #panel.border = element_blank(),axis.line = element_line(colour = "black"), #
              legend.text=element_text(face="plain", family="Times", colour="black"  #
              ),
              legend.title=element_text(face="plain", family="Times", colour="black" #
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
              ,panel.grid.major = element_blank(),   #
              panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- anova(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til)
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    if(length(comps)<7){
      p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
    }
  }
  return(p1)
}

ggplotKMCox=function(dat,title='Groups',labs=NULL,add_text=NULL,pal=NULL){
  library(ggplot2)
  library(ggsci)
  library(survival)
  library(ggpubr)
  library(survminer)
  ######
  colnames(dat)=c('time','status','groups')
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  if(is.null(pal)){
    pal=pal_lancet()(9)[c(2,4,3,1,5:6,9)]
  }
  surv=survminer::ggsurvplot(sf, data = dat
                             , palette = pal 
                             , pval = TRUE
                             , surv.median.line = 'hv'
                             , conf.int = T
                             # , linetype = "strata"
                             , xlab = "Time(years)"
                             , conf.int.style = 'step'
                             , pval.coord = c(0, 0.2)#Add p-value
                             , risk.table = TRUE
                             , ggtheme = theme_pubr(base_size=12,base_family='Times')
                             , font.family='Times'
                             , risk.table.y.text = FALSE
                             , legend.title = title
                             , legend.labs = labs)
  p1=surv$plot
  p2=surv$table
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(0.9,0.3),align = "v")
  return(g2)
}


get.IOBR.immu.format=function(tcga.t.exp.cibersort){
  tcga.t.exp.cibersort = data.frame(tcga.t.exp.cibersort)
  rownames(tcga.t.exp.cibersort) = tcga.t.exp.cibersort$ID
  tcga.t.exp.cibersort = tcga.t.exp.cibersort[, -1]
  colnames(tcga.t.exp.cibersort) = gsub('(.*)_.*', "\\1", colnames(tcga.t.exp.cibersort))
  return(tcga.t.exp.cibersort)
}
###########################
mycolors=c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#F97B72","#4B4B8F","#A5AA99")
mycolors=mycolors[-1]

################ GSE57065
load("raw_datas/GEO/gse57065.exp.RData")
load("raw_datas/GEO/gse57065.cli.RData")
all(colnames(gse57065.exp)==gse57065.cli$Acc)
################ GSE65682
load("raw_datas/GEO/gse65682.exp.RData")
load("raw_datas/GEO/gse65682.cli.RData")
all(colnames(gse65682.exp)==gse65682.cli$Acc)
################ GSE145227
load("raw_datas/GEO/gse145227.exp.RData")
load("raw_datas/GEO/gse145227.cli.RData")
all(colnames(gse145227.exp)==gse145227.cli$Acc)
################ GSE54514
load("raw_datas/GEO/gse54514.exp.RData")
load("raw_datas/GEO/gse54514.cli.RData")
all(colnames(gse54514.exp)==gse54514.cli$Acc)
################ GSE95233
load("raw_datas/GEO/gse95233.exp.RData")
load("raw_datas/GEO/gse95233.cli.RData")
all(colnames(gse95233.exp)==gse95233.cli$Acc)
######################################################################################## 
################## Step1
######################################################################################## 
library(ggplot2)
theme_custom_1=theme_classic()+
  theme(axis.text.x= element_blank(),
        # axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(r = 0, l = 0),
        legend.position = "none"
  )

theme_custom_2=theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      family = "Times"
    ),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  )
###############
load("raw_datas/Pathways_ssGSEA/gse57065.pcd.ssgsea.RData")
gse57065.pcd.ssgsea=t(gse57065.pcd.ssgsea)
gse57065.pcd.ssgsea=data.frame(gse57065.pcd.ssgsea)
dim(gse57065.pcd.ssgsea) ## 107*12
head(gse57065.cli)
gse57065.cli$`collection time`=factor(gse57065.cli$`collection time`,levels = c("healthy","0 hr","24 hr","48 hr"))
gse57065.cli$sapsii=factor(gse57065.cli$sapsii,levels = c('Low','High'))
range(gse57065.cli$age) ## adult
dim(gse57065.cli)
head(gse57065.cli)

dt=cbind(gse57065.cli,gse57065.pcd.ssgsea[rownames(gse57065.cli),])
write.table(dt,file = 'analysis/1_PCD_Clinical/gse57065.pcd.ssgsea.txt',sep = "\t",quote = F)

gse57065.cli.pcd.boxplot=get_PlotMutiBoxplot(gse57065.pcd.ssgsea[rownames(gse57065.cli),],gse57065.cli
                                             , ylab = 'ssGSEA score'
                                             , group_cols = mycolors
                                             , legend.pos = 'bl'
                                             , group.val = 'disease', xangle = 45) + labs(color = 'Group')

gse57065.cli.pcd.boxplot

###
table(gse57065.cli$disease,gse57065.cli$`collection time`)
#           healthy 0 hr 24 hr 48 hr
# healthy      25    0     0     0
# sepsis        0   28    28    26
inds=which(gse57065.cli$disease=='sepsis')
length(inds) ## 82
p2=get_PlotMutiBoxplot(gse57065.pcd.ssgsea[rownames(gse57065.cli)[inds],],gse57065.cli[inds,]
                       , ylab = 'ssGSEA score'
                       , group_cols = mycolors
                       , legend.pos = 'top'
                       , group.val = 'collection time', xangle = 45) + labs(color = 'Collection time')
p2

inds=which(gse57065.cli$disease=='sepsis')
length(inds) ## 82
table(gse57065.cli$sapsii)
# Low High 
# 42   40 
p3=get_PlotMutiBoxplot(gse57065.pcd.ssgsea[rownames(gse57065.cli)[inds],],gse57065.cli[inds,]
                       , ylab = 'ssGSEA score'
                       , group_cols = mycolors
                       , legend.pos = 'top'
                       , group.val = 'sapsii', xangle = 45) + labs(color = 'SAPSII')

p3

table(gse57065.cli$disease,gse57065.cli$sapsii)
#         Low High
# healthy   0    0
# sepsis   42   40

inds=which(gse57065.cli$sapsii=='Low')
length(inds) ## 42
get_PlotMutiBoxplot(gse57065.pcd.ssgsea[inds,],gse57065.cli[inds,]
                    , ylab = 'ssGSEA score'
                    , group_cols = mycolors
                    , legend.pos = 'top'
                    , group.val = 'collection time', xangle = 45) + labs(color = 'SAPSII-Low')

dt=data.frame(`collection time`=as.character(gse57065.cli[inds,"collection time"]),NETosis=gse57065.pcd.ssgsea[inds,'Netotic.cell.death'])
table(dt$collection.time)
mg_violin_use(data = dt
                ,melt = T
                ,group.col = mycolors
                ,xlab = "collection time"
                ,ylab = "NETosis"
                ,jitter=F
                ,test_method = 'kruskal.test'
                ,cmp_test_method = 'wilcox.test'
                ,legend.pos = NULL
                ,show_compare = T)

########################
inds=which(gse57065.cli$sapsii=='High')
length(inds) ## 40
get_PlotMutiBoxplot(gse57065.pcd.ssgsea[inds,],gse57065.cli[inds,]
                    , ylab = 'ssGSEA score'
                    , group_cols = mycolors
                    , legend.pos = 'top'
                    , group.val = 'collection time', xangle = 45) + labs(color = 'SAPSII-High')

dt=data.frame(`collection time`=as.character(gse57065.cli[inds,"collection time"]),NETosis=gse57065.pcd.ssgsea[inds,'Netotic.cell.death'])
table(dt$collection.time)
mg_violin_use(data = dt
              ,melt = T
              ,group.col = mycolors
              ,xlab = "collection time"
              ,ylab = "NETosis"
              ,jitter=F
              ,test_method = 'kruskal.test'
              ,cmp_test_method = 'wilcox.test'
              ,legend.pos = NULL
              ,show_compare = T)

########################## GSE95233
load("raw_datas/Pathways_ssGSEA/gse95233.pcd.ssgsea.RData")
gse95233.pcd.ssgsea=t(gse95233.pcd.ssgsea)
gse95233.pcd.ssgsea=data.frame(gse95233.pcd.ssgsea,check.names = F,stringsAsFactors = F)
dim(gse95233.pcd.ssgsea)
dt=cbind(gse95233.cli,gse95233.pcd.ssgsea[rownames(gse95233.cli),])
write.table(dt,file = 'analysis/1_PCD_Clinical/gse95233.pcd.ssgsea.txt',sep = "\t",quote = F)

head(gse95233.cli)
table(gse95233.cli$disease)
# healthy  sepsis 
# 42     760
gse95233.cli.pcd.boxplot=get_PlotMutiBoxplot(gse95233.pcd.ssgsea[rownames(gse95233.cli),],gse95233.cli
                                         , ylab = 'ssGSEA score'
                                         , group_cols = mycolors
                                         , legend.pos = 'bl'
                                         , group.val = 'disease', xangle = 45) + labs(color = 'Group')

gse95233.cli.pcd.boxplot
########################## GSE65682
load("raw_datas/Pathways_ssGSEA/gse65682.pcd.ssgsea.RData")
gse65682.pcd.ssgsea=t(gse65682.pcd.ssgsea)
gse65682.pcd.ssgsea=data.frame(gse65682.pcd.ssgsea,check.names = F,stringsAsFactors = F)
dim(gse65682.pcd.ssgsea)
dt=cbind(gse65682.cli,gse65682.pcd.ssgsea[rownames(gse65682.cli),])
write.table(dt,file = 'analysis/1_PCD_Clinical/gse65682.pcd.ssgsea.txt',sep = "\t",quote = F)

# dt=gse65682.cli %>% arrange(desc(Title),desc(`pneumonia diagnoses`),desc(thrombocytopenia),desc(abdominal_sepsis_and_controls))
# View(dt)

head(gse65682.cli)
table(gse65682.cli$disease)
# healthy  sepsis 
# 42     760
gse65682.disease.pcd.boxplot=get_PlotMutiBoxplot(gse65682.pcd.ssgsea[rownames(gse65682.cli),],gse65682.cli
                                                 , ylab = 'ssGSEA score'
                                                 , group_cols = mycolors
                                                 , legend.pos = 'bl'
                                                 , group.val = 'disease', xangle = 45) + labs(color = 'Group')

head(gse65682.cli)
# 
table(gse65682.cli$endotype_class)
inds=which(!is.na(gse65682.cli$endotype_class))
gse65682.Mars.pcd.boxplot=get_PlotMutiBoxplot(gse65682.pcd.ssgsea[inds,],gse65682.cli[inds,]
                                                 , ylab = 'ssGSEA score'
                                                 , group_cols = mycolors
                                                 , legend.pos = 'bl'
                                                 , group.val = 'endotype_class', xangle = 45) + labs(color = 'endotype_class')
############# Mars1
table(gse65682.cli$endotype_class,gse65682.cli$survival)
inds=which(gse65682.cli$endotype_class=='Mars1')
gse65682.Mars1.status.pcd.boxplot=get_PlotMutiBoxplot(gse65682.pcd.ssgsea[inds,],gse65682.cli[inds,]
                                               , ylab = 'ssGSEA score'
                                               , group_cols = mycolors
                                               , legend.pos = 'bl'
                                               , group.val = 'survival', xangle =  45) + labs(color = 'survival')
gse65682.Mars1.status.pcd.boxplot
############ Mars2-4
inds=which(gse65682.cli$endotype_class %in% c('Mars2','Mars3','Mars4'))
gse65682.Mars24.status.pcd.boxplot=get_PlotMutiBoxplot(gse65682.pcd.ssgsea[inds,],gse65682.cli[inds,]
                                                     , ylab = 'ssGSEA score'
                                                     , group_cols = mycolors
                                                     , legend.pos = 'bl'
                                                     , group.val = 'survival', xangle = 45) + labs(color = 'survival')

################ GSE185263
load("raw_datas/Pathways_ssGSEA/gse185263.pcd.ssgsea.RData")
gse185263.pcd.ssgsea=t(gse185263.pcd.ssgsea)
gse185263.pcd.ssgsea=data.frame(gse185263.pcd.ssgsea,check.names = F,stringsAsFactors = F)
dim(gse185263.pcd.ssgsea)

dt=cbind(gse185263.cli,gse185263.pcd.ssgsea[rownames(gse185263.cli),])
write.table(dt,file = 'analysis/1_PCD_Clinical/gse185263.pcd.ssgsea.txt',sep = "\t",quote = F)

head(gse185263.cli)
table(gse185263.cli$disease)
# healthy  sepsis 
# 44     348
gse185263.disease.pcd.boxplot=get_PlotMutiBoxplot(gse185263.pcd.ssgsea[rownames(gse185263.cli),],gse185263.cli
                                                  , ylab = 'ssGSEA score'
                                                  , group_cols = mycolors
                                                  , legend.pos = 'bl'
                                                  , group.val = 'disease', xangle = 45) + labs(color = 'Group')
gse185263.disease.pcd.boxplot
#######
table(gse185263.cli$disease,gse185263.cli$survival)
inds=which(gse185263.cli$disease=='sepsis' & !is.na(gse185263.cli$survival))
gse185263.status.pcd.boxplot=get_PlotMutiBoxplot(gse185263.pcd.ssgsea[inds,],gse185263.cli[inds,]
                       , ylab = 'ssGSEA score'
                       , group_cols = mycolors
                       , legend.pos = 'bl'
                       , group.val = 'survival', xangle = 45) + labs(color = 'survival')
gse185263.status.pcd.boxplot

library(ggpubr)
bulk.cli.pcd.boxplot=patchwork::wrap_plots(list(gse57065.cli.pcd.boxplot+rremove("x.text"),
                                                gse95233.cli.pcd.boxplot + rremove("x.text"),
                                                gse185263.disease.pcd.boxplot + rremove("x.text"),
                                                gse185263.status.pcd.boxplot,
                                                gse65682.disease.pcd.boxplot + rremove('x.text'),
                                                gse65682.Mars.pcd.boxplot + rremove('x.text'),
                                                gse65682.Mars1.status.pcd.boxplot + rremove('x.text'),
                                                gse65682.Mars24.status.pcd.boxplot), ncol = 2, nrow = 4, byrow = F)
savePDF('analysis/1_PCD_Clinical/bulk.cli.pcd.boxplot.pdf',bulk.cli.pcd.boxplot,height = 16,width = 12)
############################## KM analysis
dt=data.frame(OS.time=gse65682.cli$OS.time,OS=gse65682.cli$OS,NETosis=gse65682.pcd.ssgsea$Netotic.cell.death,stringsAsFactors = F,check.names = F)
inds=which(gse65682.cli$disease=='sepsis')
dt=dt[inds,]
dt$group=ifelse(dt$NETosis>median(dt$NETosis),'High','Low')
head(dt)

NETosis.color=mycolors[c(4,3)]
gse65682.NETosis.km=ggplotKMCox(data.frame(time = dt$OS.time
                                           , event = dt$OS
                                           , groups = dt$group)
                                , pal = NETosis.color
                                , title = 'NETosis'
                                , labs = c('High', 'Low')
                                , add_text = '') 
gse65682.NETosis.km
savePDF('analysis/1_PCD_Clinical/gse65682.NETosis.km.pdf',gse65682.NETosis.km,height = 5,width = 5)

################### GSE95233
dt=data.frame(OS.time=gse95233.cli$OS.time,OS=gse95233.cli$OS,NETosis=gse95233.pcd.ssgsea$Netotic.cell.death,stringsAsFactors = F,check.names = F)
dt=dt[which(!is.na(dt$OS)),]
res.cut=survminer::surv_cutpoint(data = dt,time = 'OS.time',event = 'OS',variables = 'NETosis')
res.cut$cutpoint$cutpoint

dt$group=ifelse(dt$NETosis>res.cut$cutpoint$cutpoint,'High','Low')
head(dt)

gse95233.NETosis.km=ggplotKMCox(data.frame(time = dt$OS.time
                                           , event = dt$OS
                                           , groups = dt$group)
                                , pal = NETosis.color
                                , title = 'NETosis'
                                , labs = c('High', 'Low')
                                , add_text = '') 
gse95233.NETosis.km
savePDF('analysis/1_PCD_Clinical/gse95233.NETosis.km.pdf',gse95233.NETosis.km,height = 5,width = 5)
