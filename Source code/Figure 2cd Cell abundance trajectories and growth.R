rm(list=ls()); require(data.table);require(dplyr);require(ggplot2);require(tidyr);
require(mgcv)
library(gratia)
require(scales)
require(lme4);require(lmerTest)

# Specify source data location 
fileloc <- ("~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/")
setwd(paste0(fileloc,"/cellAbundance/"))

# load time series of abundance for CAMA-1 and MCF-7 
mcf7dd <- data.table(read.csv(file="MCF7 U01 Cell Number afatinib ribo combo timecourse.csv"))
mcf7dd[,CellLine:="MCF-7"]
camadd <- data.table(read.csv(file="CAMA1 U01 Cell Number afatinib ribo combo timecourse.csv"))

dd<- rbind(
  mcf7dd%>%select(intersect(names(mcf7dd),names(camadd))),
  camadd%>%select(intersect(names(mcf7dd),names(camadd))) )

# annotate resistance state and order
dd[, CelltypeLab:="Resistant"]
dd[Celltype=="Sens", CelltypeLab:="Sensitive" ]
dd$CelltypeLab <- factor(dd$CelltypeLab  , levels=c("Sensitive" ,"Resistant" ))

# annotate cell line and order
dd[, CellLine_typeLab:=paste0(CellLine," ",CelltypeLab)]
dd$CellLine_typeLab <- factor(dd$CellLine_typeLab  , levels=c( "CAMA-1 Sensitive", "CAMA-1 Resistant", "MCF-7 Sensitive","MCF-7 Resistant" ))
dd$Treatment <- factor(dd$Treatment,levels=c("DMSO", "Ribociclib","Afatinib",   "Combo"   ))

# compare cell number to the mean final cell number under DMSO for each cell line
dd[, meanPredictedCellNumber:= sum(PredictedCellNumber*(Day==21)*(Treatment=="DMSO"))/sum((Day==21)*(Treatment=="DMSO")) , by=c("CellLine_typeLab")]
dd[, NormCellNumber:=PredictedCellNumber/meanPredictedCellNumber]


options(scipen=0)
ggplot(dd , aes(y=(NormCellNumber), x=Day)) + 
  geom_point(aes(col=Treatment,fill=Treatment),size=2.5)+
  geom_path(aes(col=Treatment,fill=Treatment,group=interaction(Rep, Treatment )))  +
  facet_wrap(~CellLine_typeLab, nrow=1)+#, scale="free")  +
  theme_classic(base_size=26)+
  labs(y="Cell Number \n (relative to DMSO endpoint)")+theme(aspect.ratio = 1)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
  )+
  #scale_y_continuous(labels = scientific)+
  theme(panel.border=element_blank())+
  scale_color_discrete(labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  scale_fill_discrete(labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))

SvLoc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/"
#ggsave(file=paste0(SvLoc,"Fig1/Fig1E Cell Number timecourse.pdf"), height=13, width=18, dpi=320)
#ggsave(file=paste0(SvLoc,"Fig1/Fig1E Cell Number timecourse.png"), height=13, width=18, dpi=320)

# Generate GAM/GP prediction output data.table
preds<-dd
preds[, isRiboTreated:=0]
preds[Treatment%in%c("Ribociclib","Combo"), isRiboTreated:=1]
preds[, isAfatTreated:=0]
preds[Treatment%in%c("Afatinib","Combo"), isAfatTreated:=1]
preds[, isCombTreated:=0]
preds[Treatment=="Combo", isCombTreated:=1]
preds[, isDMSO:=0]
preds[Treatment=="DMSO", isDMSO:=1]


# Contrast model variants to characterize time series
nknot<-4
gam1<-gam(log10(NormCellNumber)~ -1+CellLine_typeLab+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*(Treatment =="DMSO")),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*(Treatment =="Ribociclib")),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*(Treatment =="Afatinib")),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*(Treatment =="Combo")),k=nknot)+
           
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*(Treatment =="DMSO")),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*(Treatment =="Ribociclib")),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*(Treatment =="Afatinib")),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*(Treatment =="Combo")),k=nknot)+
            
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*(Treatment =="DMSO")),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*(Treatment =="Ribociclib")),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*(Treatment =="Afatinib")),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*(Treatment =="Combo")),k=nknot)+
            
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*(Treatment =="DMSO")),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*(Treatment =="Ribociclib")),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*(Treatment =="Afatinib")),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*(Treatment =="Combo")),k=nknot) ,
          data=preds)

gam2<-gam(log10(NormCellNumber)~ -1+CellLine_typeLab+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*isRiboTreated),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*isAfatTreated),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*isCombTreated),k=nknot)+
          
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*isRiboTreated),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*isAfatTreated),k=nknot)+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*isCombTreated),k=nknot)+
            
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*isRiboTreated),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*isAfatTreated),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*isCombTreated),k=nknot)+

            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*isRiboTreated),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*isAfatTreated),k=nknot)+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*isCombTreated),k=nknot) ,
            data=preds)

gam3 <- gam(log10(NormCellNumber)~ -1+CellLine_typeLab+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*isRiboTreated), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*isAfatTreated), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="CAMA-1 Sensitive")*Day*isCombTreated), k=nknot, bs="gp")+
           
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*isRiboTreated), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*isAfatTreated), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="CAMA-1 Resistant")*Day*isCombTreated), k=nknot, bs="gp")+
            
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*isRiboTreated), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*isAfatTreated), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="MCF-7 Sensitive")*Day*isCombTreated), k=nknot, bs="gp")+
            
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*isRiboTreated), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*isAfatTreated), k=nknot, bs="gp")+
            s(I((CellLine_typeLab=="MCF-7 Resistant")*Day*isCombTreated), k=nknot , bs="gp") ,
          data=preds)

gam4 <- gam(log10(NormCellNumber)~ -1+CellLine_typeLab+
               s(I(Day),by=Treatment, k=nknot, bs="gp"),
            data=preds,Method=REML,select = F)

gam5 <- gam(log10(NormCellNumber)~ -1+CellLine_typeLab+
              s(I(Day), k=nknot, bs="gp")+
              s(I(Day*isRiboTreated), k=nknot, bs="gp")+
              s(I(Day*isAfatTreated), k=nknot, bs="gp")+
              s(I(Day*isCombTreated), k=nknot, bs="gp"),
            data=preds,Method=REML,select = F)
  

# model comparison
AIC(gam1,gam2,gam3,gam4,gam5)         

# model statistical summaries
summary(gam3);
summary(gam4)
summary(gam5)  

# correlation of predictions
plot(predict(gam3),predict(gam4))
plot(predict(gam3),predict(gam4))
plot(predict(gam3),predict(gam5))

# Add model predictions and uncertainty 
preds$preds<-10^(predict(gam3))
predlong<-data.table( expand.grid(CellLine_typeLab=unique(preds$CellLine_typeLab) ,Treatment=unique(preds$Treatment)  ,Day=seq(min(preds$Day),max(preds$Day))))
predlong[, isRiboTreated:=0]
predlong[Treatment%in%c("Ribociclib","Combo"), isRiboTreated:=1]
predlong[, isAfatTreated:=0]
predlong[Treatment%in%c("Afatinib","Combo"), isAfatTreated:=1]
predlong[, isCombTreated:=0]
predlong[Treatment=="Combo", isCombTreated:=1]

predlong$preds<-10^(predict(gam3,newdata=predlong))
predlong$ucl<-10^(predict(gam3,newdata=predlong)+predict(gam3,newdata=predlong, se.fit=T)$se.fit)
predlong$lcl<-10^(predict(gam3,newdata=predlong)-predict(gam3,newdata=predlong, se.fit=T)$se.fit)


# AUC metric of effect
smooth_valsR <- smooth_estimates(gam5, "s(I(Day * isRiboTreated))")
integralR <- sum(smooth_valsR$.estimate * diff(smooth_valsR$`I(Day * isRiboTreated)`)[1] )
smooth_valsA <- smooth_estimates(gam5, "s(I(Day * isAfatTreated))")
integralA <- sum(smooth_valsA$.estimate * diff(smooth_valsA$`I(Day * isAfatTreated)`)[1] )
smooth_valsC <- smooth_estimates(gam5, "s(I(Day * isCombTreated))")
integralC <- sum(smooth_valsC$.estimate * diff(smooth_valsC$`I(Day * isCombTreated)`)[1] )

# visualize
ggplot(preds, aes(y= NormCellNumber, x=Day, col=Treatment, fill=Treatment, group=Treatment) )+
  geom_point(size=3)+
  geom_path(data=predlong, aes(y=preds),linewidth=1.5)+
  geom_ribbon(data=predlong, aes(y=preds, ymax=ucl, ymin=lcl), col=NA, alpha=0.3)+
  theme_classic(base_size=26)+
  theme(aspect.ratio = 1)+
  facet_wrap(~CellLine_typeLab,nrow=1)+ # scale_y_continuous(breaks=log10(1000*c(0.5,1,2,4,8,16,32,64,128)), labels=1000*c(0.5,1,2,4,8,16,32,64,128))+
  labs(y="Cell number \n (relative to DMSO endpoint)", x="Day")+
  scale_color_discrete(name="Treatment",labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  scale_fill_discrete(name="Treatment" ,labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black") )

#ggsave(file=paste0(SvLoc,"Fig1/Fig1C Cell Number trajectories.pdf"), height=12, width=18, dpi=320)
#ggsave(file=paste0(SvLoc,"Fig1/Fig1C Cell Number trajectories.png"), height=12, width=18, dpi=320)
ggplot(preds, aes(y= NormCellNumber, x=Day, col=Treatment, fill=Treatment, group=Treatment) )+
  geom_point(size=3)+
  geom_path(data=predlong, aes(y=preds),linewidth=1.5)+
  geom_ribbon(data=predlong, aes(y=preds, ymax=ucl, ymin=lcl), col=NA, alpha=0.3)+
  theme_classic(base_size=26)+
  theme(aspect.ratio = 1)+
  facet_wrap(~CellLine_typeLab,nrow=1)+ # scale_y_continuous(breaks=log10(1000*c(0.5,1,2,4,8,16,32,64,128)), labels=1000*c(0.5,1,2,4,8,16,32,64,128))+
  labs(y="Cell Number \n (relative to DMSO endpoint)", x="Day")+
  scale_color_discrete(name="Treatment",labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  scale_fill_discrete(name="Treatment" ,labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  theme(text = element_blank(),legend.text = element_blank(),legend.title = element_blank()  )
#ggsave(file=paste0(SvLoc,"Fig1/BLANK_Fig1C Cell Number trajectories.pdf"), height=12, width=18, dpi=320)
#ggsave(file=paste0(SvLoc,"Fig1/BLANK_Fig1C Cell Number trajectories.png"), height=12, width=18, dpi=320)


### Growth rate analysis
# duplicate data.table and determine intial cell number at day 0
preds2<-data.table(preds)
preds2[,Measurement0:=sum(PredictedCellNumber*(Day==0)), by=c("Rep", "CellLine_typeLab", "Treatment")]

# calculate relative growth rate (rgr)
preds2[, rgr:=(log(PredictedCellNumber)-log(Measurement0))/(Day-0)]

ggplot(preds2[Day==max(Day)] , aes(y= rgr,x=CellLine_typeLab, col=Treatment,fill=Treatment,group=Treatment ) )+
  geom_point(position=position_dodge(width=0.75),aes(group=interaction(CellLine_typeLab,Treatment)),size=3)+
  geom_boxplot(position=position_dodge(),aes(group=interaction(CellLine_typeLab,Treatment)),alpha=0.4)+
  theme_classic(base_size=26)+
  theme(aspect.ratio = 1)+
  scale_color_discrete(name="Treatment")+scale_fill_discrete(name="Treatment")+
  labs(y="Cancer growth rate \n (cells/cell/hour)",x="Cell line")

# calcaulte median rgr across replicates under DMSO for each cell line
plotnorm<- data.table( preds2[Day==max(Day)] %>%group_by(Treatment,CellLine_typeLab)%>%mutate(medCellNumber=median(PredictedCellNumber),medrgr= median(rgr)) )
plotnorm[, meanrgrNoTreat:=sum(medrgr*(Treatment=="DMSO"))/sum(Treatment=="DMSO"), by=c("CellLine_typeLab")]
# calculate change in rgr relative to DMSO and the median across replicates
plotnorm[,rgrdiff:=(rgr) - (meanrgrNoTreat)]
tmp <- plotnorm[(Treatment!="DMSO") ] %>% group_by(CellLine_typeLab) %>% dplyr::summarise(rgrdiffDMSO= median(rgrdiff))
plotnorm1 <- merge(plotnorm,tmp, by="CellLine_typeLab")#[Treatment!="DMSO"]

# Calculate the additive expectation of the two monotherapies
expect<-plotnorm[!Treatment%in%c("DMSO","Combo") ]%>%group_by(CellLine_typeLab,Treatment)%>%summarise(mueff=mean(rgr-meanrgrNoTreat))%>%group_by(CellLine_typeLab)%>%summarise(expect=sum(mueff))
expect<-data.table(plotnorm%>%
                     group_by(CellLine_typeLab, Treatment, meanrgrNoTreat)%>%
                     summarise(mueff=mean(rgr))) [!Treatment%in%c("DMSO","Combo") ]%>%
  group_by(CellLine_typeLab,meanrgrNoTreat)%>%
  summarise(expect=sum(mueff) -2*mean(meanrgrNoTreat))
                                                                                                                                                            

ggplot(plotnorm , aes(y= rgr-meanrgrNoTreat, x=Treatment, col=Treatment, fill=Treatment, group=Treatment ) )+
  geom_hline(data=expect,aes(yintercept=expect) , linetype="dashed")+
  geom_point(size=3)+
  geom_boxplot(aes(group=interaction(CellLine_typeLab,Treatment)),alpha=0.4)+
  #geom_errorbar(aes(group=interaction(CellLine_typeLab,Treatment)),alpha=0.4)+
  #stat_boxplot(geom= 'errorbar', width = 0.6)+
  theme_classic(base_size=26)+
  theme(aspect.ratio = 1, legend.position = "bottom")+
  facet_wrap(~CellLine_typeLab ,scale="free", nrow=1 )+ # scale_y_continuous(breaks=log10(1000*c(0.5,1,2,4,8,16,32,64,128)), labels=1000*c(0.5,1,2,4,8,16,32,64,128))+
  #scale_color_manual(name="FGF2 \n(ng/ml)",values=cls[c(1,3)])+scale_fill_manual(name="FGF2 \n(ng/ml)",values=cls[c(1,3)])+
  scale_color_discrete(name="",labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  scale_fill_discrete(name="" ,labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  scale_x_discrete(labels=rep("",4))+
  labs(y="Cancer growth rate \n (compared to DMSO)",x="Treatment")+
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black") )
#ggsave(file=paste0(SvLoc,"Fig1/Fig1D growth rate.pdf"), height=12, width=18, dpi=320)
#ggsave(file=paste0(SvLoc,"Fig1/Fig1D growth rate.png"), height=12, width=18, dpi=320)



ggplot(plotnorm , aes(y= rgr-meanrgrNoTreat, x=Treatment, col=Treatment, fill=Treatment, group=Treatment ) )+
  geom_hline(data=expect,aes(yintercept=expect) , linetype="dashed")+
  geom_point(size=3)+
  geom_boxplot(aes(group=interaction(CellLine_typeLab,Treatment)),alpha=0.4)+
  #geom_errorbar(aes(group=interaction(CellLine_typeLab,Treatment)),alpha=0.4)+
  #stat_boxplot(geom= 'errorbar', width = 0.6)+
  theme_classic(base_size=26)+
  theme(aspect.ratio = 1, legend.position = "bottom")+
  facet_wrap(~CellLine_typeLab ,scale="free", nrow=1 )+ # scale_y_continuous(breaks=log10(1000*c(0.5,1,2,4,8,16,32,64,128)), labels=1000*c(0.5,1,2,4,8,16,32,64,128))+
  #scale_color_manual(name="FGF2 \n(ng/ml)",values=cls[c(1,3)])+scale_fill_manual(name="FGF2 \n(ng/ml)",values=cls[c(1,3)])+
  scale_color_discrete(name="",labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  scale_fill_discrete(name="" ,labels=c("DMSO", "Ribociclib","Afatinib", "Combination"))+
  scale_x_discrete(labels=rep("",4))+
  labs(y="Cancer growth rate \n (compared to DMSO)",x="Treatment")+
  theme(text = element_blank(),legend.text = element_blank(),legend.title = element_blank()  )
#ggsave(file=paste0(SvLoc,"Fig1/BLANK_Fig1D growth rate.pdf"), height=12, width=18, dpi=320)
#ggsave(file=paste0(SvLoc,"Fig1/BLANK_Fig1D growth rate.png"), height=12, width=18, dpi=320)


# Statistics summarizing growth rate differences
# Linear mixed effects model partitioning the 
#m1<-lmer(rgr-meanrgrNoTreat ~ Treatment+ (1+Treatment| CellLine_typeLab), data=plotnorm)
m1<-lmer(rgr-meanrgrNoTreat ~ -1+isRiboTreated +isAfatTreated+ isCombTreated + (0+isRiboTreated +isAfatTreated+ isCombTreated| CellLine_typeLab), data=plotnorm,REML=F)
m2<-lmer(rgr-meanrgrNoTreat ~ -1+isRiboTreated +isAfatTreated+ (0+isRiboTreated +isAfatTreated+ isCombTreated| CellLine_typeLab), data=plotnorm,REML=F)

#model comparison (with vs without synergy) using LRT and AIC
anova(m1,m2,test="chisq")
AIC(m1,m2)

# extract summary statistics
summary(m1)
plot( predict(m1 ), (plotnorm$rgr -plotnorm$meanrgrNoTreat) )

#cell line specific linear model
summary(lm(rgr-meanrgrNoTreat ~ Treatment, data=plotnorm[CellLine_typeLab=="CAMA-1 Resistant"]))
summary(lm(rgr-meanrgrNoTreat ~ Treatment, data=plotnorm[CellLine_typeLab=="CAMA-1 Sensitive"]))
summary(lm(rgr-meanrgrNoTreat ~ Treatment, data=plotnorm[CellLine_typeLab=="MCF-7 Resistant"]))
summary(lm(rgr-meanrgrNoTreat ~ Treatment, data=plotnorm[CellLine_typeLab=="MCF-7 Sensitive"]))

#growth suppression synergy::CAMA-1 Resistant:Est=-0.10,sd=0.0028,t=-37.95,p=7.2e-14; CAMA-1 Sensitive:Est=-0.21,sd=0.0058,t=-36.45,p=1.17e-13;MCF-7 Resistant:Est=-0.12,sd=0.0037,t=-31.47,p=6.7e-13; MCF-7 Sensitive:Est=-0.083,sd=0.0024,t=-34.38,p=2.34e-13

