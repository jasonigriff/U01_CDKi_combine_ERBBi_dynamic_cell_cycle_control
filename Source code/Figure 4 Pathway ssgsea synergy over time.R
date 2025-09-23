rm(list=ls())
require(msigdbr);require(data.table); require(dplyr); require(ggplot2);require(ggsci);require(colorspace)
require(R.utils)
require(umap)
require(Rdimtools)
require(openxlsx)
require(lme4)
require(lmerTest)
require(parallel)
library( readxl )
library(metap);require(ggsci)

# load timecourse rnaseq data analysis results from synergy gam model for Resistant and Sensitive cell lines
fileloc <- ("~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/")
setwd(fileloc)
nsheets<-length( excel_sheets( "sh10050_MCF7_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx" ) )
dd0mcfS <- rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx",sheet=sheet)
}))
dd0camaS <- rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx",sheet=sheet)
}))

Rloc <-paste0(fileloc,"reverseComparison/")
dd0mcfR <- fread(paste0(Rloc,"sh10050_MCF7_H.C2.C5.C6_gam_short.term_parametric_table_synergy.txt"))
dd0mcfR[grepl( "Sen", p.covariate  ),State:="Sens"]
dd0camaR <- fread(paste0(Rloc,"sh10050_CAMA1_H.C2.C5.C6_gam_short.term_parametric_table_synergy.txt"))
dd0camaR[grepl( "Sen", p.covariate  ),State:="Sens"]

# common pathways across datasets
intrcols<-intersect(intersect(names(dd0camaR),names(dd0mcfR)), intersect(names(dd0camaS),names(dd0mcfS)) )
intrGS<-intersect(intersect(unique(dd0camaR$Gene_Set),unique(dd0mcfR$Gene_Set) ), intersect(unique(dd0camaS$Gene_Set),unique(dd0mcfS$Gene_Set)  ) )

# join data
dd0ResS <- data.table(rbind(data.table(CellLine="CAMA1",ReferenceState="Sensitive",dd0camaS%>%dplyr::select(intrcols )),
                            data.table(CellLine="MCF7", ReferenceState="Sensitive",dd0mcfS%>%dplyr::select( intrcols))) )[!grepl("_DN",Gene_Set)][ 
                              grepl("BIOCARTA_",Gene_Set)|grepl("HALLMARK_",Gene_Set) ][Gene_Set%in%intrGS][Timescale=="LongTerm"][State=="Sens"] #|grepl("KEGG_",Gene_Set)|grepl("REACTOME_",Gene_Set)

dd0ResR <- data.table(rbind(data.table(CellLine="CAMA1",ReferenceState="Resistant",dd0camaR%>%dplyr::select(intrcols )),
                            data.table(CellLine="MCF7", ReferenceState="Resistant",dd0mcfR%>%dplyr::select( intrcols ))) )[!grepl("_DN",Gene_Set)][ 
                              grepl("BIOCARTA_",Gene_Set)|grepl("HALLMARK_",Gene_Set) ][Gene_Set%in%intrGS][Timescale=="LongTerm"][State=="RiboR"] #|grepl("KEGG_",Gene_Set)|grepl("REACTOME_",Gene_Set)

dd0Res <- rbind(dd0ResS, dd0ResR)

# set levels
dd0Res$Treatment <- factor(dd0Res$Treatment,c("DMSO", "AfatMain","RiboMain", "CombSynergy"))
dd0Res$State <- factor(dd0Res$State,c("Sens","RiboR"))
dd0Res$CellLine <- factor(dd0Res$CellLine,c("CAMA1", "MCF7"))   #dd0Res[, Day:=(Hour-as.numeric(as.character(timeunderdrug)))/24 ]
dd0Res$Gene_SetLab <- gsub("_"," ",dd0Res$Gene_Set)
dd0Res$ReferenceState <- factor(dd0Res$ReferenceState,c("Sensitive", "Resistant"))   #dd0Res[, Day:=(Hour-as.numeric(as.character(timeunderdrug)))/24 ]

parlist<-dd0Res$p.covariate%>%unique()
selpars<-parlist[grepl("CumDay",parlist)& grepl("CombTreated",parlist)]

# Subset synergey param data
dd0Res<-dd0Res[p.covariate%in%selpars]

# # reactiv estimate and p value (FDR) across timepoints
dd0Res[,syner:= sum( Estimate ) , by=c("CellLine","Gene_Set","State","Treatment") ]
dd0Res[,FDRsyner:= sumlog(FDR)$p , by=c("CellLine","Gene_Set","State","Treatment") ]
dd0Res[,consdir_synerDay:=  2*abs(mean(Estimate>0)) -1 , by= c("Gene_Set","Treatment","Day") ]
dd0Res[,mean_FDRsynerDay:= sumlog(FDR)$p , by= c("Gene_Set","Treatment","Day") ]
dd0Res[,mean_EstsynerDay:= mean(Estimate) , by= c("Gene_Set","Treatment","Day") ]

 #dd0Res[State=="Sens"][Treatment=="RiboMain"][grep("HALLMARK_E2",Gene_Set)]
# 
# # compare cell lines within resistance state (mean reactivation, reactivation direction consistency (-1 ->1)  )
dd0Res[,mean_syner:= mean(syner), by=c("Gene_Set","Treatment","State") ]
dd0Res[,consdir_syner:=  2*abs(mean(syner>0)) -1 , by= c("Gene_Set","Treatment","State") ]
dd0Res[,signif_syner:= FDRsyner<0.05 ]
dd0Res[,mean_FDRsyner:= sumlog(unique(FDRsyner))$p , by= c("Gene_Set","Treatment","State","Day") ]
dd0Res[,consist_signif_syner:= mean(   signif_syner * consdir_syner) , by=c("Gene_Set","Treatment","State") ]
 
# # compare across all cell lines
dd0Res[,Overallmean_syner:= mean(syner), by=c("Gene_Set","Treatment") ]
dd0Res[,Overallconsdir_syner:=  2*abs(mean(syner>0)) -1 , by= c("Gene_Set","Treatment") ]
dd0Res[,Overallsignif_syner:= FDRsyner<0.05 ]
dd0Res[,Overallmean_FDRsyner:= sumlog(unique(FDRsyner))$p , by= c("Gene_Set","Treatment","Day") ]
dd0Res[,Overallconsist_signif_syner:= mean(   Overallsignif_syner * Overallconsdir_syner) , by=c("Gene_Set","Treatment") ]

# vars used to filter: FDR; mean_FDRsynerDay; consdir_synerDay;mean_EstsynerDay
# requires:Estimate
data.table(dd0Res[Treatment!="DMSO"][Day==6][FDR<0.05][abs(consdir_synerDay)>0.5]%>%group_by(Gene_SetLab)%>%summarise(n=n()))[n>1]
GSplot6<-data.table(dd0Res[Treatment!="DMSO"][Day==6][FDR<0.05][mean_FDRsynerDay<0.05][abs(consdir_synerDay)==1][abs(mean_EstsynerDay)>.5]%>%group_by(Gene_SetLab)%>%summarise(mean_EstsynerDay=mean(mean_EstsynerDay),n=n()))[n>1][order(-mean_EstsynerDay)]$Gene_SetLab%>%unique()
GSplot12<-data.table(dd0Res[Treatment!="DMSO"][Day==12][FDR<0.05][mean_FDRsynerDay<0.05][abs(consdir_synerDay)==1][abs(mean_EstsynerDay)>.5]%>%group_by(Gene_SetLab)%>%summarise(mean_EstsynerDay=mean(mean_EstsynerDay),n=n()))[n>1][order(-mean_EstsynerDay)]$Gene_SetLab%>%unique()
GSplot18<-data.table(dd0Res[Treatment!="DMSO"][Day==18][FDR<0.05][mean_FDRsynerDay<0.05][abs(consdir_synerDay)==1][abs(mean_EstsynerDay)>.5]%>%group_by(Gene_SetLab)%>%summarise(mean_EstsynerDay=mean(mean_EstsynerDay),n=n()))[n>1][order(-mean_EstsynerDay)]$Gene_SetLab%>%unique()
# At each timepoint we identified pathways showing consistent synergistic combination treatment effects across cell lines (Fisher’s combined FDR<0.05; mean(effect)>0.5; consistent sign(effect);  FDR<0.05 in multiple cell lines
#GSplot6<-dd0Res[Treatment!="DMSO"][Day==6][FDR<0.05][mean_FDRsynerDay<0.05][abs(consdir_synerDay)>0.75][order(mean_EstsynerDay)]$Gene_SetLab%>%unique()#%>%sort()
#GSplot12<-dd0Res[Treatment!="DMSO"][Day==12][FDR<0.05][mean_FDRsynerDay<0.05][abs(consdir_synerDay)>0.75][order(mean_EstsynerDay)]$Gene_SetLab%>%unique()#%>%sort()
#GSplot18<-dd0Res[Treatment!="DMSO"][Day==18][FDR<0.05][mean_FDRsynerDay<0.05][abs(consdir_synerDay)>0.75][order(mean_EstsynerDay)]$Gene_SetLab%>%unique()#%>%sort()

GSplot<-unique(c(GSplot6,GSplot12,GSplot18))

tt<- dd0Res[Day==6][Gene_SetLab%in%GSplot6][order(-mean_EstsynerDay)]
tt$Gene_SetLab<- factor(tt$Gene_SetLab, levels=GSplot6)
ggplot(tt%>%group_by(Gene_SetLab,CellLine,State,Day)%>%dplyr::summarise(ActEst=sum(Estimate)), aes(y=ActEst, x=Gene_SetLab, col=State,
                                                                                                          group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(~Day)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=6)+
  theme(aspect.ratio=2.5)+
  labs(y="Synergy", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))

tt<- dd0Res[Day==12][Gene_SetLab%in%GSplot12]
tt$Gene_SetLab<- factor(tt$Gene_SetLab, levels=GSplot12)
ggplot(tt%>%group_by(Gene_SetLab,CellLine,State,Day)%>%dplyr::summarise(ActEst=sum(Estimate)), aes(y=ActEst, x=Gene_SetLab, col=State,
                                                                                                   group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(~Day)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=6)+
  theme(aspect.ratio=2.5)+
  labs(y="Synergy", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))


tt<- dd0Res[Day==18][Gene_SetLab%in%GSplot18]
tt$Gene_SetLab<- factor(tt$Gene_SetLab, levels=GSplot18)
ggplot(tt%>%group_by(Gene_SetLab,CellLine,State,Day)%>%dplyr::summarise(ActEst=sum(Estimate)), aes(y=ActEst, x=Gene_SetLab, col=State,
                                                                                                                                       group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(~Day)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=6)+
  theme(aspect.ratio=2.5)+
  labs(y="Synergy", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))



#Bad :GSplot <- dd0Res[Treatment!="DMSO"][abs(consist_signif_syner)>=0.5][mean_FDRsyner<0.05][(abs(Overallconsist_signif_syner)>0.25) & (abs(Overallmean_syner)>0.)]$Gene_SetLab%>%unique()%>%sort()
# GSplot
twotypes <- dd0Res [Gene_SetLab %in% GSplot] 
twotypes[Gene_Set=="HALLMARK_E2F_TARGETS"][p.covariate=="CumDay6:CombTreated"]
twotypes[Gene_Set=="HALLMARK_MITOTIC_SPINDLE"]
dd0Res[Treatment!="DMSO"][Gene_Set=="BIOCARTA_MITOCHONDRIA_PATHWAY"][CellLine=="CAMA1"]

#twotypes[,SensEst:=sum(Estimate*(State=="Sens")), by=c("CellLine","Gene_Set","Day")]
#twotypes[,ActEst:=SensEst+ Estimate*(State!="Sens")]
twotypes[,ActEst:=sum(Estimate), by=c("CellLine","Gene_Set","Day","State")]

twotypes$DayLab <- paste0("Day:",twotypes$Day)
twotypes[DayLab=="Day:6",DayLab:="Day 0-6"]
twotypes[DayLab=="Day:12",DayLab:="Day 6-12"]
twotypes[DayLab=="Day:18",DayLab:="Day 12-18"]
twotypes$DayLab <- factor(twotypes$DayLab , levels=unique(twotypes$DayLab ) )


twotypes$Gene_SetLab <- factor(twotypes$Gene_SetLab , levels=rev(unique( twotypes[order(Overallmean_syner)]$Gene_SetLab )))
#twotypes$Gene_SetLab <- factor(twotypes$Gene_SetLab , levels=rev(GSplot))


ggplot(twotypes%>%group_by(Gene_SetLab,CellLine,State)%>%dplyr::summarise(ActEst=sum(Estimate)), aes(y=ActEst, x=Gene_SetLab, col=State,
                                                                                                          group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  #facet_grid(~DayLab)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=5)+coord_flip()+
  theme_classic(base_size=26)+
  theme(aspect.ratio=2.5)+
  labs(y="Synergistic effect on pathway", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))+ 
  theme(text = element_text(colour = "black"),
        axis.text = element_text(colour = "black" ))

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Synergy effect R vs Sv3.png", width=16, height=16, dpi=320)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Synergy effect R vs Sv3.pdf", width=16, height=16, dpi=320)


ggplot(twotypes%>%
         group_by(Gene_SetLab,CellLine,State,DayLab)%>%
         dplyr::summarise(ActEst=sum(Estimate)), aes(y=ActEst, x=Gene_SetLab, col=State,
                                                                                                   group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(~DayLab)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=5)+coord_flip()+
  theme_classic(base_size=26)+
  theme(aspect.ratio=2.5)+
  labs(y="Synergistic effect on pathway", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))+
  theme(text = element_text(colour = "black"),
        axis.text = element_text(colour = "black" ))


#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Synergy effect accumulation over time R vs Sv3.png", width=20, height=20, dpi=320)
#write.csv( twotypes, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/ SItable I Fig4a Combination synergy effects over time in R and S CAMA1 and MCF7.csv")

subsetGSplt<-c("HALLMARK E2F TARGETS","HALLMARK G2M CHECKPOINT","HALLMARK MITOTIC SPINDLE","HALLMARK MYC V1","HALLMARK MYC V2",#"BIOCARTA MPR PATHWAY","BIOCARTA EFP PATHWAY",
               "BIOCARTA MCM PATHWAY","BIOCARTA G1 PATHWAY","BIOCARTA CELLCYCLE PATHWAY","BIOCARTA G2 PATHWAY","BIOCARTA RACCYCD PATHWAY",
               "HALLMARK APOPTOSIS","HALLMARK P53 PATHWAY","BIOCARTA DEATH PATHWAY","BIOCARTA MITOCHONDRIA PATHWAY")
#write.csv( twotypes, file= )
           
ggplot(twotypes[Gene_SetLab%in%subsetGSplt]%>%
         group_by(Gene_SetLab,CellLine,State,DayLab)%>%
         dplyr::summarise(ActEst=sum(Estimate)), aes(y=ActEst, x=Gene_SetLab, col=State,
                                                   group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(~DayLab)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=5)+coord_flip()+
  theme_classic(base_size=6)+
  theme(aspect.ratio=1.5)+
  #geom_vline(xintercept=4.5, col="grey")+
  labs(y="Synergistic effect on pathway", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig4a Synergy effect accumulation over time R vs Sv2 trimmed.png"), width=16, height=16, dpi=600)
#ggsave(file=paste0(MSfigureoutput,"Fig4a Synergy effect accumulation over time R vs Sv2 trimmed.pdf"), width=16, height=16, dpi=600)

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Synergy effect accumulation over time R vs Sv2 trimmed.png", width=16, height=16, dpi=320)


## conclusion= Cell cycle, MYC and DNA damage synergistically turned off at t=6 and further at t=12 but no further by t=18

subsetGS<-c("HALLMARK E2F TARGETS","HALLMARK G2M CHECKPOINT",#"HALLMARK MITOTIC SPINDLE", 
            "BIOCARTA DEATH PATHWAY","BIOCARTA MITOCHONDRIA PATHWAY")
dynamSyn<- data.table(twotypes[Gene_SetLab%in%subsetGS]%>%
  group_by(Gene_SetLab,CellLine,State,DayLab,Day)%>%
  dplyr::summarise(ActEst=sum(ActEst)))
dynamSyn$Gene_SetLab <- factor(dynamSyn$Gene_SetLab, levels=subsetGS)
ggplot(dynamSyn, aes(y=ActEst, x=Day))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(~Gene_SetLab)+
  geom_smooth(method="gam",formula=y~s(x, k=3))+
  geom_point(aes(shape=CellLine, col=State,
                 group=interaction(Day,CellLine)),position=position_dodge(width=0.5),size=3)+
  theme_classic(base_size=6)+
  theme(aspect.ratio=1)+
  labs(y="Synergy", x="Day")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))

ggplot(dynamSyn, aes(y=ActEst, x=DayLab,group=Gene_SetLab, fill=Gene_SetLab,col=Gene_SetLab))+
  geom_hline(yintercept = 0, linetype="dashed")+
  #facet_grid(~Gene_SetLab)+
 # geom_smooth(method="gam",formula=y~s(x, k=3))+
  geom_boxplot(col="black",aes(group=interaction(DayLab,Gene_SetLab)))+
  geom_point(pch=21,col="black",aes(fill=Gene_SetLab,
                 group=interaction(DayLab,Gene_SetLab)),position=position_dodge(width=0.75),size=3)+
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+
  labs(y="Synergistic effect", x="")+
  scale_color_manual(name="Pathway", values=c("blue","skyblue","red","darkred"),labels=c("Hallmark E2F Targets", "Hallmark G2M Checkpoint","Biocarta Death Pathway", "Biocarta Mitochondria Pathway")) +
  scale_fill_manual(name="Pathway", values=c("blue","skyblue","red","darkred"),labels=c("Hallmark E2F Targets", "Hallmark G2M Checkpoint","Biocarta Death Pathway", "Biocarta Mitochondria Pathway"))+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Synergy effect vs day four pathwaysB.png", width=11, height=11, dpi=320)


MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig4b Synergy effect vs day four pathways.png"), width=11, height=11, dpi=600)
#ggsave(file=paste0(MSfigureoutput,"Fig4b Synergy effect vs day four pathways.pdf"), width=11, height=11, dpi=600)

summary( lm(ActEst~-1+DayLab,   data= dynamSyn[Gene_SetLab=="BIOCARTA DEATH PATHWAY"]) )
summary( lm(ActEst~-1+DayLab,   data= dynamSyn[Gene_SetLab=="BIOCARTA MITOCHONDRIA PATHWAY"]) )
summary( lm(ActEst~-1+DayLab,   data= dynamSyn[Gene_SetLab=="HALLMARK E2F TARGETS"]) )
summary( lm(ActEst~-1+DayLab,   data= dynamSyn[Gene_SetLab=="HALLMARK G2M CHECKPOINT"]) )


# Biocarta death pathway:: 
# Day=0-6:Est=-0.10,sd=0.21, t=-0.46,p=0.66,   
# Day=6-12:Est=1.07,sd=0.21,t=5.01,p=0.0007,  
# Day=12-18:Est=0.92,sd=0.21,t=4.30,p=0.002 
# 
# Biocarta mitochondria pathway:: 
# Day=0-6:Est=0.11,sd=0.17, t=0.64,p=0.54,   
# Day=6-12:Est=0.38,sd=0.17,t=2.17,p=0.058,  
# Day=12-18:Est=1.05,sd=0.17,t=6.03,p=0.0002 
# 
# 
# Hallmark E2F targets:: 
# Day=0-6:Est=-0.59,sd=0.22, t=-2.72,p=0.02,   
# Day=6-12:Est=-0.89,sd=0.22,t=-4.06,p=0.002,  
# Day=12-18:Est=-0.22,sd=0.22,t=-0.99,p=0.35 
# 
# Biocarta G2M checkpoint:: 
# Day=0-6:Est=-0.61,sd=0.22, t=-2,80,p=0.02,   
# Day=6-12:Est=-0.87,sd=0.22,t=-3.99,p=0.003,  
# Day=12-18:Est=-0.23,sd=0.22,t=-1.06,p=0.32 
# 
# 















unique(dd0Res[grepl("HALLMARK",Gene_Set)|grepl("BIOCARTA",Gene_Set)][p.covariate%in%selpars ][FDR <0.05]$Gene_Set)


plotThis <- dd0Res[grepl("HALLMARK",Gene_Set)|grepl("BIOCARTA",Gene_Set)][p.covariate%in%selpars ][FDR <0.05][][order(-Estimate)]
plotThis$Gene_SetLab <- gsub("_"," ",plotThis$Gene_Set)


orderplot<-data.table(plotThis %>% group_by(Gene_SetLab,p.covariate)%>%dplyr::summarise(mu=mean(Estimate),consistdir=-1+sum(Estimate>0), n=n() ))[order(mu)]
orderplot2<-data.table(dd0Res[grepl("HALLMARK",Gene_Set)|grepl("BIOCARTA",Gene_Set)][p.covariate%in%selpars]  %>% group_by(Gene_SetLab)%>%dplyr::summarise(mu=mean(Estimate),consistdir=-1+sum(Estimate>0), n=n() ))[order(mu)]

plotThis2 <- plotThis[ Gene_SetLab%in% orderplot[n>1][consistdir!=0]
                       [#(!grepl("StateRiboR",p.covariate))
                         #&grepl("CumDay6",p.covariate)
                       ]$Gene_SetLab][
                         #  (!grepl("StateRiboR",p.covariate))&grepl("CumDay6",p.covariate)
                       ]
unique(plotThis2$Gene_SetLab)
#ggplot(plotThis2,aes(x=Estimate, y=FDR))+geom_point() +facet_wrap(State~Day)# %>% hist()

plotThis2[State=="Sens"&Day==6]
plotThis2[State=="Sens"&Day==12]
plotThis2[State=="Sens"&Day==18]

plotThis2[Day==6]
plotThis2[Day==12]
plotThis2[Day==18]
plotThis2$Gene_Set%>%unique()

plotThis2[]$Gene_Set%>%unique()
# plotThis2$Gene_Set%>%unique()
#  Cell cycle:  "HALLMARK_E2F_TARGETS"           "HALLMARK_G2M_CHECKPOINT"       "HALLMARK_MITOTIC_SPINDLE"      "BIOCARTA_CELLCYCLE_PATHWAY"      
#                 "BIOCARTA_RACCYCD_PATHWAY"         "BIOCARTA_G2_PATHWAY"      "BIOCARTA_MCM_PATHWAY"
#  MYC:         "HALLMARK_MYC_TARGETS_V1"        ,      "HALLMARK_MYC_TARGETS_V2" 
#  DNA damage: "BIOCARTA_ARF_PATHWAY"                 
#  Death: "BIOCARTA_MITOCHONDRIA_PATHWAY"             "BIOCARTA_DEATH_PATHWAY"    "HALLMARK_INFLAMMATORY_RESPONSE"
#  Other:  "BIOCARTA_MPR_PATHWAY"    "BIOCARTA_EFP_PATHWAY"  "BIOCARTA_NKCELLS_PATHWAY"      "HALLMARK_HYPOXIA" 


ggplot(dd0Res[p.covariate%in%selpars ][Gene_Set%in% unique(plotThis2[]$Gene_Set) ] [], aes(y=Estimate, x=Gene_SetLab, fill=Estimate,
                                                                                           group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(.~Day)+
  geom_point(aes(shape=State),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2.5)+ labs(y="Resistant cell pathway activition \n (pre-tx compared to sensitive cells)", x="Gene set")+
  scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF7"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))

# 
# ggplot(dd0Res[p.covariate%in%selpars ][Gene_Set%in% unique(plotThis2[]$Gene_Set) ] [State=="Sens"], aes(y=Estimate, x=Gene_SetLab, fill=Estimate,
#                                                                                                         group=interaction(Gene_SetLab,CellLine)))+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   facet_grid(State~Day)+
#   geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
#   theme_classic(base_size=12)+
#   theme(aspect.ratio=2.5)+ labs(y="Resistant cell pathway activition \n (pre-tx compared to sensitive cells)", x="Gene set")+
#   scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF7"))+#,values=cols[-1])+
#   scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))
# 
# ggplot(dd0Res[p.covariate%in%selpars ][Gene_Set%in% unique(plotThis2[State!="Sens"]$Gene_Set) ] [State=="Sens"], aes(y=Estimate, x=Gene_SetLab, fill=Estimate,
#                                                                                                                      group=interaction(Gene_SetLab,CellLine)))+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   facet_grid(State~Day)+
#   geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
#   theme_classic(base_size=12)+
#   theme(aspect.ratio=2.5)+ labs(y="Resistant cell pathway activition \n (pre-tx compared to sensitive cells)", x="Gene set")+
#   scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF7"))+#,values=cols[-1])+
#   scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))

#twotypes<-dd0Res[p.covariate%in%selpars ][Gene_Set%in% c( unique(plotThis2[State!="Sens"]$Gene_Set) , "BIOCARTA_ATRBRCA_PATHWAY") ]
twotypes <- dd0Res[p.covariate %in% selpars ][Gene_Set %in% c( unique(plotThis2[]$Gene_Set)) ]#, "BIOCARTA_ATRBRCA_PATHWAY") ]
twotypes[Gene_Set=="HALLMARK_E2F_TARGETS"][p.covariate=="CumDay6:CombTreated"]
twotypes[Gene_Set=="HALLMARK_E2F_TARGETS"][CellLine=="CAMA1"][Day==6]

#twotypes[,SensEst:=sum(Estimate*(State=="Sens")), by=c("CellLine","Gene_Set","Day")]
#twotypes[,ActEst:=SensEst+ Estimate*(State!="Sens")]
twotypes[,ActEst:=sum(Estimate), by=c("CellLine","Gene_Set","Day","State")]

twotypes$DayLab <- paste0("Day:",twotypes$Day)
twotypes[DayLab=="Day:6",DayLab:="Day 0-6"]
twotypes[DayLab=="Day:12",DayLab:="Day 6-12"]
twotypes[DayLab=="Day:18",DayLab:="Day 12-18"]
twotypes$DayLab <- factor(twotypes$DayLab , levels=unique(twotypes$DayLab ) )




# 
# twotypes$Gene_SetLab <- factor(twotypes$Gene_SetLab , levels=rev( c("HALLMARK E2F TARGETS","HALLMARK G2M CHECKPOINT" ,"HALLMARK MITOTIC SPINDLE","BIOCARTA CELLCYCLE PATHWAY" ,
#                                                                     "BIOCARTA G2 PATHWAY","BIOCARTA MCM PATHWAY",  
#                                                                     "BIOCARTA DEATH PATHWAY" ,"BIOCARTA MITOCHONDRIA PATHWAY" ,
#                                                                     "HALLMARK MYC TARGETS V1","HALLMARK MYC TARGETS V2",
#                                                                     "BIOCARTA ARF PATHWAY","BIOCARTA ATRBRCA PATHWAY",
#                                                                     "BIOCARTA EFP PATHWAY", "BIOCARTA MPR PATHWAY","HALLMARK HYPOXIA", "HALLMARK INFLAMMATORY RESPONSE" ,"BIOCARTA NKCELLS PATHWAY" ) ))
# 
# twotypes$Gene_SetLab <- factor(twotypes$Gene_SetLab , levels=rev( c("BIOCARTA G2 PATHWAY", "BIOCARTA CELLCYCLE PATHWAY" ,"HALLMARK MITOTIC SPINDLE",
#                                                                     "HALLMARK E2F TARGETS","HALLMARK G2M CHECKPOINT" ,
#                                                                     "BIOCARTA MCM PATHWAY",   "BIOCARTA RACCYCD PATHWAY",
#                                                                     "BIOCARTA DEATH PATHWAY" ,"BIOCARTA MITOCHONDRIA PATHWAY" ,
#                                                                     "HALLMARK MYC TARGETS V1","HALLMARK MYC TARGETS V2",
#                                                                     "BIOCARTA ARF PATHWAY","BIOCARTA ATRBRCA PATHWAY",
#                                                                     "BIOCARTA EFP PATHWAY", "BIOCARTA MPR PATHWAY","HALLMARK HYPOXIA" , "HALLMARK INFLAMMATORY RESPONSE" ,"BIOCARTA NKCELLS PATHWAY" ) ))
# 
# twotypes$Gene_SetLab <- factor(twotypes$Gene_SetLab , levels=rev( c("BIOCARTA G2 PATHWAY", "BIOCARTA CELLCYCLE PATHWAY" ,"HALLMARK MITOTIC SPINDLE",
#                                                                     "HALLMARK E2F TARGETS","HALLMARK G2M CHECKPOINT" ,
#                                                                     "BIOCARTA MCM PATHWAY",   "BIOCARTA RACCYCD PATHWAY",
#                                                                     
#                                                                     "HALLMARK MYC TARGETS V1","HALLMARK MYC TARGETS V2",
#                                                                     "BIOCARTA ARF PATHWAY","BIOCARTA ATRBRCA PATHWAY",
#                                                                     "BIOCARTA EFP PATHWAY", "BIOCARTA MPR PATHWAY",
#                                                                     "BIOCARTA DEATH PATHWAY" ,"BIOCARTA MITOCHONDRIA PATHWAY" ,
#                                                                     "HALLMARK INFLAMMATORY RESPONSE" ,"HALLMARK HYPOXIA" , "BIOCARTA NKCELLS PATHWAY" ) ))


ggplot(twotypes, aes(y=ActEst, x=Gene_SetLab, col=State,
                     group=interaction(Gene_SetLab)))+
  geom_boxplot(fill="grey",outlier.color =NA, show.legend = F)+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(~DayLab)+
  #geom_point(aes(group=interaction(Gene_SetLab,CellLine), shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  #geom_jitter(aes(group=interaction(Gene_SetLab,CellLine), shape=CellLine),width=0.5,size=3)+coord_flip()+
  geom_point(aes(group=interaction(Gene_SetLab,State,CellLine), shape=CellLine),position=position_dodge(width=0.25),size=3)+coord_flip()+
  theme_classic(base_size=10)+
  theme(aspect.ratio=1.5)+
  labs(y="Cumulative synergistic pathway modulation during combination treatment \n (additional to individual treatment effects)", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7")) 
#subtwotypes<- twotypes[!Gene_SetLab%in%c("HALLMARK HYPOXIA","")]


#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Synergy effect accumulation over time R vs S.png", width=16, height=16, dpi=320)





twotypes$Gene_SetLab <- factor(twotypes$Gene_SetLab , levels=rev(unique(orderplot2[n>1][Gene_SetLab%in%unique(twotypes$Gene_SetLab)]$Gene_SetLab) ))

ggplot(twotypes, aes(y=ActEst, x=Gene_SetLab, col=State,
                     group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(~DayLab)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2.5)+
  labs(y="Cumulative synergistic pathway modulation \n (compared to individual treatment effects)", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))

## conclusion= Cell cycle, MYC and DNA damage synergistically turned off at t=6 and further at t=12 but no further by t=18
ggplot(twotypes[abs(Overallconsdir_syner)>0.5][abs(Overallconsist_signif_syner)>0.05]%>%group_by(Gene_SetLab,CellLine,State)%>%dplyr::summarise(ActEst=sum(ActEst)), aes(y=ActEst, x=Gene_SetLab, col=State,
                                                                                                   group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  #facet_grid(~DayLab)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=16)+
  theme(aspect.ratio=2.5)+
  labs(y="Synergy", x="Gene set")+
  scale_color_manual(name="Resistance", values=c("blue","red"),labels=c("Sensitive", "Resistant"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Synergy effect overall over time R vs S.png", width=16, height=16, dpi=320)


ggplot(plotThis2[] , aes(y=Estimate, x=Gene_SetLab, fill=Estimate,
                         group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  facet_grid(State~Day)+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2.5)+ labs(y="Resistant cell pathway activition \n (pre-tx compared to sensitive cells)", x="Gene set")+
  scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF7"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))

ggplot(plotThis2[] , aes(y=Day, x=Gene_SetLab, fill=Estimate,
                         group=interaction(Gene_SetLab,CellLine)))+
  geom_tile()+
  #geom_hline(yintercept = 0, linetype="dashed")+
  facet_wrap(~interaction(CellLine,State), nrow=1)+
  #geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+
  coord_flip()+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2.5)#+# labs(y="Resistant cell pathway activition \n (pre-tx compared to sensitive cells)", x="Gene set")+
#  scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF7"))+#,values=cols[-1])+
# scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))




plotThis$Gene_SetLab<- factor(plotThis$Gene_SetLab, levels=unique(orderplot$Gene_Set))
nplot<-10
wchplot<-orderplot[c(1:nplot,((nrow(orderplot)-(nplot-1)):nrow(orderplot)))]$Gene_SetLab
subplotThis <- plotThis[Gene_SetLab%in%wchplot][Gene_SetLab%in%orderplot[consistdir!=0]$Gene_SetLab  ]
subplotThis$Gene_SetLab<- factor(subplotThis$Gene_SetLab, levels=unique(wchplot))


subplotThis2<-merge(subplotThis,orderplot, by="Gene_SetLab")
ggplot(subplotThis , aes(y=Estimate, x=Gene_SetLab, col=CellLine,
                         group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2.5)+ labs(y="Resistant cell pathway activition \n (pre-tx compared to sensitive cells)", x="Gene set")+
  scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF7"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))







rm(list=ls())
require(msigdbr);require(data.table); require(dplyr); require(ggplot2);require(ggsci);require(colorspace)
require(R.utils)
require(umap)
require(Rdimtools)
require(openxlsx)
require(lme4)
require(lmerTest)
require(parallel)
library( readxl )

# load timecourse rnaseq data
fileloc <- ("~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/")
setwd(fileloc)
nsheets<-length( excel_sheets( "sh10050_MCF7_H.C2.C5.C6_gam_short.term_smooths_table_synergy.xlsx" ) )

dd0mcf <-rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_smooths_table_synergy.xlsx",sheet=sheet)
}))

#read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_smooths_table_synergy.xlsx")
dd0cama<-rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_smooths_table_synergy.xlsx",sheet=sheet)
}))


dd0Res <- data.table(rbind(data.table(CellLine="CAMA1",dd0cama%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )),
                           data.table(CellLine="MCF7", dd0mcf%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ))) )[!grepl("_DN",Gene_Set)][ 
                             grepl("BIOCARTA_",Gene_Set)|grepl("HALLMARK_",Gene_Set)]# |grepl("KEGG_",Gene_Set)|grepl("REACTOME_",Gene_Set)]

dd0Res$Treatment <- factor(dd0Res$Treatment,c("DMSO", "AfatMain","RiboMain", "CombSynergy"))
dd0Res$State <- factor(dd0Res$State,c("Sens","RiboR"))
dd0Res$CellLine <- factor(dd0Res$CellLine,c("CAMA1", "MCF7"))
#dd0Res[, Day:=(Hour-as.numeric(as.character(timeunderdrug)))/24 ]

dd0Res[Treatment=="CombSynergy"][FDR <0.05]

# Ribo treatment leads to cell cycle reactivation !
# Cell cycle and MYC reactivation
dd0Res[CellLine=="CAMA1"][State=="Sens"][Treatment=="RiboMain"][Hour==24][p.value <0.05][order(-(effect))]
dd0Res[CellLine=="MCF7"][State=="Sens"][Treatment=="RiboMain"][Hour==24][p.value <0.05][order(-(effect))]
# not reactivated under ribo in riboR as never turned off in first place
dd0Res[CellLine=="CAMA1"][State=="RiboR"][Treatment=="RiboMain"][Hour==24][p.value <0.05][order(-(effect))]
dd0Res[CellLine=="MCF7"][State=="RiboR"][Treatment=="RiboMain"][Hour==24][p.value <0.05][order(-(effect))]

dd0Res[CellLine=="CAMA1"][State=="Sens"][Treatment=="AfatMain"][Hour==24][p.value <0.05][order(-(effect))]
dd0Res[CellLine=="MCF7"][State=="Sens"][Treatment=="AfatMain"][Hour==24][p.value <0.05][order(-(effect))]
dd0Res[CellLine=="CAMA1"][State=="RiboR"][Treatment=="AfatMain"][Hour==24][p.value <0.05][order(-(effect))]
dd0Res[CellLine=="MCF7"][State=="RiboR"][Treatment=="AfatMain"][Hour==24][p.value <0.05][order(-(effect))]






commonIncreaseR <- intersect( 
  dd0Res[CellLine=="CAMA1"][State=="Sens"][Treatment=="RiboMain"][Hour==24][FDR <0.05][effect>0]$Gene_Set,
  dd0Res[CellLine=="MCF7"][State=="Sens"][Treatment=="RiboMain"][Hour==24][FDR <0.05][effect>0]$Gene_Set
)
commonDecreaseR <- intersect( 
  dd0Res[CellLine=="CAMA1"][State=="Sens"][Treatment=="RiboMain"][Hour==24][FDR <0.05][effect<0]$Gene_Set,
  dd0Res[CellLine=="MCF7"][State=="Sens"][Treatment=="RiboMain"][Hour==24][FDR <0.05][effect<0]$Gene_Set
)


commonIncreaseA<-intersect( 
  dd0Res[CellLine=="CAMA1"][State=="Sens"][Treatment=="AfatMain"][Hour==24][FDR <0.05][effect>0]$Gene_Set,
  dd0Res[CellLine=="MCF7"][State=="Sens"][Treatment=="AfatMain"][Hour==24][FDR <0.05][effect>0]$Gene_Set
)
commonDecreaseA<-intersect( 
  dd0Res[CellLine=="CAMA1"][State=="Sens"][Treatment=="AfatMain"][Hour==24][FDR <0.05][effect<0]$Gene_Set,
  dd0Res[CellLine=="MCF7"][State=="Sens"][Treatment=="AfatMain"][Hour==24][FDR <0.05][effect<0]$Gene_Set
)


commonIncreaseC<-intersect( 
  dd0Res[CellLine=="CAMA1"][State=="Sens"][Treatment=="CombSynergy"][Hour==24][FDR <0.05][effect>0]$Gene_Set,
  dd0Res[CellLine=="MCF7"][State=="Sens"][Treatment=="CombSynergy"][Hour==24][FDR <0.05][effect>0]$Gene_Set
)
commonDecreaseC<-intersect( 
  dd0Res[CellLine=="CAMA1"][State=="Sens"][Treatment=="CombSynergy"][Hour==24][FDR <0.05][effect<0]$Gene_Set,
  dd0Res[CellLine=="MCF7"][State=="Sens"][Treatment=="CombSynergy"][Hour==24][FDR <0.05][effect<0]$Gene_Set
)


lu2<-rbind( data.table(Treatment="RiboMain",Direction= "Increase", Gene_Set=commonIncreaseR),
            data.table(Treatment="RiboMain",Direction= "Decrease",Gene_Set=commonDecreaseR),
            data.table(Treatment="AfatMain",Direction= "Increase", Gene_Set=commonIncreaseA),
            data.table(Treatment="AfatMain",Direction= "Decrease",Gene_Set=commonDecreaseA),
            data.table(Treatment="CombSynergy",Direction= "Increase", Gene_Set=commonIncreaseC),
            data.table(Treatment="CombSynergy",Direction= "Decrease",Gene_Set=commonDecreaseC)
)
lu2finalized<- data.table(lu2)
lu2finalized$Treatment<- gsub("Main","",lu2finalized$Treatment)
saveloc<- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/Lab U01/SynergyPlasticitySSGSEArewiring/"
#write.csv(lu2finalized, file=paste0(saveloc,"H.C2.C5.C6 plasticity rewiring by treatment.csv" ))

Res1finalized <-  dd0Res[Treatment!="DMSO"][Gene_Set%in%lu2finalized$Gene_Set][Hour==24]#[p.value <0.05][State=="Sens"]
#write.csv(Res1finalized, file=paste0(saveloc,"H.C2.C5.C6Results plasticity rewiring by treatment.csv" ))
Res1finalized[p.value <0.05][State=="Sens"][order(Treatment,-effect)]

plotThis<- data.table( Res1finalized)
plotThis$Gene_SetLab <- gsub("_"," ",plotThis$Gene_Set)
orderplot<-data.table(plotThis[State=="Sens"][Treatment=="RiboMain"]%>%group_by(Gene_Set)%>%summarise(mu=mean(effect)))[order(mu)]
plotThis$Gene_SetLab<- factor(plotThis$Gene_SetLab, levels=gsub("_"," ",orderplot$Gene_Set))

#gsub("_"," ",Gene_Set)

plotThis[,StateLab:="Sensitive"];plotThis[State=="RiboR",StateLab:="Resistant"]
plotThis$StateLab<- factor(plotThis$StateLab, levels=c("Sensitive","Resistant"))
plotThis[,TreatmentLab:="Combination"];plotThis[Treatment=="AfatMain",TreatmentLab:="Afatinib"];plotThis[Treatment=="RiboMain",TreatmentLab:="Ribociclib"]
plotThis$TreatmentLab<- factor(plotThis$TreatmentLab, levels=c("Afatinib","Ribociclib","Combination"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(4)

ggplot(plotThis[] , aes(y=effect, x=Gene_SetLab, col=TreatmentLab, group=interaction(Gene_SetLab,CellLine)))+geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=26)+
  facet_wrap(~StateLab) + theme(aspect.ratio=2)+ labs(y="Rewiring effect", x="Gene set")+
  scale_color_manual(name="Treatment",values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Rewiring effect of R vs S under drug.png", width=16, height=16, dpi=320)

ggplot(plotThis[grep("HALLMARK",Gene_Set)] , aes(y=effect, x=Gene_SetLab, col=TreatmentLab, group=interaction(Gene_SetLab,CellLine)))+geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=26)+
  facet_wrap(~StateLab) + theme(aspect.ratio=1.8)+ labs(y="Rewiring effect", x="Gene set")+
  scale_color_manual(name="Treatment",values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Rewiring effect of R vs S under drug HALLMARK.png", width=15, height=14, dpi=320)


# load timecourse rnaseq data
list.files(fileloc)
#dd0mcf <-read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx")
#dd0cama<-read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx")
dd0mcf <-rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx",sheet=sheet)
}))

#read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_smooths_table_synergy.xlsx")
dd0cama<-rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx",sheet=sheet)
}))

dd0 <- data.table(rbind(dd0cama%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ),dd0mcf%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )))

dd0[,Lineage:=paste0(CellLine," ",State)]
dd0$Lineage <- factor(dd0$Lineage  , levels=c("CAMA1 Sen" ,"MCF7 Sen", "CAMA1 RiboR"  ,  "MCF7 RiboR"  ))

dd0lu <- dd0[Gene_Set%in%lu2$Gene_Set]

dd0lu[,scaleSSGSEA:= scale(Enrichment_Score), by=c("Gene_Set","CellLine","State")]
dd0[,scaleSSGSEAB:= scale(Enrichment_Score), by=c("Gene_Set","CellLine")]

dd0lu[,SensExpression.SSGSEA := sum(Enrichment_Score*(State=="Sen"&Treatment=="DMSO")), by=c("Gene_Set","CellLine","Hour")]
dd0lu[,SensExpression.SSGSEAB:= sum(Enrichment_Score*(Treatment=="DMSO" )), by=c("Gene_Set","CellLine","Hour","State")]
dd0lu[,NormSensExpression.SSGSEA:=Enrichment_Score -SensExpression.SSGSEA, by=c("Gene_Set","CellLine","Hour", "Treatment")]
dd0lu[,NormSensExpression.SSGSEAB:=Enrichment_Score -SensExpression.SSGSEAB, by=c("Gene_Set","CellLine","Hour", "Treatment")]
dd0lu$Treatment <- factor(dd0lu$Treatment,c("DMSO", "Afat","Ribo", "Comb"))
dd0lu$Hour%>%unique()
dd0lu$Gene_Set%>%unique()
p1<-ggplot(dd0lu[Treatment!="DMSO"],
           aes(x=as.factor(Hour), y=Gene_Set, fill=NormSensExpression.SSGSEA ))+
  theme_classic(base_size=16)+
  geom_tile()+
  facet_grid(Treatment~Lineage,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_fill_gradient2(name="Expression \n log2FC \n vs \n sensitive cells \n under DMSO" , low="blue", high="darkred", mid="white")+
  theme(aspect.ratio=3)+#labs(title=gsub("_"," ",x), x="Hour")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p1

p1<- ggplot(dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1" )],
  aes(x=as.factor(Hour), y=gsub("_"," ",Gene_Set), fill=NormSensExpression.SSGSEA ))+
  theme_classic(base_size=16)+
  geom_tile()+
  facet_grid(Treatment~Lineage,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_fill_gradient2(name="ssGSEA pathway activity \n vs \n sensitive cells \n under DMSO" , low="blue", high="darkred", mid="white")+
  theme(aspect.ratio=3)+labs(y="Pathway", x="Hour")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p1


require(ggsci)
p2<- ggplot(dd0lu[Treatment!="DMSO"],# [
            #  Gene_Set%in%c( "HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2" )],
            aes(x=Hour, y=NormSensExpression.SSGSEA, col= gsub("_"," ", gsub("HALLMARK_","",Gene_Set)), group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=26)+
  geom_point(size=2.5)+ geom_line()+
  facet_grid(Treatment~Lineage)+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_color_npg(name="Pathway", )+
  theme(aspect.ratio=3)#+#labs(title=gsub("_"," ",x), x="Hour")#+
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2

egdatTx0 <- dd0lu[timeunderdrug==0] [ Gene_Set%in%c( "HALLMARK_E2F_TARGETS" )][Lineage%in% c("CAMA1 RiboR","CAMA1 Sen" ) ]
egdatTx0$Day<-as.numeric(as.character(egdatTx0$day)); egdatTx0[Day==1, Day:=0]
ggplot(egdatTx0[Lineage%in% c("CAMA1 RiboR")],
       aes(x=Hour/24, y=NormSensExpression.SSGSEAB, fill=Treatment,col= Treatment, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0)/24, linetype=2,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=6.5, aes(alpha=1/(1+(Day))))+ geom_line()+
  #facet_grid(.~State ,scales="free")+
  #scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  #scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title="Combined tx effects", y="Pathway modulation
                             \n (relative to DMSO)", x="Day")+
  scale_color_discrete(name="Treatment",labels=c("DMSO","Afatinib","Ribociclib","Combination") )+
  scale_fill_discrete(name="Treatment",labels=c("DMSO","Afatinib","Ribociclib","Combination") )+
  scale_x_continuous(breaks=c(0,6,12,18))+
  scale_alpha( range = c(0.4, 1),guide = "none")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic tx effects.png", width=8, height=8, dpi=320)


egdatTx <- dd0lu[] [ Gene_Set%in%c( "HALLMARK_E2F_TARGETS" )][Lineage%in% c("CAMA1 RiboR","CAMA1 Sen" ) ]
egdatTx$Day<-as.numeric(as.character(egdatTx$day)); egdatTx[Day==1, Day:=0]
ggplot(egdatTx[Lineage%in% c("CAMA1 RiboR")],
       aes(x=Hour/24, y=NormSensExpression.SSGSEAB, fill=Treatment,col= Treatment, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0)/24, linetype=2,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=6.5, aes(alpha=1/(1+(Day))))+ geom_line()+
  #facet_grid(.~State ,scales="free")+
  #scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  #scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title="Combined tx effects", y="Pathway modulation
                             \n (relative to DMSO)", x="Day")+
  scale_color_discrete(name="Treatment",labels=c("DMSO","Afatinib","Ribociclib","Combination") )+
  scale_fill_discrete(name="Treatment",labels=c("DMSO","Afatinib","Ribociclib","Combination") )+
  scale_x_continuous(breaks=c(0,6,12,18))+
  scale_alpha( range = c(0.4, 1),guide = "none")
ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic tx effects2.png", width=8, height=8, dpi=320)



egdat <- dd0lu[Treatment=="Ribo"][Lineage%in% c("CAMA1 Sen","CAMA1 RiboR" ) ] [ Gene_Set%in%c( "HALLMARK_E2F_TARGETS" )]
egdat$Day<-as.numeric(as.character(egdat$day)); egdat[Day==1, Day:=0]
#egdat$Day<- factor(egdat$Day, levels=c("0","6","12","18"))
ggplot(egdat,
       aes(x=Hour/24, y=NormSensExpression.SSGSEA, fill=State,col= State, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0, 6,12,18)/24, linetype=2,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=6.5, aes(alpha=1/(1+(Day))))+ geom_line()+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title="HALLMARK E2F TARGETS", y="Pathway activity", x="Day")+
  scale_x_continuous(breaks=c(0,6,12,18))+
  scale_alpha( range = c(0.4, 1),guide = "none")

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic input data.png", width=8, height=8, dpi=320)

ggplot(egdat[Hour==0],
       aes(x=State, y=NormSensExpression.SSGSEA, fill=State,col= State, group=interaction(State) ))+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=6.5)+ 
  # geom_segment(data=egdat[Hour<=24][Hour>0],aes(x = 0, 
  #                                               y = Hr0score, 
  #                                               xend = Hour, 
  #                                               yend = NormSensExpression.SSGSEA
  # ),
  # arrow = arrow(length = unit(0.5, "cm") , type="closed") ,show.legend = F)+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1, legend.position = "none")+labs(title="Baseline ", y="Pathway activity", x="")+
  scale_x_discrete(labels=c("Resistant", "Sensitive"))

egdat[,Hr0score:=sum(NormSensExpression.SSGSEA*(Hour==0)), by="State"]
ggplot(egdat[Hour<=24],
       aes(x=Hour, y=NormSensExpression.SSGSEA, fill=State,col= State, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0)/24, linetype=2,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=6.5, aes(alpha=1/(1+(Day))))+ 
  geom_segment(data=egdat[Hour<=24][Hour>0],aes(x = 0, 
                                                y = Hr0score, 
                                                xend = Hour, 
                                                yend = NormSensExpression.SSGSEA
  ),
  arrow = arrow(length = unit(0.5, "cm") , type="closed") ,show.legend = F)+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title="HALLMARK E2F TARGETS", y="Pathway activity", x="Hour")+
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  scale_alpha( range = c(0.4, 1),guide = "none")

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic input data P1.png", width=8, height=8, dpi=320)

ggplot(egdat[Hour<=24],
       aes(x=Hour, y=NormSensExpression.SSGSEA, fill=State,col= State, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0)/24, linetype=2,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=11, aes(alpha=1/(1+(Day))))+ 
  geom_segment(data=egdat[Hour<=24][Hour>0],aes(x = 0, 
                                                y = Hr0score, 
                                                xend = Hour, 
                                                yend = NormSensExpression.SSGSEA
  ),
  arrow = arrow(length = unit(0.75, "cm") , type="closed") ,show.legend = F)+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title= "    Baseline difference & \n Initial tx response (Day 0)", y="Pathway activity", x="Hour")+
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  scale_alpha( range = c(0.4, 1),guide = "none")+theme(legend.position="none")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic input data P1.png", width=7, height=7, dpi=320)

inpt2<- egdat[timeunderdrug==0]
inpt2[ ,lagy:=lag(NormSensExpression.SSGSEA), by=State]
inpt2[ ,lagx:=lag(Hour), by=State]
inpt2[ ,lagX:=lag(Day), by=State]

ggplot(inpt2,
       aes(x=Day, y=NormSensExpression.SSGSEA, fill=State,col= State, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0,6,12,18)/24, linetype=2,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=11, aes(alpha=1/(1+(Day))))+ 
  geom_segment(data=inpt2,aes(x = lagX, 
                              y = lagy, 
                              xend = Day, 
                              yend = NormSensExpression.SSGSEA
  ),
  arrow = arrow(length = unit(0.75, "cm") , type="closed") ,show.legend = F)+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title="Long-term trend \n     (Day 0-18)", y="Pathway activity", x="Day")+
  scale_x_continuous(breaks=c(0,6,12,18))+
  scale_alpha( range = c(0.4, 1),guide = "none")+theme(legend.position="none")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic input data P2.png", width=7, height=7, dpi=320)



egdat3 <- dd0lu[Treatment=="Ribo"][Lineage%in% c("CAMA1 Sen","CAMA1 RiboR" ) ] [ Gene_Set%in%c( "HALLMARK_E2F_TARGETS" )]
egdat3[,HrDrug0score:=sum(NormSensExpression.SSGSEA*(timeunderdrug==0)), by=c("State","day")]
egdat3$Day<-as.numeric(as.character(egdat3$day)); egdat3[Day==1, Day:=0]


ggplot(egdat3[],
       aes(x=as.numeric(timeunderdrug), y=NormSensExpression.SSGSEA, fill=State,col= State, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0)/24, linetype=2,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=11, aes(alpha=1/(1+(Day))))+ 
  geom_segment(data=egdat3[][as.numeric(timeunderdrug)>0],aes(x = 0, 
                                                              y = HrDrug0score, 
                                                              xend = as.numeric(timeunderdrug), 
                                                              yend = NormSensExpression.SSGSEA
  ),
  arrow = arrow(length = unit(0.75, "cm") , type="closed") ,show.legend = F)+
  facet_grid(.~Day,scales="free")+
  scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title="Initial tx response (Day 0)", y="Pathway activity", x="Hour")+
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  scale_alpha( range = c(0.4, 1),guide = "none")

ggplot(egdat3[timeunderdrug==24],
       aes(x=Day+1, y=(24)*(NormSensExpression.SSGSEA-HrDrug0score)/as.numeric(timeunderdrug), fill=State,col= State, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0,6,12,18)/24, linetype=2,col="grey")+
  geom_hline(yintercept=0, linetype=3,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=11, aes(alpha=1/(1+(Day))))+ 
  geom_segment(data=egdat3[timeunderdrug==24][as.numeric(timeunderdrug)>0],aes(x = Day, 
                                                                               y = 0, 
                                                                               xend = Day+1, 
                                                                               yend = (24)*(NormSensExpression.SSGSEA-HrDrug0score)/as.numeric(timeunderdrug)
  ),
  arrow = arrow(length = unit(0.75, "cm") , type="closed") ,show.legend = F)+
  #facet_grid(.~timeunderdrug,scales="free")+
  scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title="Plasticity of tx responses \n           (Day 0-18)", y="Tx response (activation/day)", x="Day")+
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  scale_alpha( range = c(0.4, 1),guide = "none")+theme(legend.position="none")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic input data P3.png", width=7, height=7, dpi=320)


dd0lu[,Treatmentlab:="Combination"]
dd0lu[Treatment=="Afat",Treatmentlab:="Afatinib"]
dd0lu[Treatment=="Ribo",Treatmentlab:="Ribociclib"]
dd0lu[Treatment=="DMSO",Treatmentlab:="DMSO"]
dd0lu$Treatmentlab <- factor(dd0lu$Treatmentlab,  levels=c("DMSO" ,"Afatinib", "Ribociclib","Combination"  ))

p2<- ggplot(dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1" )],
  aes(x=Hour, y=NormSensExpression.SSGSEA, col= gsub("_"," ", gsub("HALLMARK_","",Gene_Set)), group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=26)+
  geom_hline(linetype="dashed", aes(yintercept=0))+
  geom_point(size=2.5)+ geom_line()+
  # facet_grid(Treatmentlab~Lineage)+
  facet_grid(Lineage~Treatmentlab)+
  scale_color_npg(name="Pathway", )+
  theme(aspect.ratio=1)+#+#labs(title=gsub("_"," ",x), x="Hour")+
  labs(y= "ssGSEA pathway activity \n (vs  sensitive cell under DMSO)")
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2
sloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/GeneExpressionTimeCourse/Plasticity/"
#ggsave(p2,file=paste0(sloc,file="PlasticitySSGSEA by treatment_.pdf"), dpi=320,height=24,width=18)


p2<- ggplot(dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("BIOCARTA_G1_PATHWAY","BIOCARTA_G2_PATHWAY","BIOCARTA_MCM_PATHWAY", "BIOCARTA_RACCYCD_PATHWAY" ,"BIOCARTA_CELLCYCLE_PATHWAY")],
  aes(x=Hour, y=NormSensExpression.SSGSEA, col= gsub("_"," ", gsub("HALLMARK_","",Gene_Set)), group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=26)+
  geom_hline(linetype="dashed", aes(yintercept=0))+
  geom_point(size=2.5)+ geom_line()+
  # facet_grid(Treatmentlab~Lineage)+
  facet_grid(Lineage~Treatmentlab)+
  scale_color_npg(name="Pathway", )+
  theme(aspect.ratio=1)+#+#labs(title=gsub("_"," ",x), x="Hour")+
  labs(y= "ssGSEA pathway activity \n (vs  sensitive cell under DMSO)")
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2
sloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/GeneExpressionTimeCourse/Plasticity/"
#ggsave(p2,file=paste0(sloc,file="PlasticitySSGSEA by treatment_.pdf"), dpi=320,height=24,width=18)

p2<- ggplot(dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("BIOCARTA_MPR_PATHWAY","BIOCARTA_CARM_ER_PATHWAY","BIOCARTA_EFP_PATHWAY")],
  aes(x=Hour, y=NormSensExpression.SSGSEA, col= gsub("_"," ", gsub("HALLMARK_","",Gene_Set)), group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=26)+
  geom_hline(linetype="dashed", aes(yintercept=0))+
  geom_point(size=2.5)+ geom_line()+
  # facet_grid(Treatmentlab~Lineage)+
  facet_grid(Lineage~Treatmentlab)+
  scale_color_npg(name="Pathway", )+
  theme(aspect.ratio=1)+#+#labs(title=gsub("_"," ",x), x="Hour")+
  labs(y= "ssGSEA pathway activity \n (vs  sensitive cell under DMSO)")
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2
sloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/GeneExpressionTimeCourse/Plasticity/"
#ggsave(p2,file=paste0(sloc,file="PlasticitySSGSEA by treatment_.pdf"), dpi=320,height=24,width=18)

p2<- ggplot(dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("BIOCARTA_ARF_PATHWAY","BIOCARTA_ATRBRCA_PATHWAY", "BIOCARTA_ATM_PATHWAY" ,"BIOCARTA_P53_PATHWAY")],
  aes(x=Hour, y=NormSensExpression.SSGSEA, col= gsub("_"," ", gsub("HALLMARK_","",Gene_Set)), group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=26)+
  geom_hline(linetype="dashed", aes(yintercept=0))+
  geom_point(size=2.5)+ geom_line()+
  # facet_grid(Treatmentlab~Lineage)+
  facet_grid(Lineage~Treatmentlab)+
  scale_color_npg(name="Pathway", )+
  theme(aspect.ratio=1)+#+#labs(title=gsub("_"," ",x), x="Hour")+
  labs(y= "ssGSEA pathway activity \n (vs  sensitive cell under DMSO)")
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2
sloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/GeneExpressionTimeCourse/Plasticity/"
#ggsave(p2,file=paste0(sloc,file="PlasticitySSGSEA by treatment_.pdf"), dpi=320,height=24,width=18)



p2<- ggplot(dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","BIOCARTA_ATRBRCA_PATHWAY", "HALLMARK_MYC_TARGETS_V1" )],
  aes(x=Hour, y=NormSensExpression.SSGSEA, col= Treatmentlab, group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=26)+
  geom_hline(linetype="dashed", aes(yintercept=0))+
  geom_point(size=2.5)+ geom_line()+
  facet_grid(gsub("_"," ", gsub("HALLMARK_","",Gene_Set))~Lineage)+
  scale_color_npg(name="Treatment", )+
  theme(aspect.ratio=1)+#+#labs(title=gsub("_"," ",x), x="Hour")+
  labs(y= "ssGSEA pathway activity \n (vs  sensitive cell under DMSO)")
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2



PlotTHIS<- dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","BIOCARTA_ATRBRCA_PATHWAY", "HALLMARK_MYC_TARGETS_V1" )]
PlotTHIS$Gene_SetLab<-  gsub("_"," ",PlotTHIS$Gene_Set)

PlotTHIS[Gene_SetLab=="HALLMARK E2F TARGETS",Gene_SetLab:="Hallmark \n E2F Targets"]
PlotTHIS[Gene_SetLab=="HALLMARK G2M CHECKPOINT",Gene_SetLab:="Hallmark \n G2M Checkpoint"]
PlotTHIS[Gene_SetLab=="HALLMARK MYC TARGETS V1",Gene_SetLab:="Hallmark \n MYC Targets V1"]
PlotTHIS[Gene_SetLab=="BIOCARTA ATRBRCA PATHWAY",Gene_SetLab:="Biocarta \n ATRBRCA pathway"]
PlotTHIS$Gene_SetLab<- factor(PlotTHIS$Gene_SetLab, labels=c("Hallmark \n E2F Targets","Hallmark \n G2M Checkpoint","Hallmark \n MYC Targets V1" ,"Biocarta \n ATRBRCA pathway"  ))

PlotTHIS[   ,LineageLab:="Sensitive \n CAMA-1"]
PlotTHIS[Lineage=="MCF7 Sen",LineageLab:="Sensitive \n MCF7"]
PlotTHIS[Lineage=="CAMA1 RiboR",LineageLab:="Resistant \n CAMA-1"]
PlotTHIS[Lineage=="MCF7 RiboR",LineageLab:="Resistant \n MCF7"]
PlotTHIS$LineageLab<- factor(PlotTHIS$LineageLab, labels=c("Sensitive \n CAMA-1","Sensitive \n MCF7","Resistant \n CAMA-1" ,"Resistant \n MCF7"  ))

p2<- ggplot(PlotTHIS,
            aes(x=Hour, y=NormSensExpression.SSGSEA, col= Treatmentlab, group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=26)+
  geom_hline(linetype="dashed", aes(yintercept=0))+
  geom_point(size=2.5)+ geom_line()+
  facet_grid( Gene_SetLab~LineageLab)+
  scale_color_npg(name="Treatment", )+
  theme(aspect.ratio=1)+#+#labs(title=gsub("_"," ",x), x="Hour")+
  labs(y= "ssGSEA pathway activity \n (vs  sensitive cell under DMSO)")
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2
sloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/GeneExpressionTimeCourse/Plasticity/"
#ggsave(p2,file=paste0(sloc,file="PlasticitySSGSEA by treatment_v3.pdf"), dpi=320,height=24,width=18)



p2<- ggplot(dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2","HALLMARK_MTORC1_SIGNALING" )],
  aes(x=Hour, y=NormSensExpression.SSGSEA, col= gsub("_"," ", gsub("HALLMARK_","",Gene_Set)), group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=16)+
  geom_point()+ geom_line()+
  facet_grid(Treatment~Lineage,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_color_npg(name="Pathway", )+
  theme(aspect.ratio=3)#+#labs(title=gsub("_"," ",x), x="Hour")#+
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2

