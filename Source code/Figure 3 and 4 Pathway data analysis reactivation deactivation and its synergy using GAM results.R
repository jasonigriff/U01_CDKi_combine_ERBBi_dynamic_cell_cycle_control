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

# Specify source data location 
fileloc <- ("~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/")
setwd(fileloc)

## Load parametric inferences from GAM model of time course RNAseq data for CAMA-1 and MCF-7
# Primary comparison (sensitive cells as reference)
nsheets <- length( excel_sheets( "sh10050_MCF7_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx" ) )
dd0mcfS <- rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx",sheet=sheet)
}))
dd0camaS <- rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx",sheet=sheet)
}))

# Reverse comparison (resistant cells as reference)
Rloc <-paste0(fileloc,"reverseComparison/")
dd0mcfR <- fread(paste0(Rloc,"sh10050_MCF7_H.C2.C5.C6_gam_short.term_parametric_table_synergy.txt"))
dd0mcfR[grepl( "Sen", p.covariate  ),State:="Sens"]
dd0camaR <- fread(paste0(Rloc,"sh10050_CAMA1_H.C2.C5.C6_gam_short.term_parametric_table_synergy.txt"))
dd0camaR[grepl( "Sen", p.covariate  ),State:="Sens"]

# Identify common pathways across data sets
intrcols <- intersect(intersect(names(dd0camaR),names(dd0mcfR)), intersect(names(dd0camaS),names(dd0mcfS)) )
intrGS <- intersect(intersect(unique(dd0camaR$Gene_Set),unique(dd0mcfR$Gene_Set) ), intersect(unique(dd0camaS$Gene_Set),unique(dd0mcfS$Gene_Set)  ) )

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

# get scores per treatment from additive + interaction effects #dd0Res[State=="Sens"] [grep("HALLMARK_E2",Gene_Set)][CellLine=="CAMA1"][Day==18]
dd0Res[,EstimateAbs:= sum(Estimate*(Treatment!="DMSO")) , by=c("CellLine","State","Day","Gene_Set")]
dd0Res[Treatment!="CombSynergy",EstimateAbs:= Estimate ]

# Define the early vs late phase
dd0Res[,Phase:="Early"]
dd0Res[Day>=12,Phase:="Late"]

# reactiv estimate and p value (FDR) across timepoints
dd0Res[,reactiv:= sum( (-(Day==min(Day))*Estimate)+((Day!=min(Day))*Estimate) ) , by=c("CellLine","Gene_Set","State","Treatment") ]
dd0Res[,FDRreactiv:= sumlog(FDR)$p , by=c("CellLine","State","Gene_Set","Treatment") ]

# compare cell lines within resistance state (mean reactivation, reactivation direction consistency (-1 ->1)  )
dd0Res[,mean_reactiv:= mean(reactiv), by=c("Gene_Set","Treatment","State") ]
dd0Res[,consdir_reactiv:=  2*abs(mean(reactiv>0)) -1 , by= c("Gene_Set","Treatment","State") ]
dd0Res[,signif_reactiv:= FDRreactiv<0.05 ]
dd0Res[,mean_FDRreactiv:= sumlog(unique(FDRreactiv))$p , by= c("Gene_Set","Treatment","State","Day") ]
dd0Res[,consist_signif_reactiv:= mean(   signif_reactiv * consdir_reactiv) , by=c("Gene_Set","Treatment","State") ]

# compare across all cell lines
dd0Res[,Overallmean_reactiv:= mean(reactiv), by=c("Gene_Set","Treatment") ]
dd0Res[,Overallconsdir_reactiv:=  2*abs(mean(reactiv>0)) -1 , by= c("Gene_Set","Treatment") ]
dd0Res[,Overallsignif_reactiv:= FDRreactiv<0.05 ]
dd0Res[,Overallmean_FDRreactiv:= sumlog(unique(FDRreactiv))$p , by= c("Gene_Set","Treatment","Day") ]
dd0Res[,Overallconsist_signif_reactiv:= mean(   Overallsignif_reactiv * Overallconsdir_reactiv) , by=c("Gene_Set","Treatment") ]

#abs(consist_signif_reactiv)==1
#mean_FDRreactiv<0.05
#(abs(Overallconsist_signif_reactiv)>0) & (abs(Overallmean_reactiv)>0.5)

#dd0Res[Treatment!="DMSO"][Gene_SetLab%in%"HALLMARK TNFA SIGNALING VIA NFKB"][State=="Sens"][Treatment=="RiboMain"]
# view results
dd0Res[abs(consist_signif_reactiv)==1][mean_FDRreactiv<0.05][Treatment!="DMSO"][Phase=="Early"][Treatment=="RiboMain"][order(-mean_reactiv)]
dd0Res[abs(consist_signif_reactiv)==1][mean_FDRreactiv<0.05][Treatment!="DMSO"][Phase=="Early"][Treatment=="AfatMain"][order(-mean_reactiv)]
dd0Res[abs(consist_signif_reactiv)==1][mean_FDRreactiv<0.05][Treatment!="DMSO"][Phase=="Early"][Treatment=="CombSynergy"][order(-mean_reactiv)]

# Select consistently and significantly rewiring pathways in either resistant or sensitive cells
GSplot<- dd0Res[Treatment!="DMSO"][Phase=="Early"][abs(consist_signif_reactiv)==1][mean_FDRreactiv<0.05][(abs(Overallconsist_signif_reactiv)>0) & (abs(Overallmean_reactiv)>0.5)]$Gene_SetLab%>%unique()%>%sort()

# Plotting results
plotthis <- dd0Res[Phase=="Early"][Treatment!="DMSO"][Gene_SetLab%in%GSplot]
ordgs2 <- rev(unique(rev(dd0Res[Gene_SetLab%in%GSplot][CellLine=="CAMA1"][][Phase=="Early"][State=="Sens"][Treatment=="RiboMain"|Treatment=="AfatMain"] [order(-Overallmean_reactiv)]$Gene_SetLab)))

# order pathways and treatments and label cell lines
plotthis$Gene_SetLab<- factor(plotthis$Gene_SetLab , levels=rev(ordgs2))
plotthis$Treatment<- factor(plotthis$Treatment , levels=c("RiboMain","AfatMain","CombSynergy"))
plotthis[,CellLineLab:="CAMA-1"]
plotthis[CellLine=="MCF7",CellLineLab:="MCF-7"]


# extended results plot
topN <- 25
whcplot<-ordgs2[c(1:topN, (length(ordgs2)-topN+1): length(ordgs2))]
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(4)
pl1dd<-plotthis#[Gene_SetLab%in%whcplot] 
ggplot( pl1dd , aes(y=(Gene_SetLab),x=reactiv,col=Treatment,shape=CellLineLab, group=interaction(Gene_SetLab,CellLine)))+
  theme_classic(base_size=26)+theme(aspect.ratio=2)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_point(position=position_dodge(width=0.5),size=3)+facet_grid(.~ReferenceState)+
  labs(y="Gene set", x="Reactivation/deactivation \n (-ve=turning off, +ve=turning on)") +
  scale_color_manual(name="Treatment \n effect",labels= c("Ribociclib","Afatinib","Synergy"),values=cols[-1])+# , labels=c("Ribociclib","Afatinib","Combination"))
  scale_fill_manual(name="Treatment \n effect",labels= c("Ribociclib","Afatinib","Synergy"),values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/FigS8a Reactiv or deactiv drug effect on R and S.png", width=22.5, height=18.5, dpi=320)

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Long term rewiring/Fig4 Rewiring effect of R and S under drug.png", width=22, height=18, dpi=320)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Long term rewiring/Fig4 Rewiring effect of R and S under drug.pdf", width=22, height=18, dpi=320)

#write.csv( pl1dd, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/ SourceData I Fig4 Rewiring effect of R and S under drug.csv")

# Save output table
dataoutput <- unique( plotthis%>%select(CellLine,ReferenceState,Treatment,Gene_Set,
                                     # Estimate, Std.Error, t_value, p.value, 
                                      reactiv,FDRreactiv,mean_reactiv,mean_FDRreactiv,Overallmean_reactiv,Overallmean_FDRreactiv))[order(Treatment,Overallmean_reactiv,CellLine,ReferenceState)]
setnames( dataoutput , 
          old=c("ReferenceState", "reactiv", "FDRreactiv", "mean_reactiv", "mean_FDRreactiv" ,"Overallmean_reactiv","Overallmean_FDRreactiv"),
          new=c("State", "Rewiring effect", "FDR Rewiring effect", 
                "State mean Rewiring effect", "FDR State mean Rewiring effect",
                "Overall mean Rewiring effect", "Overall State mean Rewiring effect"))

#write.csv( dataoutput, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Table all significant output for Rewiring.csv")


# Specific plot of recurrently influential pathways across analyses
sspaths<-c("BIOCARTA_CELLCYCLE_PATHWAY"  ,    
           "BIOCARTA_MCM_PATHWAY"      ,                  
           "HALLMARK_MITOTIC_SPINDLE"  ,      
           "HALLMARK_E2F_TARGETS"     ,       
           "HALLMARK_G2M_CHECKPOINT"  ,     
           "HALLMARK_MYC_TARGETS_V1"   ,   
           "HALLMARK_MYC_TARGETS_V2" ,     
           "BIOCARTA_ATRBRCA_PATHWAY" ,
           "BIOCARTA_DEATH_PATHWAY" ,
           "BIOCARTA_MITOCHONDRIA_PATHWAY" )
SIfig5a <- dataoutput[Gene_Set%in%sspaths] #%>%select(CellLine:Hour,EstimateAbs, Phase,reactiv,consist_signif_reactiv,mean_FDRreactiv,Overallconsist_signif_reactiv,Overallmean_reactiv )
#write.csv( SIfig5a, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/SI Table 5A reactivation.csv")

ggplot( pl1dd[Gene_Set%in%sspaths] , aes(y=(Gene_SetLab),x=reactiv,col=Treatment,shape=CellLineLab, group=interaction(Gene_SetLab,CellLine)))+
  theme_classic(base_size=26)+theme(aspect.ratio=1.6)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_point(position=position_dodge(width=0.5),size=3)+facet_grid(.~ReferenceState)+
  labs(y="Gene set", x="Reactivation/deactivation \n (-ve=turning off, +ve=turning on)") +
  scale_color_manual(name="Treatment \neffect",labels=c("Ribociclib","Afatinib","Synergy"), values=cols[-1])+# , labels=c("Ribociclib","Afatinib","Combination"))
  scale_fill_manual(name="Treatment \neffect",labels=c("Ribociclib","Afatinib","Synergy"),values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7")) +
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig5a Selected10pathways Reactiv or Deactiv of R and S pathways under drug.png"), width=15, height=15, dpi=320)
#ggsave(file=paste0(MSfigureoutput,"Fig5a Selected10pathways Reactiv or Deactiv of R and S pathways under drug.pdf"), width=15, height=15, dpi=320)
##ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Long term rewiring/Fig4 Selected10pathways Rewiring effect of R and S under drug.png", width=15, height=15, dpi=320)
##ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Long term rewiring/Fig4 Selected10pathways Rewiring effect of R and S under drug.pdf", width=15, height=15, dpi=320)

ggplot( pl1dd[Gene_Set%in%sspaths] , aes(y=(Gene_SetLab),x=reactiv,col=Treatment,shape=CellLineLab, group=interaction(Gene_SetLab,CellLine)))+
  theme_classic(base_size=26)+theme(aspect.ratio=1.6)+
  geom_vline(xintercept=0, linetype="dashed")+
  geom_point(position=position_dodge(width=0.5),size=3)+facet_grid(.~ReferenceState)+
  labs(y="Gene set", x="Reactivation/deactivation \n (-ve=turning off, +ve=turning on)") +
  scale_color_manual(name="Treatment \neffect",labels=c("Ribociclib","Afatinib","Synergy"), values=cols[-1])+# , labels=c("Ribociclib","Afatinib","Combination"))
  scale_fill_manual(name="Treatment \neffect",labels=c("Ribociclib","Afatinib","Synergy"),values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7")) +
  theme( text = element_blank(), legend.text = element_blank(), legend.title = element_blank(), )

#ggsave(file=paste0(MSfigureoutput,"BLANK_Fig5a Selected10pathways Reactiv or Deactiv of R and S pathways under drug.png"), width=10, height=10, dpi=320)
#ggsave(file=paste0(MSfigureoutput,"BLANK_Fig5a Selected10pathways Reactiv or Deactiv of R and S pathways under drug.pdf"), width=10, height=10, dpi=320)


# load timecourse rnaseq data
list.files(fileloc)
#dd0mcf <-rbindlist(lapply(1:nsheets,function(sheet){
#  read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx",sheet=sheet)
#}))
#write.csv(dd0mcf, file="sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
#dd0cama<-rbindlist(lapply(1:nsheets,function(sheet){
#  read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx",sheet=sheet)
#}))
#write.csv(dd0cama, file="sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
#setwd to "~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/"
dd0mcf <- fread(file="sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
dd0cama <- fread(file="sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
dd0 <- data.table(rbind(dd0cama%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ),dd0mcf%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )))

#Annotate and order lineage and add gente set label
dd0[,Lineage:=paste0(CellLine," ",State)]
dd0$Lineage <- factor(dd0$Lineage  , levels=c("CAMA1 Sen" ,"MCF7 Sen", "CAMA1 RiboR"  ,  "MCF7 RiboR"  ))
dd0$Gene_SetLab <- gsub("_"," ",dd0$Gene_Set)

# Look up data for focal gene sets
dd0lu <- dd0[Gene_SetLab%in%GSplot]

# Scale and normalize
dd0lu[,scaleSSGSEA:= scale(Enrichment_Score), by=c("Gene_Set","CellLine","State")]
dd0lu[,scaleSSGSEAB:= scale(Enrichment_Score), by=c("Gene_Set","CellLine")]
dd0lu[,SensExpression.SSGSEA := sum(Enrichment_Score*(State=="Sen"&Treatment=="DMSO")), by=c("Gene_Set","CellLine","Hour")]
dd0lu[,SensExpression.SSGSEAB:= sum(Enrichment_Score*(Treatment=="DMSO" )), by=c("Gene_Set","CellLine","Hour","State")]
dd0lu[,NormSensExpression.SSGSEA:=Enrichment_Score -SensExpression.SSGSEA, by=c("Gene_Set","CellLine","Hour", "Treatment")]
dd0lu[,NormSensExpression.SSGSEAB:=Enrichment_Score -SensExpression.SSGSEAB, by=c("Gene_Set","CellLine","Hour", "Treatment")]

# order treatments and add ordered label column
dd0lu$Treatment <- factor(dd0lu$Treatment,c("DMSO", "Afat","Ribo", "Comb"))
dd0lu[,Treatmentlab:="Combination"]
dd0lu[Treatment=="Afat",Treatmentlab:="Afatinib"]
dd0lu[Treatment=="Ribo",Treatmentlab:="Ribociclib"]
dd0lu[Treatment=="DMSO",Treatmentlab:="DMSO"]
dd0lu$Treatmentlab <- factor(dd0lu$Treatmentlab,  levels=c("DMSO" ,"Afatinib", "Ribociclib","Combination"  ))

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


# Select four representative pathways reflecting cell cycle and apoptosis 
PlotTHIS<- dd0lu[Treatment!="DMSO"][
  Gene_Set%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","BIOCARTA_MITOCHONDRIA_PATHWAY",  "BIOCARTA_DEATH_PATHWAY")] 
# Add gene set annotation and order
PlotTHIS$Gene_SetLab<-  gsub("_"," ",PlotTHIS$Gene_Set)
PlotTHIS[Gene_SetLab=="HALLMARK E2F TARGETS",Gene_SetLab:="Hallmark \n E2F Targets"]
PlotTHIS[Gene_SetLab=="HALLMARK G2M CHECKPOINT",Gene_SetLab:="Hallmark \n G2M Checkpoint"]
PlotTHIS[Gene_SetLab=="HALLMARK MYC TARGETS V1",Gene_SetLab:="Hallmark \n MYC Targets V1"]
PlotTHIS[Gene_SetLab=="BIOCARTA ATRBRCA PATHWAY",Gene_SetLab:="Biocarta \n ATRBRCA pathway"]
PlotTHIS[Gene_SetLab=="BIOCARTA MITOCHONDRIA PATHWAY",Gene_SetLab:="Biocarta \n Mitochondria pathway"]
PlotTHIS[Gene_SetLab=="BIOCARTA DEATH PATHWAY",Gene_SetLab:="Biocarta \n Death pathway"]
PlotTHIS$Gene_SetLab<- factor(PlotTHIS$Gene_SetLab, levels=c("Hallmark \n E2F Targets","Hallmark \n G2M Checkpoint","Biocarta \n Death pathway","Biocarta \n Mitochondria pathway"   ))

#Add lineage and treatment label and order
PlotTHIS[   ,LineageLab:="CAMA-1 \n Sensitive"]
PlotTHIS[Lineage=="MCF7 Sen",LineageLab:="MCF-7 \n Sensitive"]
PlotTHIS[Lineage=="CAMA1 RiboR",LineageLab:="CAMA-1 \n Resistant"]
PlotTHIS[Lineage=="MCF7 RiboR",LineageLab:="MCF-7 \n Resistant"]
PlotTHIS$LineageLab<- factor(PlotTHIS$LineageLab, levels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant" ,"MCF-7 \n Sensitive","MCF-7 \n Resistant"  ))
PlotTHIS$Treatmentlab<- factor(PlotTHIS$Treatmentlab, levels=c("Ribociclib","Afatinib","Combination"  ))

p2<- ggplot(PlotTHIS,
            aes(x=Hour/24, y=NormSensExpression.SSGSEA, col= Treatmentlab, group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=26)+
  geom_hline(linetype="dashed", aes(yintercept=0))+
  geom_point(size=3)+ geom_line(linewidth=1.2)+
  facet_grid( Gene_SetLab~LineageLab, scales="free_y")+
  scale_color_manual(name="Treatment",values=cols[-1])+# , labels=c("Ribociclib","Afatinib","Combination"))
  scale_fill_manual(name="Treatment",values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))+
  theme(aspect.ratio=1)+#+#labs(title=gsub("_"," ",x), x="Hour")+
  labs(y= "ssGSEA pathway activity \n (vs  sensitive cell under DMSO)", x="Day")+
  scale_x_continuous(breaks=c(0,6,12,18) )+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

#theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#scale_y_discrete(guide=guide_axis(n.dodge=1)) 
p2
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/PlasticitySSGSEA by treatment_v4.1.png", width=18.5, height=18.5, dpi=320)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/PlasticitySSGSEA by treatment_v4.1.pdf", width=18.5, height=18.5, dpi=320)

MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig5b Selected4pathway timcourse Reactiv or Deactiv of R and S pathways under drug.png"), width=18.5, height=18.5, dpi=320)
#ggsave(file=paste0(MSfigureoutput,"Fig5b Selected4pathway timcourse Reactiv or Deactiv of R and S pathways under drug.pdf"), width=18.5, height=18.5, dpi=320)

p2 +  theme( text = element_blank(), legend.text = element_blank(), legend.title = element_blank(), )
#ggsave(file=paste0(MSfigureoutput,"BLANK_Fig5b Selected4pathway timcourse Reactiv or Deactiv of R and S pathways under drug.png"), width=15, height=15, dpi=320)
#ggsave(file=paste0(MSfigureoutput,"BLANK_Fig5b Selected4pathway timcourse Reactiv or Deactiv of R and S pathways under drug.pdf"), width=15, height=15, dpi=320)









# Heatmaps of pathway activity over time
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



### Timecourse plots for groups of pathways of interest
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
p2
sloc <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/GeneExpressionTimeCourse/Plasticity/"
#ggsave(p2,file=paste0(sloc,file="PlasticitySSGSEA by treatment_.pdf"), dpi=320,height=24,width=18)





### Example timecourse plots used to build GAM model schematic
egdatTx0 <- dd0lu[timeunderdrug==0] [ Gene_Set%in%c( "HALLMARK_E2F_TARGETS" )][Lineage%in% c("CAMA1 RiboR","CAMA1 Sen" ) ]
egdatTx0$Day<-as.numeric(as.character(egdatTx0$day)); egdatTx0[Day==1, Day:=0]
ggplot(egdatTx0[Lineage%in% c("CAMA1 RiboR")],
       aes(x=Hour/24, y=NormSensExpression.SSGSEAB, fill=Treatment,col= Treatment, group=interaction(State,Gene_Set,Treatment) ))+
  geom_vline(xintercept=24*c(0)/24, linetype=2,col="grey")+
  theme_classic(base_size=26)+
  geom_point(pch=21,size=6.5, aes(alpha=1/(1+(Day))))+ geom_line()+
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
   theme(aspect.ratio=1)+labs(title="Combined tx effects", y="Pathway modulation
                             \n (relative to DMSO)", x="Day")+
  scale_color_discrete(name="Treatment",labels=c("DMSO","Afatinib","Ribociclib","Combination") )+
  scale_fill_discrete(name="Treatment",labels=c("DMSO","Afatinib","Ribociclib","Combination") )+
  scale_x_continuous(breaks=c(0,6,12,18))+
  scale_alpha( range = c(0.4, 1),guide = "none")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic tx effects2.png", width=8, height=8, dpi=320)



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
  scale_color_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  scale_fill_npg(name="Cell type",labels=c("Resistant","Sensitive") )+
  theme(aspect.ratio=1)+labs(title="Plasticity of tx responses \n           (Day 0-18)", y="Tx response (activation/day)", x="Day")+
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  scale_alpha( range = c(0.4, 1),guide = "none")+theme(legend.position="none")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Schematic input data P3.png", width=7, height=7, dpi=320)




