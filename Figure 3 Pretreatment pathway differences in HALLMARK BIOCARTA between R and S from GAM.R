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

# Specify source data location 
fileloc <- ("~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/")
setwd(fileloc)

# Load inferences from GAM model of time course RNAseq data for CAMA-1 and MCF-7
nsheets<-length( excel_sheets( "sh10050_MCF7_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx" ) )

dd0mcf <-rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx",sheet=sheet)
}))

dd0cama<-rbindlist(lapply(1:nsheets,function(sheet){
  read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_parametric_table_synergy.xlsx",sheet=sheet)
}))

# join GAM data for CAMA-1 and MCF-7
dd0Res <- data.table(rbind(data.table(CellLine="CAMA1",dd0cama%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )),
                           data.table(CellLine="MCF7", dd0mcf%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ))) )[!grepl("_DN",Gene_Set)][ 
                             grepl("BIOCARTA_",Gene_Set)|grepl("HALLMARK_",Gene_Set) |grepl("KEGG_",Gene_Set)|grepl("REACTOME_",Gene_Set)]

# list unique pathways
unique(dd0Res[   grepl("BIOCARTA_",Gene_Set)|grepl("HALLMARK_",Gene_Set)  ]$Gene_Set)

# set order of treatment, cell line and resistance state factors
dd0Res$Treatment <- factor(dd0Res$Treatment,c("DMSO", "AfatMain","RiboMain", "CombSynergy"))
dd0Res$State <- factor(dd0Res$State,c("Sens","RiboR"))
dd0Res$CellLine <- factor(dd0Res$CellLine,c("CAMA1", "MCF7"))
#dd0Res[, Day:=(Hour-as.numeric(as.character(timeunderdrug)))/24 ]

# list the different parameters of the model and subset those linked to pre-treatment pathway differences between resistant and sensitive cells
parlist<-dd0Res$p.covariate%>%unique()
parlist[!grepl("CumDay",parlist)& !grepl("isHour",parlist)]

# Subset the differences between resistant and sensitive cells in Hallmark and Biocarta pathways with FDR<0.05 and order by effect size
plotThis <-dd0Res[grepl("HALLMARK",Gene_Set)|grepl("BIOCARTA",Gene_Set)][p.covariate=="StateRiboR" ][FDR <0.05][order(-Estimate)]
# Generate pathway label annotation column
plotThis$Gene_SetLab <- gsub("_"," ",plotThis$Gene_Set)
# Summarize mean difference in pathway activity between resistant and sensitive cells when comparing across cell lines and the consistency of the direction of the effect
orderplot<-data.table(plotThis %>% group_by(Gene_SetLab)%>%summarise(mu=mean(Estimate),consistdir=-1+sum(Estimate>0) ))[order(mu)]
orderplot<-data.table(dd0Res[grepl("HALLMARK",Gene_Set)|grepl("BIOCARTA",Gene_Set)][p.covariate=="StateRiboR" ] %>% group_by(Gene_Set)%>%summarise(mu=mean(Estimate),consistdir=-1+sum(Estimate>0) ))[order(mu)]
# add pathway label annotation column and set order 
orderplot$Gene_SetLab <- gsub("_"," ",orderplot$Gene_Set)
plotThis$Gene_SetLab<- factor(plotThis$Gene_SetLab, levels=orderplot$Gene_SetLab)

# Select Top 100 increased and decreased pathways in resistant cells
nplot<-100
wchplot<-orderplot[c(1:nplot,((nrow(orderplot)-(nplot-1)):nrow(orderplot)))]$Gene_SetLab
subplotThis <- plotThis[Gene_SetLab%in%wchplot][Gene_SetLab%in%orderplot[consistdir!=0]$Gene_SetLab  ]
subplotThis <- plotThis[][Gene_SetLab%in%orderplot[consistdir!=0]$Gene_SetLab  ]
subplotThis$Gene_SetLab<- factor(subplotThis$Gene_SetLab, levels=unique(wchplot))
subplotThis[,nlines:=length(FDR), by=Gene_Set]
subplotThis[nlines==2]
#write.csv(subplotThis[nlines==2] , file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/tables/Pre tx resistant cell pathway activity.csv")

# plot differential activity of the 100 increased and decreased in resistant cells
ggplot(subplotThis , aes(y=Estimate, x=Gene_SetLab, col=CellLine,
                         group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=3)+coord_flip()+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2.5)+ labs(y="Resistant cell pathway activition \n (pre-tx compared to sensitive cells)", x="Gene set")+
  scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF-7"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Resistant initial diff R vs S HALLMARK BIOCARTA.png", width=16, height=16, dpi=320)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Resistant initial diff R vs S HALLMARK BIOCARTA.pdf", width=16, height=16, dpi=320)

# plot differential activity of specific cell cycle and MYC pathways in resistant cells
subplotThisB<-subplotThis[Gene_Set %in%c( "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1" ,"HALLMARK_MYC_TARGETS_V2" )] 
subplotThisB[,Gene_SetLab:=gsub("HALLMARK","HALLMARK \n",Gene_SetLab)]
subplotThisB$Gene_SetLab<- factor(subplotThisB$Gene_SetLab, levels= c("HALLMARK \n MYC TARGETS V2",  "HALLMARK \n MYC TARGETS V1","HALLMARK \n G2M CHECKPOINT", "HALLMARK \n E2F TARGETS"))
subplotThisB[, ucl:=Estimate+1.96*Std.Error]
subplotThisB[, lcl:=Estimate-1.96*Std.Error]
ggplot(subplotThisB, aes(y=Estimate, x=Gene_SetLab, col=CellLine,
                         group=interaction(Gene_SetLab,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_errorbar(aes(ymax=ucl, ymin=lcl),width=0.3,position=position_dodge(width=0.5))+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=4)+coord_flip()+
  theme_classic(base_size=26)+
  theme(aspect.ratio=1)+ labs(y="Resistant cell pathway activition \n (pre-tx compared to sensitive cells)", x="Gene set")+
  scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF7"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Simplified Resistant initial diff R vs S HALLMARKcellcycleMYC.png", width=10, height=10, dpi=320)



# load timecourse rnaseq data
list.files(fileloc)
# If first time using code, translate .xlsx file to .csv format for ease of reading 
#dd0mcf <-read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx")
#dd0cama<-read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx")
#dd0mcf <-rbindlist(lapply(1:nsheets,function(sheet){
#  read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx",sheet=sheet)
#}))
#dd0cama<-rbindlist(lapply(1:nsheets,function(sheet){
# read.xlsx(xlsxFile = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.xlsx",sheet=sheet)
#}))
#write.csv(dd0cama,file="sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
#write.csv(dd0mcf,file="sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
#read.xlsx(xlsxFile = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_smooths_table_synergy.xlsx")
dd0cama <- fread(file="sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
dd0mcf <- fread(file="sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")

# Join CAMA-1 and MCF-7 data
dd0 <- data.table(rbind(dd0cama %>% dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ),
                        dd0mcf %>% dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )))
# Annotate lineage and set order of factor names
dd0[,Lineage:=paste0(CellLine," ",State)]
dd0$Lineage <- factor(dd0$Lineage  , levels=c("CAMA1 Sen" ,"MCF7 Sen", "CAMA1 RiboR"  ,  "MCF7 RiboR"  ))

# Extract pre-treatment data for pathways showing differential activity between resistant and sensitive cells that is consistent across cell lines  
dd0lu <- dd0[Treatment=="DMSO"][Gene_Set%in%subplotThis[nlines==2]$Gene_Set][day==1][][Hour==0]
# Scale and center ssGSEA data fro treated samples
dd0lu[,scaleSSGSEA:= scale(Enrichment_Score), by=c("Gene_Set","CellLine","State")]
dd0lu[,scaleSSGSEAB:= scale(Enrichment_Score), by=c("Gene_Set","CellLine")]
# Contrast pathway activity to that of sensitive cell  
dd0lu[,SensExpression.SSGSEA := sum(Enrichment_Score*(State=="Sen"&Treatment=="DMSO")), by=c("Gene_Set","CellLine","Hour")]
dd0lu[,SensExpression.SSGSEAB:= sum(Enrichment_Score*(Treatment=="DMSO" )), by=c("Gene_Set","CellLine","Hour","State")]
dd0lu[,NormSensExpression.SSGSEA:=Enrichment_Score -SensExpression.SSGSEA, by=c("Gene_Set","CellLine","Hour", "Treatment")]
dd0lu[,NormSensExpression.SSGSEAB:=Enrichment_Score -SensExpression.SSGSEAB, by=c("Gene_Set","CellLine","Hour", "Treatment")]

# order factors
dd0lu$Treatment <- factor(dd0lu$Treatment,c("DMSO", "Afat","Ribo", "Comb"))
dd0lu$Gene_SetLab <- gsub("_"," ",dd0lu$Gene_Set)
dd0lu$Gene_SetLab<- factor(dd0lu$Gene_SetLab, levels=unique(wchplot))

# Extract data for resistant cells as we are comparing to sensitive cells
# then summarise the differential expression and use this to order pathways
ordplot2<-  data.table(dd0lu[Lineage%in%c("CAMA1 RiboR","MCF7 RiboR")] %>% group_by(Gene_SetLab)%>%summarise(mu=mean(NormSensExpression.SSGSEA) ) )[order(mu)]

# Identify the top 15 and visualize
nplot<-15
wchplot2<-ordplot2[c(1:nplot,((nrow(ordplot2)-(nplot-1)):nrow(ordplot2)))]$Gene_SetLab
dd0lu2 <-dd0lu[Gene_SetLab%in%wchplot2]
dd0lu2$Gene_SetLab<- factor(dd0lu2$Gene_SetLab, levels=ordplot2$Gene_SetLab)

ggplot(dd0lu2[State!="Sen"] , aes(y=NormSensExpression.SSGSEA, x=Gene_SetLab, col=CellLine,
                                 group=interaction(Gene_SetLab ,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=5)+coord_flip()+
  theme_classic(base_size=28)+
  theme(aspect.ratio=1.6)+ labs(y="Pre-tx resistant cell pathway activity \n (compared to pre-tx sensitive cells)", x="Gene set")+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))+
  scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF-7"))+
  theme(
    axis.text = element_text(color="black"),
    axis.ticks = element_line(color = "black")
  )
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Resistant initial diff compare R vs S HALLMARK BIOCARTA.png", width=16, height=16, dpi=320)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Fig3 Resistant initial diff compare R vs S HALLMARK BIOCARTA.pdf", width=16, height=16, dpi=320)


ggplot(dd0lu2[State!="Sen"] , aes(y=NormSensExpression.SSGSEA, x=Gene_SetLab, col=CellLine,
                                  group=interaction(Gene_SetLab ,CellLine)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(aes(shape=CellLine),position=position_dodge(width=0.5),size=5)+coord_flip()+
  theme_classic(base_size=28)+
  theme(aspect.ratio=1.6)+ labs(y="Pre-tx resistant cell pathway activity \n (compared to pre-tx sensitive cells)", x="Gene set")+
  #scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF7"))+#,values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF-7"))+
  scale_color_jco(name="Cell line", labels=c("CAMA-1", "MCF-7"))+
  theme(
    text = element_blank(),legend.text = element_blank(),legend.title = element_blank(),
    axis.ticks = element_line(color = "black")
  )
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/BLANK_Fig3 Resistant initial diff compare R vs S HALLMARK BIOCARTA.png", width=12, height=12, dpi=320)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/BLANK_Fig3 Resistant initial diff compare R vs S HALLMARK BIOCARTA.pdf", width=16, height=16, dpi=320)


#range(subplotThis[nlines==2][Gene_SetLab%in% unique(dd0lu2[State!="Sen"]$Gene_SetLab)]$FDR)
# Sample size assessments for figure legends
# dd0[Gene_Set=="HALLMARK_ADIPOGENESIS"][day==1][Hour==0][Treatment!="DMSO"]
# dd0[,SamplDuplic:=0]
# dd0[day==1&Hour==0&Treatment!="DMSO",SamplDuplic:=1]
# dd0[Gene_Set=="HALLMARK_ADIPOGENESIS"][SamplDuplic==0]
# dd0[Gene_Set=="HALLMARK_ADIPOGENESIS"][SamplDuplic==0][day==1]%>%nrow()
# 180
# 3*2*
# # 2 paired resistant/sensitive cell lines, 4 Treatments measured on 4 days at 3 hours
# # 2 paired resistant/sensitive cell lines, 
# # 3 drug treatments measured on first day of treatment at 2 hour marks and on 3 later days at 3 hour marks
# # DMSO treatment measured on 4 day 3 hour marks
# 4*(3*1*2+
#   3*3*3+
#   1*4*3)
# 4*(3*1*2+1*1*3)
# dd0Res$Gene_Set%>%unique()%>%length()
# dd0Res$Hour%>%unique() 
# dd0$day%>%unique() 


