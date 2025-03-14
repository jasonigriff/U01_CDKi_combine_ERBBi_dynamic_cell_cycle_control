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
fileloc <- "~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/"
fileloc <- ("~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/")
list.files(fileloc)
setwd(fileloc)
dd0mcf <-read.csv(file = "sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
dd0cama<-read.csv(file = "sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")

dd0 <- data.table(rbind(dd0cama%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ),
                        dd0mcf%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )))

dd0[,Lineage:=paste0(CellLine," ",State)]
dd0$Lineage <- factor(dd0$Lineage  , levels=c("CAMA1 Sen" ,"CAMA1 RiboR"  ,"MCF7 Sen",   "MCF7 RiboR"  ))

dd0[ grepl("BIOCARTA_",Gene_Set)|grepl("HALLMARK_",Gene_Set)|grepl("KEGG_",Gene_Set)|grepl("REACTOME_",Gene_Set)]$Gene_Set%>%unique()

dd0[ grepl("APOPTOSIS",Gene_Set)|grepl("DEATH",Gene_Set)|grepl("CASPASE",Gene_Set)]$Gene_Set%>%unique()
apoptosis_pathways<-dd0[ grepl("BIOCARTA_",Gene_Set)|grepl("HALLMARK_",Gene_Set)|grepl("KEGG_",Gene_Set)|grepl("REACTOME_",Gene_Set)][ 
  grepl("APOPTOSIS",Gene_Set)|grepl("DEATH",Gene_Set)|grepl("CASPASE",Gene_Set)]$Gene_Set%>%unique()

lu <- data.table(Gene_Set=c("BIOCARTA_DEATH_PATHWAY",
                            "BIOCARTA_MITOCHONDRIA_PATHWAY",
                            "KEGG_APOPTOSIS",
                            "HALLMARK_APOPTOSIS",
                            #"HALLMARK_DNA_REPAIR",
                            "REACTOME_APOPTOSIS_INDUCED_DNA_FRAGMENTATION",
                            "REACTOME_APOPTOSIS",
                            "REACTOME_DEATH_RECEPTOR_SIGNALING",
                            "REACTOME_PROGRAMMED_CELL_DEATH",
                            "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CASPASE_ACTIVATORS_AND_CASPASES",
                            "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
                            "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DEATH_RECEPTORS_AND_LIGANDS",
                            "REACTOME_CASPASE_ACTIVATION_VIA_DEATH_RECEPTORS_IN_THE_PRESENCE_OF_LIGAND",
                            "REACTOME_CASPASE_ACTIVATION_VIA_DEPENDENCE_RECEPTORS_IN_THE_ABSENCE_OF_LIGAND",
                            "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
                            "REACTOME_CASPASE_MEDIATED_CLEAVAGE_OF_CYTOSKELETAL_PROTEINS",
                            "REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES",
                            "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
                            "REACTOME_NRIF_SIGNALS_CELL_DEATH_FROM_THE_NUCLEUS"))

dd0lu <- dd0[Gene_Set%in%lu$Gene_Set]

dd0lu[,scaleSSGSEA:= scale(Enrichment_Score), by=c("Gene_Set","CellLine","State")]
dd0lu[,scaleSSGSEAB:= scale(Enrichment_Score), by=c("Gene_Set","CellLine")]

dd0lu[,SensExpression.SSGSEA := sum(Enrichment_Score*(State=="Sen"&Treatment=="DMSO")), by=c("Gene_Set","CellLine","Hour")]
dd0lu[,SensExpression.SSGSEAB:= sum(Enrichment_Score*(Treatment=="DMSO" )), by=c("Gene_Set","CellLine","Hour","State")]
dd0lu[,NormSensExpression.SSGSEA:=Enrichment_Score -SensExpression.SSGSEA, by=c("Gene_Set","CellLine","Hour", "Treatment")]
dd0lu[,NormSensExpression.SSGSEAB:=Enrichment_Score -SensExpression.SSGSEAB, by=c("Gene_Set","CellLine","Hour", "Treatment")]
dd0lu$Treatment <- factor(dd0lu$Treatment,c("DMSO", "Afat","Ribo", "Comb"))
dd0lu$Hour%>%unique()
dd0lu$Gene_Set%>%unique()

FinDay<-dd0lu[#Treatment!="DMSO"][
  Hour>=(max(Hour)-24)]
tmpdd<-FinDay[Gene_Set=="HALLMARK_APOPTOSIS"]

tmp<-lmer(scaleSSGSEAB~1+RiboTreated +AfatTreated+ CombTreated +(1|Hour)+(1|Lineage) ,data=FinDay[Gene_Set=="HALLMARK_APOPTOSIS"])
#tmp<-lmer(scaleSSGSEAB~1+RiboTreated +AfatTreated+ CombTreated +(1|Hour)+(1+RiboTreated +AfatTreated+ CombTreated|Lineage) ,data=FinDay[Gene_Set=="HALLMARK_APOPTOSIS"])
tmpdd$pred<-predict(tmp,re.form=NA)
ggplot(tmpdd,
       aes(x=Treatment,  y=scaleSSGSEAB ))+
  theme_classic(base_size=6)+
  geom_point()+geom_point(size=5,aes(y=pred))


tmp<-lmer(NormSensExpression.SSGSEAB~-1+RiboTreated +AfatTreated+ CombTreated +(1|Hour)+
            (0+RiboTreated|Lineage) +
            (0+AfatTreated|Lineage) +  (0+CombTreated|Lineage)  ,data=FinDay[Gene_Set=="HALLMARK_APOPTOSIS"])
tmp<-lmer(NormSensExpression.SSGSEAB~-1+RiboTreated +AfatTreated+ CombTreated +(1|Hour)+(1|Lineage) ,data=FinDay[Gene_Set=="HALLMARK_APOPTOSIS"])
#tmp<-lmer(scaleSSGSEAB~1+RiboTreated +AfatTreated+ CombTreated +(1|Hour)+(1+RiboTreated +AfatTreated+ CombTreated|Lineage) ,data=FinDay[Gene_Set=="HALLMARK_APOPTOSIS"])
tmpdd$pred<-predict(tmp,re.form=NA)
tmpdd$pred<-predict(tmp)
ggplot(tmpdd,
       aes(x=Treatment,  y=NormSensExpression.SSGSEAB ))+
  theme_classic(base_size=6)+
  geom_point()+geom_point(size=5,aes(y=pred))

summary(tmp)


finFULL2<- rbindlist(mclapply(dd0lu$Gene_Set%>%unique(), function(x){
  pred<-FinDay[Gene_Set== x];  cat(x)
  if(length(unique(pred$CellLine))>1){
    m1 <- lmer(NormSensExpression.SSGSEAB~ -1 +  RiboTreated + AfatTreated + CombTreated +
                 (0+RiboTreated|Lineage) +
                 (0+AfatTreated|Lineage) +
                 (0+CombTreated|Lineage) +
                 (1|Hour) , data=pred)
    outmat <- cbind(Gene_Set=x,Contrast="All_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )
    return(outmat)
  }else{ NULL }
}, mc.cores =detectCores()-1 ))
setnames(finFULL2,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
fincomp2 <- finFULL2[rn%in%c("CombTreated")]
ggplot(finFULL2,aes(x=Estimate, y=pvalue, col=pvalue<0.05))+geom_point()+
  facet_wrap(~rn)
finFULL2[pvalue<0.05]$Gene_Set
finFULL2[pvalue<0.05][rn=="CombTreated"][order(Estimate)]
orderbyt<-finFULL2[rn=="CombTreated"][order(Estimate)]$Gene_Set
#orderbyt<-finFULL2[rn=="CombTreated"][order(tvalue)]$Gene_Set
finFULL2$Gene_Set<- factor( finFULL2$Gene_Set, levels=orderbyt)
finFULL2$Gene_SetLab<- gsub("_"," ",finFULL2$Gene_Set)
finFULL2$Gene_SetLab<- factor( finFULL2$Gene_SetLab , levels=gsub("_"," ",orderbyt))
finFULL2$rn<- factor(finFULL2$rn,  levels=c("RiboTreated","AfatTreated","CombTreated"  ))


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ggplot(finFULL2,aes(x=Estimate, y=Gene_Set, col=rn))+geom_point()+
  labs(x="Estimate",y="Pathway")
ggplot(finFULL2,aes(y=Estimate,x=Gene_SetLab,fill=rn))+  theme_classic(base_size=26)+
  #geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  labs(y="Estimate",x="Pathway")+
  #facet_wrap(~CellLine)+
  coord_flip() +
  theme(aspect.ratio=1)+
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSetMainTreatmentAndSynergyByEnd.pdf", dpi=320, height=6, width=20)

finFULL2[rn=="CombTreated"]
#write.csv(finFULL2[rn=="CombTreated"], file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/SI stats pathway apoptosis synergy.csv")

FinDay$Gene_Set<- factor( FinDay$Gene_Set, levels=orderbyt)

dd0lu$Gene_SetLab<- gsub("_"," ",dd0lu$Gene_Set)
dd0lu$Gene_SetLab<- factor( dd0lu$Gene_SetLab , levels=gsub("_"," ",orderbyt))
dd0lu[,CelltypeScaledNormSensExpression.SSGSEAB:=scale(NormSensExpression.SSGSEAB), by= c("Gene_Set","CellLine","State")]
dd0lu[   ,LineageLab:="CAMA-1 \n Sensitive"]
dd0lu[CellLine=="MCF7"& State=="Sen",LineageLab:="MCF-7 \n Sensitive"]
dd0lu[CellLine=="CAMA1"& State=="RiboR",LineageLab:="CAMA-1 \n Resistant"]
dd0lu[CellLine=="MCF7"& State=="RiboR",LineageLab:="MCF-7 \n Resistant"]
dd0lu$LineageLab<- factor(dd0lu$LineageLab, levels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant" ,"MCF-7 \n Sensitive","MCF-7 \n Resistant"  ))
dd0lu[,Treatmentlab:="Combination"]
dd0lu[Treatment=="Afat",Treatmentlab:="Afatinib"]
dd0lu[Treatment=="Ribo",Treatmentlab:="Ribociclib"]
dd0lu[Treatment=="DMSO",Treatmentlab:="DMSO"]
dd0lu$Treatmentlab <- factor(dd0lu$Treatmentlab,  levels=c("DMSO" ,"Ribociclib","Afatinib","Combination"  ))

# apopcor<-tidyr::spread(dd0lu%>%
#                 select(Lineage,Treatment,Hour,Gene_Set, NormSensExpression.SSGSEAB),
#               Gene_Set, NormSensExpression.SSGSEAB)%>%
#   select(-c(Lineage,Treatment,Hour))%>%
#   cor()
# unname(apopcor)%>%
#   corrplot::corrplot.mixed()
# 
# ggplot(tidyr::spread(dd0lu[Treatment!="DMSO"]%>%
#                        select(Lineage,Treatment,Hour,Gene_Set, NormSensExpression.SSGSEAB),
#                      Gene_Set, NormSensExpression.SSGSEAB),
#        aes(BIOCARTA_MITOCHONDRIA_PATHWAY,REACTOME_APOPTOSIS,col=Lineage))+geom_point()

ggplot(dd0lu[Treatment!="DMSO"],#[Gene_Set=="HALLMARK_APOPTOSIS"],
       aes(x=as.factor(Hour), y=Gene_SetLab, fill=CelltypeScaledNormSensExpression.SSGSEAB ))+
  theme_classic(base_size=16)+
  geom_tile()+
  facet_grid(Treatmentlab~LineageLab,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_fill_gradient2(name="Gene set enrichment \n scaled(compared to DMSO)" , low="blue", high="darkred", mid="white")+
  theme(aspect.ratio=1)+#labs(title=gsub("_"," ",x), x="Hour")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")+
  scale_y_discrete(guide=guide_axis(n.dodge=1)) +
  labs(x="Time (hours)",y="Pathway")
ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/ApoptosisGeneSetsOverTime.png", width=20, height=8, dpi=320)




gamOrd<-c("REACTOME APOPTOSIS","REACTOME PROGRAMMED CELL DEATH",                                        
 "REACTOME TP53 REGULATES TRANSCRIPTION OF CELL DEATH GENES",             
 "REACTOME INTRINSIC PATHWAY FOR APOPTOSIS",                              
 "KEGG APOPTOSIS",                                                        
 "REACTOME NRIF SIGNALS CELL DEATH FROM THE NUCLEUS",                     
 "HALLMARK APOPTOSIS",                                                    
 "REACTOME FOXO MEDIATED TRANSCRIPTION OF CELL DEATH GENES",              
 "BIOCARTA MITOCHONDRIA PATHWAY",                                         
 "REACTOME CASPASE ACTIVATION VIA EXTRINSIC APOPTOTIC SIGNALLING PATHWAY",
 "BIOCARTA DEATH PATHWAY")           

dd0lu$Gene_SetLab<- factor( dd0lu$Gene_SetLab , levels=gsub("_"," ",gamOrd))


ggplot(dd0lu[Treatment!="DMSO"],#[Gene_Set=="HALLMARK_APOPTOSIS"],
       aes(x=as.factor(Hour), y=Gene_SetLab, fill=CelltypeScaledNormSensExpression.SSGSEAB ))+
  theme_classic(base_size=16)+
  geom_tile()+
  facet_grid(Treatmentlab~LineageLab,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_fill_gradient2(name="Gene set enrichment \n scaled(compared to DMSO)" , low="blue", high="darkred", mid="white")+
  theme(aspect.ratio=1)+#labs(title=gsub("_"," ",x), x="Hour")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")+
  scale_y_discrete(guide=guide_axis(n.dodge=1)) +
  labs(x="Time (hours)",y="Pathway")
ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/ApoptosisGeneSetsOverTimeGAMordered.png", width=20, height=8, dpi=320)








# 
# 
# ggplot(dd0lu,#[Gene_Set=="HALLMARK_APOPTOSIS"],
#        aes(x=as.factor(Hour), y=Gene_SetLab, fill=scaleSSGSEAB ))+
#   theme_classic(base_size=6)+
#   geom_tile()+
#   facet_grid(Treatment~LineageLab,scales="free")+
#   #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
#   scale_fill_gradient2(name="Gene set enrichment \n (compared to DMSO)" , low="blue", high="darkred", mid="white")+
#   theme(aspect.ratio=3)+#labs(title=gsub("_"," ",x), x="Hour")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#   scale_y_discrete(guide=guide_axis(n.dodge=1)) +
#   labs(x="Time (hours)",y="Pathway")
# 
# 
ggplot(FinDay[Treatment!="DMSO"]%>%group_by(Gene_Set,Treatment,Lineage)%>%
         summarise(NormSensExpression.SSGSEAB= mean(NormSensExpression.SSGSEAB))
       ,#[Gene_Set=="HALLMARK_APOPTOSIS"],
       aes(x=Lineage, y=Gene_Set, fill=NormSensExpression.SSGSEAB ))+
  theme_classic(base_size=6)+
  geom_tile()+
  facet_grid(~Treatment,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_fill_gradient2(name="Expression \n log2FC \n vs \n sensitive cells \n under DMSO" , low="blue", high="darkred", mid="white")+
  theme(aspect.ratio=3)+#labs(title=gsub("_"," ",x), x="Hour")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_y_discrete(guide=guide_axis(n.dodge=1))

# 
# ggplot(FinDay,#[Gene_Set=="HALLMARK_APOPTOSIS"],
#        aes(x=as.factor(Hour), y=Gene_Set, fill=NormSensExpression.SSGSEAB ))+
#   theme_classic(base_size=6)+
#   geom_tile()+
#   facet_grid(Treatment~Lineage,scales="free")+
#   #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
#   scale_fill_gradient2(name="Expression \n log2FC \n vs \n sensitive cells \n under DMSO" , low="blue", high="darkred", mid="white")+
#   theme(aspect.ratio=3)+#labs(title=gsub("_"," ",x), x="Hour")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#   scale_y_discrete(guide=guide_axis(n.dodge=1)) 
# 
# ggplot(dd0lu[Treatment!="DMSO"],
#        aes(x=as.factor(Hour), y=Gene_Set, fill=scaleSSGSEAB ))+
#   theme_classic(base_size=6)+
#   geom_tile()+
#   facet_grid(Treatment~Lineage,scales="free")+
#   #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
#   scale_fill_gradient2(name="Expression \n log2FC \n vs \n sensitive cells \n under DMSO" , low="blue", high="darkred", mid="white")+
#   theme(aspect.ratio=3)+#labs(title=gsub("_"," ",x), x="Hour")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#   scale_y_discrete(guide=guide_axis(n.dodge=1)) 
# 
# 
# 
# 
# 
# 
# 
# p1<-ggplot(dd0lu[Treatment!="DMSO"],
#            aes(x=as.factor(Hour), y=Gene_Set, fill=NormSensExpression.SSGSEAB ))+
#   theme_classic(base_size=6)+
#   geom_tile()+
#   facet_grid(Treatment~Lineage,scales="free")+
#   #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
#   scale_fill_gradient2(name="Expression \n log2FC \n vs \n sensitive cells \n under DMSO" , low="blue", high="darkred", mid="white")+
#   theme(aspect.ratio=3)+#labs(title=gsub("_"," ",x), x="Hour")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#   scale_y_discrete(guide=guide_axis(n.dodge=1)) 
# p1
# 
# 
# 
# 
# p1<- ggplot(dd0lu[Treatment!="DMSO"][
#   Gene_Set%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1" )],
#   aes(x=as.factor(Hour), y=gsub("_"," ",Gene_Set), fill=NormSensExpression.SSGSEA ))+
#   theme_classic(base_size=16)+
#   geom_tile()+
#   facet_grid(Treatment~Lineage,scales="free")+
#   #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
#   scale_fill_gradient2(name="ssGSEA pathway activity \n vs \n sensitive cells \n under DMSO" , low="blue", high="darkred", mid="white")+
#   theme(aspect.ratio=3)+labs(y="Pathway", x="Hour")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#   scale_y_discrete(guide=guide_axis(n.dodge=1)) 
# p1
