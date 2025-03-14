rm(list=ls())
require(msigdbr);require(data.table); require(dplyr); require(ggplot2);require(ggsci);require(colorspace)
require(R.utils)
require(umap)
require(Rdimtools)
require(openxlsx)
require(lme4)
require(lmerTest)
require(parallel)
require(scales)

# load timecourse rnaseq data
fileloc <- ("~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse")
setwd(fileloc)
dd0mcf <-read.csv(file = "sh10050_MCF7_genes_gam_short.term_inputData_synergy.csv")
dd0cama<-read.csv(file = "sh11141_CAMA1_genes_gam_short.term_inputData_synergy.csv")
dd0 <- data.table(rbind(dd0cama%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ),dd0mcf%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )))
dd0$Treatment <- factor(dd0$Treatment,c("DMSO", "Afat","Ribo", "Comb"))
dd0$State <- factor(dd0$State,c("Sen","RiboR"))
dd0$CellLine <- factor(dd0$CellLine,c("CAMA1", "MCF7"))
dd0[, Day:=(Hour-as.numeric(as.character(timeunderdrug)))/24 ]
# Scale each genes expression data for each cell line (line & resistance)
dd0[,scalelog2cpm:= scale(Expression.log2cpm), by=c("Gene","CellLine","State")]
dd0[,scalelog2cpmB:= scale(Expression.log2cpm), by=c("Gene","CellLine")]

# Scale cell line expression of gene g relative to sensitive cells at time t under DMSO
dd0[,SensExpression.log2cpm := sum(Expression.log2cpm*(State=="Sen"&Treatment=="DMSO")), by=c("Gene","CellLine","Hour")]
dd0[,SensExpression.log2cpmB:= sum(Expression.log2cpm*(Treatment=="DMSO" )), by=c("Gene","CellLine","Hour","State")]
dd0[,NormSensExpression.log2cpm:=Expression.log2cpm -SensExpression.log2cpm, by=c("Gene","CellLine","Hour", "Treatment")]
dd0[,NormSensExpression.log2cpmB:=Expression.log2cpm -SensExpression.log2cpmB, by=c("Gene","CellLine","Hour", "Treatment")]

dd0[,Lineage:=paste0(CellLine," ",State)]
dd0$Lineage <- factor(dd0$Lineage  , levels=c("CAMA1 Sen" ,"MCF7 Sen", "CAMA1 RiboR"  ,  "MCF7 RiboR"  ))
# list all genes
selectgenes<- unique(dd0$Gene)


#load ssGSEA pathwaygeneset metadata
ssGSEAmeta <- data.table( msigdbr(species = "Homo sapiens") )
ssGSEAmeta[gs_cat%in%c("C2","H")][gs_subcat%in%c("","CP:BIOCARTA","CP:REACTOME","CP:KEGG" )]$gs_name%>%unique()
ssGSEAmeta[gs_cat%in%c("C2","H")][gs_subcat%in%c("CP:BIOCARTA")]$gs_name%>%unique()
ssGSEAmeta[,ngenes:=length(unique(gene_symbol)), by=gs_name]

gamOrd<-c("REACTOME APOPTOSIS",
          "REACTOME PROGRAMMED CELL DEATH",                                        
 "REACTOME TP53 REGULATES TRANSCRIPTION OF CELL DEATH GENES",             
 "REACTOME INTRINSIC PATHWAY FOR APOPTOSIS",                              
 "KEGG APOPTOSIS",                                                        
 "REACTOME NRIF SIGNALS CELL DEATH FROM THE NUCLEUS",                     
 "HALLMARK APOPTOSIS",                                                    
 "REACTOME FOXO MEDIATED TRANSCRIPTION OF CELL DEATH GENES",              
 "BIOCARTA MITOCHONDRIA PATHWAY",                                         
 "REACTOME CASPASE ACTIVATION VIA EXTRINSIC APOPTOTIC SIGNALLING PATHWAY",
 "BIOCARTA DEATH PATHWAY")           

apopGSmeta <- ssGSEAmeta[gs_name%in%gsub(" ","_",gamOrd)][]
apopGSmeta$gs_name <- factor(apopGSmeta$gs_name,levels=gsub(" ","_",gamOrd))
apopGSmeta <- apopGSmeta[order(gs_name)]
apopGSmeta$gs_name%>%unique()
apopGSmeta$gs_description%>%unique()

apopGSmeta[gs_name=="BIOCARTA_DEATH_PATHWAY"]$gene_symbol
apopGSmeta[gs_name=="REACTOME_APOPTOSIS"]$gene_symbol

presabsMatdd<- unique(apopGSmeta%>%
         dplyr::select(gs_name,gene_symbol))%>%
  mutate(present=1)%>%
  tidyr::spread(gs_name,present,fill=0)
presabsMat<-as.matrix(presabsMatdd[,-1])
rownames(presabsMat)<-presabsMatdd[,1]
colnames(presabsMat)<- substr(colnames(presabsMat), start = 1, stop = 15)
presabs <- data.table(unique(apopGSmeta%>%
     dplyr::select(gs_name,gene_symbol))%>%
     mutate(present=1)%>%
     tidyr::spread(gs_name,present,fill=0)%>%
     tidyr::gather(gs_name,present,-gene_symbol))
presabs$gs_name <- factor(presabs$gs_name,levels=gsub(" ","_",gamOrd))

pheatmap::pheatmap(t(presabsMat),cex=0.5,cluster_rows=F,
                   scale="none",
                   clustering_distance_cols= "binary",
                   clustering_distance_rows= "binary")
pheatmap::pheatmap(t(presabsMat),cex=0.45,cluster_rows=T,
                   scale="none",
                   clustering_distance_cols= "minkowski",
                   clustering_distance_rows= "minkowski")
pheatmap::pheatmap(t(presabsMat[,3:11]),cex=0.45,cluster_rows=T,
                   scale="none",
                   clustering_distance_cols= "minkowski",
                   clustering_distance_rows= "minkowski")

library(UpSetR)
library("ComplexHeatmap")
library("enrichplot")

#cnetplot
ll1<-lapply(unique(apopGSmeta$gs_name), function(i){apopGSmeta[gs_name==i]$gene_symbol})
names(ll1)<-gsub("_"," ",unique(apopGSmeta$gs_name))
set.seed(1234)
cnetplot(ll1, node_label="category",
         layout="fr",
         cex.params = list(category_node = 0.25,category_label = 0.35,gene_label=0.1),#
         showCategory =length(unique(apopGSmeta$gs_name)))

#UpSet plot
m1<-make_comb_mat(presabsMat,mode ="distinct" )
ss1 = set_size(m1)
cs1 = comb_size(m1)
UpSet(m1,  set_order = order(-ss1),
      comb_order = order(comb_degree(m1), -cs1))

minimaldd<-unique(apopGSmeta %>% dplyr::select(gs_name,gene_symbol))


REACTAPandPROGgenes <-data.table(presabsMatdd)[REACTOME_APOPTOSIS==1     & REACTOME_PROGRAMMED_CELL_DEATH==1 &
                        # REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES==1&
                         REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS==0      &                       
                         KEGG_APOPTOSIS==0                       &                                  
REACTOME_NRIF_SIGNALS_CELL_DEATH_FROM_THE_NUCLEUS==0   &                   
HALLMARK_APOPTOSIS==0                   &                                  
 REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES==0   &            
BIOCARTA_MITOCHONDRIA_PATHWAY==0         &                          
REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY==0& 
BIOCARTA_DEATH_PATHWAY==0    
                         ]$gene_symbol

REACT_TP53genes <-data.table(presabsMatdd)[REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES==1     &
                                                 #REACTOME_APOPTOSIS==0& 
                                                 #REACTOME_PROGRAMMED_CELL_DEATH== 0&
                                                 REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS==0      &                       
                                                 KEGG_APOPTOSIS==0                       &                                  
                                                 REACTOME_NRIF_SIGNALS_CELL_DEATH_FROM_THE_NUCLEUS==0   &                   
                                                 HALLMARK_APOPTOSIS==0                   &                                  
                                                 REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_CELL_DEATH_GENES==0   &            
                                                 BIOCARTA_MITOCHONDRIA_PATHWAY==0         &                          
                                                 REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY==0& 
                                                 BIOCARTA_DEATH_PATHWAY==0    
]$gene_symbol



BiocartMitoDeathgenes <-data.table(presabsMatdd)[(BIOCARTA_DEATH_PATHWAY==1&BIOCARTA_MITOCHONDRIA_PATHWAY==1 
                     # &REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY
                            ) &
                           REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES==0 
]$gene_symbol

# More specific cell cycle components
plotthis<-dd0[Gene%in%paste0("CASP",1:10)]
plotthis[   ,LineageLab:="CAMA-1 \n Sensitive"]
plotthis[CellLine=="MCF7"& State=="Sen",LineageLab:="MCF-7 \n Sensitive"]
plotthis[CellLine=="CAMA1"& State=="RiboR",LineageLab:="CAMA-1 \n Resistant"]
plotthis[CellLine=="MCF7"& State=="RiboR",LineageLab:="MCF-7 \n Resistant"]
plotthis$LineageLab<- factor(plotthis$LineageLab, levels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant" ,"MCF-7 \n Sensitive","MCF-7 \n Resistant"  ))
plotthis[,Treatmentlab:="Combination"]
plotthis[Treatment=="Afat",Treatmentlab:="Afatinib"]
plotthis[Treatment=="Ribo",Treatmentlab:="Ribociclib"]
plotthis[Treatment=="DMSO",Treatmentlab:="DMSO"]
plotthis$Treatmentlab <- factor(plotthis$Treatmentlab,  levels=c("DMSO" ,"Ribociclib","Afatinib","Combination"  ))

plotthis$Gene<- factor(plotthis$Gene,
                       levels=c(unique(plotthis$Gene)[unique(plotthis$Gene)%in%paste0("CASP",1:9)],
                                unique(plotthis$Gene)[unique(plotthis$Gene)%in%paste0("CASP",10)]))
plotthis[,DMSOscalelog2cpm:= scalelog2cpm*(Treatment=="DMSO"), by=c("Gene","CellLine","Hour","State")]
plotthis[,Normscalelog2cpm:=scalelog2cpm-DMSOscalelog2cpm, by=c("Gene","CellLine","Hour","State")]
plotthis[,scaleNormSensExpression.log2cpmB:=scale(NormSensExpression.log2cpmB,center =F), by=c("Gene")]
plotthis$Treatmentlab<- factor(plotthis$Treatmentlab,
                       levels=c("DMSO","Ribociclib","Afatinib","Combination" ))
plotthis$Lineage <- factor(plotthis$Lineage  , levels=c("CAMA1 Sen" ,"CAMA1 RiboR"  ,"MCF7 Sen",   "MCF7 RiboR"  ))


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

p2<-plotthis[Treatment!="DMSO"][Gene=="CASP2"]
gam2<- mgcv::gam(scaleNormSensExpression.log2cpmB~-1+
                   s(I((Hour)*(Treatment=="Ribo")),k=3)+
                   s(I((Hour)*(Treatment=="Afat")),k=3)+
                   s(I((Hour)*(Treatment=="Comb")),k=3)
                   ,data=plotthis[Treatment!="DMSO"][Gene=="CASP2"] )
p2$pred<- predict(gam2)
ggplot( plotthis[Treatment!="DMSO"]#[Gene%in%c("CDK1","CDK2","CDK4","CDK6")]#[grep("CDK",)] 
        , aes(y=NormSensExpression.log2cpmB, x=Hour,#shape= Lineage,
              fill=Treatmentlab,col=Treatmentlab, group=interaction(Gene,Lineage,State,CellLine,Treatment)))+
  #geom_path()+
  geom_smooth(method="gam",formula=y~s(x,k=4),aes(group=interaction(Gene,Treatmentlab)) )+
  geom_point(size=1)+
  facet_wrap(~Gene, scales="free_y",nrow=2)+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+
  scale_shape_manual(name="Resistance",values=c(22,23,24,25), labels=c("CAMA-1 Sensitive","MCF-7 Sensitive", "CAMA-1 Resistant","MCF-7 Resistant"))+
  scale_color_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  scale_fill_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  labs(x="Time (Hours)",y="Gene expression \n (log2FC vs DMSO)")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/SynergyComboCASP_ExpressionTimecoursevsDMSO.png", dpi=320, height=6, width=17)

ggplot( plotthis[Treatment!="DMSO"][Day==max(Day)]
        , aes(y=NormSensExpression.log2cpmB, x=Gene,#shape= Lineage,
              fill=Treatmentlab,col=Treatmentlab, group=interaction(Gene,Lineage,State,CellLine,Treatment)))+
  geom_hline(yintercept=0, linetype="dashed")+
  #geom_smooth(method="gam",formula=y~s(x,k=4),aes(group=interaction(Gene,Treatment)) )+
  geom_boxplot(position=position_dodge(width=1),col=1,aes(group=interaction(Gene,Lineage,State,CellLine,Treatment)))+
  geom_point(#aes(shape=Lineage),
             position=position_dodge(width=1),size=1, col=1,pch=21)+
  #facet_wrap(~Gene, scales="free_y",nrow=2)+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+
  scale_shape_manual(name="Resistance",values=c(22,23,24,25), labels=c("CAMA-1 Sensitive","MCF-7 Sensitive", "CAMA-1 Resistant","MCF-7 Resistant"))+
  scale_color_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  scale_fill_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  labs(x="Gene",y="Gene expression \n (log2FC vs DMSO)")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/SynergyComboCASP_ExpressionEndDayvsDMSO.png", dpi=320, height=6, width=17)
ggplot( plotthis[Treatment!="DMSO"][Day==max(Day)]
        , aes(y=NormSensExpression.log2cpmB, x=Gene,#shape= Lineage,
              fill=Treatmentlab,col=Treatmentlab, group=interaction(Gene,Lineage,State,CellLine,Treatment)))+
  geom_hline(yintercept=0, linetype="dashed")+
  #geom_smooth(method="gam",formula=y~s(x,k=4),aes(group=interaction(Gene,Treatment)) )+
  geom_boxplot(position=position_dodge(width=1),col=1,aes(group=interaction(Gene,
                                                                            #Lineage,State,CellLine,
                                                                            Treatmentlab)))+
  geom_point(aes(shape=Lineage),
    position=position_dodge(width=1),size=1, col=1)+#,pch=21)+
  #facet_wrap(~Gene, scales="free_y",nrow=2)+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+
  scale_shape_manual(name="Cell line",values=c(22,23,24,25), labels=c("CAMA-1 Sensitive","MCF-7 Sensitive", "CAMA-1 Resistant","MCF-7 Resistant"))+
  scale_color_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  scale_fill_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  labs(x="Gene",y="Gene expression \n (log2FC vs DMSO)")
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/SynergyComboCASP_ExpressionEndDayvsDMSO_B.png", dpi=320, height=6, width=17)

caspord<-c("CASP9","CASP10" ,"CASP4",  "CASP8"  ,"CASP7",  "CASP3" , "CASP2" , "CASP6" )
plotthisCas<-plotthis[Treatment!="DMSO"][Day==max(Day)]
plotthisCas$Gene<- factor(plotthisCas$Gene, levels=caspord)
ggplot( plotthisCas
        , aes(y=NormSensExpression.log2cpmB, x=Gene,#shape= Lineage,
              fill=Treatmentlab,col=Treatmentlab, group=interaction(Gene,Lineage,State,CellLine,Treatmentlab)))+
  geom_hline(yintercept=0, linetype="dashed")+
  #geom_smooth(method="gam",formula=y~s(x,k=4),aes(group=interaction(Gene,Treatment)) )+
    stat_boxplot(col=1,aes(group=interaction(Gene,
                                     #Lineage,State,CellLine,
                                     Treatmentlab)),position=position_dodge(width=1),geom ='errorbar') + 
  geom_boxplot(outlier.color =NA,position=position_dodge(width=1),col=1,aes(group=interaction(Gene,
                                                                                              #Lineage,State,CellLine,
                                                                                              Treatmentlab)))+
  geom_point(aes(shape=Lineage),
             position=position_dodge(width=1),size=1, col=1)+#,pch=21)+
  #facet_wrap(~Gene, scales="free_y",nrow=2)+
  theme_classic(base_size=26)+theme(aspect.ratio=1)+
  scale_shape_manual(name="Cell line",values=c(22,23,24,25), labels=c("CAMA-1 Sensitive","MCF-7 Sensitive", "CAMA-1 Resistant","MCF-7 Resistant"))+
  scale_color_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  scale_fill_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  labs(x="Gene",y="Gene expression \n (log2FC vs DMSO)")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/SynergyComboCASP_ExpressionEndDayvsDMSO_C.png", dpi=320, height=6, width=17)


# Statistics
finFULLCasp<- rbindlist(mclapply(caspord, function(x){
  pred<-dd0[Gene%in% x][Day ==max(Day)];  cat(x)
  if(length(unique(pred$CellLine))>1){
    m1 <- lmer(NormSensExpression.log2cpmB~ -1 +  RiboTreated + AfatTreated + CombTreated +
                # (0+RiboTreated|Lineage) +
                # (0+AfatTreated|Lineage) +
                # (0+CombTreated|Lineage) +
                 (1|Hour) + (1|Lineage), data=pred)
    outmat <- cbind(Gene=x,Contrast="All_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )
    return(outmat)
  }else{ NULL }
}, mc.cores =detectCores()-1 ))
setnames(finFULLCasp, old=c("rn","Std. Error", "t value", "Pr(>|t|)" ), new=c("Parameter","Std.Error", "tvalue", "pvalue" ))

finFULLCasp[Gene%in%unique(plotthis$Gene)][Parameter=="CombTreated"][order(-abs(pvalue))]
finFULLCasp[Gene%in%unique(plotthis$Gene)][Parameter=="AfatTreated"][order(-Estimate)]
finFULLCasp[Gene%in%unique(plotthis$Gene)][Parameter=="RiboTreated"][order(-Estimate)]


finFULLCasp$Gene<- factor(finFULLCasp$Gene, levels=rev(caspord))

finFULLCasp$Parameter<- factor(finFULLCasp$Parameter, levels=c("RiboTreated","AfatTreated","CombTreated"))

ggplot(finFULLCasp, aes(y=Estimate, x=Gene, fill=Parameter,col=Parameter,
                                         group=interaction(Gene,Parameter)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  #facet_grid(~DayLab)+
  #geom_point(position=position_dodge(width=1),size=5)+
  coord_flip()+
  theme_classic(base_size=26)+
  theme(aspect.ratio=2.5)+
  theme(aspect.ratio=1)+
  geom_bar(#aes(alpha=pvalue<0.05),
           position="dodge", stat="identity")+
  labs(y="Effect on pathway activity", x="Pathway")+
  scale_color_manual(name="Treatment \n effect", 
                     values=gg_color_hue(4)[-1],
                     labels=c("Ribociclib","Afatinib", "Synergy"))+
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))


#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/EffectSizesSynergyComboCASP_ExpressionEndDayvsDMSO_C.png", dpi=320, height=6, width=17)






# Statistics
bclord <- c(
            rev(c("BCL2A1","BCL2L10","BCL2L1","BCL2L2","BCL2")),
            rev(c( "BIM","PUMA", "NOXA","BMF","BIK",
            "BAK1","BAX","BOK"      ,"HRK","BAD","BID"  ))            )


finFULLBCLfam<- rbindlist(mclapply(bclord, function(x){
  pred<-dd0[Gene%in% x][Day ==max(Day)];  cat(x)
  if(length(unique(pred$CellLine))>1){
    m1 <- lmer(NormSensExpression.log2cpmB~ -1 +  RiboTreated + AfatTreated + CombTreated +
                 # (0+RiboTreated|Lineage) +
                 # (0+AfatTreated|Lineage) +
                 # (0+CombTreated|Lineage) +
                 (1|Hour) + (1|Lineage), data=pred)
    outmat <- cbind(Gene=x,Contrast="All_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )
    return(outmat)
  }else{ NULL }
}, mc.cores =detectCores()-1 ))
setnames(finFULLBCLfam, old=c("rn","Std. Error", "t value", "Pr(>|t|)" ), new=c("Parameter","Std.Error", "tvalue", "pvalue" ))

finFULLBCLfam[][Parameter=="CombTreated"][order(-abs(pvalue))]
finFULLBCLfam[][Parameter=="AfatTreated"][order(-Estimate)]
finFULLBCLfam[][Parameter=="RiboTreated"][order(-Estimate)]


finFULLBCLfam$Gene<- factor(finFULLBCLfam$Gene, levels=rev(bclord))

finFULLBCLfam$Parameter<- factor(finFULLBCLfam$Parameter, levels=c("RiboTreated","AfatTreated","CombTreated"))

ggplot(finFULLBCLfam, aes(y=Estimate, x=Gene, fill=Parameter,col=Parameter,
                        group=interaction(Gene,Parameter)))+
  geom_hline(yintercept = 0, linetype="dashed")+
  #facet_grid(~DayLab)+
  #geom_point(position=position_dodge(width=1),size=5)+
  coord_flip()+
  theme_classic(base_size=26)+
  theme(aspect.ratio=2.5)+
  theme(aspect.ratio=1)+
  geom_bar(aes(alpha=pvalue<0.05),
    position="dodge", stat="identity")+
  labs(y="Effect on pathway activity", x="Pathway")+
  scale_color_manual(name="Treatment \n effect", 
                     values=gg_color_hue(4)[-1],
                     labels=c("Ribociclib","Afatinib", "Synergy"))+
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))


#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/EffectSizesSynergyComboCASP_ExpressionEndDayvsDMSO_C.png", dpi=320, height=6, width=17)





ggplot( plotthis[Treatment!="DMSO"][Day==max(Day)]
        , aes(y=NormSensExpression.log2cpmB, x=Gene,#shape= Lineage,
              fill=Treatmentlab,col=Treatmentlab, group=interaction(Gene,Lineage,State,CellLine,Treatmentlab)))+
  geom_hline(yintercept=0, linetype="dashed")+
  #geom_smooth(method="gam",formula=y~s(x,k=4),aes(group=interaction(Gene,Treatment)) )+
  geom_boxplot(position=position_dodge(width=1),col=1,aes(group=interaction(Gene,#Lineage,State,CellLine,
                                                                            Treatmentlab)))+
  geom_point(aes(shape=Lineage),
    position=position_dodge(width=1),size=1, col=1)+#,pch=21)+
  facet_grid(~Treatmentlab,scales="free_x")+
  theme_classic(base_size=12)+theme(aspect.ratio=1)+
  scale_shape_manual(name="Resistance",values=c(22,23,24,25), labels=c("CAMA-1 Sensitive","MCF-7 Sensitive", "CAMA-1 Resistant","MCF-7 Resistant"))+
  scale_color_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  scale_fill_manual(name="Treatment", values=gg_color_hue(4)[-1],labels=c( "Ribociclib","Afatinib",  "Combination"))+
  labs(x="Gene",y="Gene expression \n (log2FC vs DMSO)")




ggplot(plotthis[Treatment!="DMSO"],#[Gene_Set=="HALLMARK_APOPTOSIS"],
       aes(x=as.factor(Hour), y=Gene, fill=scaleNormSensExpression.log2cpmB ))+
  theme_classic(base_size=16)+
  geom_tile()+
  facet_grid(Treatmentlab~LineageLab,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  scale_fill_gradient2(name="Gene set enrichment \n scaled(compared to DMSO)" , low="blue", high="darkred", mid="white")+
  theme(aspect.ratio=1)+#labs(title=gsub("_"," ",x), x="Hour")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")+
  scale_y_discrete(guide=guide_axis(n.dodge=1)) +
  labs(x="Time (hours)",y="Pathway")



                

ssGSEAmeta[gs_name %in%c("BIOCARTA_DEATH_PATHWAY","BIOCARTA_MITOCHONDRIA_PATHWAY")]
ssGSEAmeta[gs_name %in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")]
apgenes <- (ssGSEAmeta[gs_name %in%c("BIOCARTA_DEATH_PATHWAY","BIOCARTA_MITOCHONDRIA_PATHWAY"
                                     )][]$gene_symbol)%>%unique()%>%sort()
#ccgenes <- (ssGSEAmeta[gs_name %in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT"
#                                     )]$gene_symbol)%>%unique()%>%sort()
ccgenes <-data.table(unique(ssGSEAmeta[gs_name %in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT",
                                           "HALLMARK_MITOTIC_SPINDLE","BIOCARTA_RACCYCD_PATHWAY",
                                           "BIOCARTA_MCM_PATHWAY",
                                           "BIOCARTA_CELLCYCLE_PATHWAY",
                                           "BIOCARTA_G1_PATHWAY","BIOCARTA_G2_PATHWAY"
)]%>%select(gs_name,gene_symbol))%>%group_by(gene_symbol)%>%summarise(n=n()))$gene_symbol%>%sort()

apgenes <- (ssGSEAmeta[gs_name %in%c("BIOCARTA_DEATH_PATHWAY","BIOCARTA_MITOCHONDRIA_PATHWAY"#,
                                     #"HALLMARK_APOPTOSIS"#, 
                                     #"HALLMARK_P53_PATHWAY"
)][]$gene_symbol)%>%unique()%>%sort()


subcc<-finFULL[Gene%in%c(ccgenes)]
setnames(subcc,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
subap<-finFULL[Gene%in%c(apgenes)]
setnames(subap,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))

ggplot(subcc[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][pvalue<0.05],aes(y=Estimate,x=rn))+geom_boxplot()+geom_point()+
  facet_wrap(~Contrast)+coord_flip() #+geom_violin()
ggplot(subap[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][pvalue<0.05],aes(y=Estimate,x=rn))+geom_boxplot()+
  geom_point()+
  facet_wrap(~Contrast)+coord_flip() 


finFULLEffs<-data.table(finFULL[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%apgenes]%>%group_by(Gene)%>%
                          mutate(                                 maxeff=max(Estimate),
                                 mineff=max(Estimate),
                                 medeff=median(Estimate),
                                 meaneff=mean(Estimate),
                                 meansynergeff=mean((rn=="CombTreated")*Estimate)
                                 ))[
                                   order(meaneff)]

finFULLEffs$Gene<- factor(finFULLEffs$Gene,levels=unique(finFULLEffs$Gene))
finFULLEffs[, CellLine :="Resistant"]
finFULLEffs[Contrast=="Sen_Comb_vs_DMSO", CellLine :="Sensitive"]
finFULLEffs$CellLine<- factor(finFULLEffs$CellLine,levels=c("Sensitive","Resistant"))
finFULLEffs$rn<- factor(finFULLEffs$rn,levels=c("RiboTreated","AfatTreated","CombTreated"))
setnames(finFULLEffs,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))


ggplot(finFULLEffs[c(1:(10*6) ,(nrow(finFULLEffs)-(c(1:(10*6) )-1)) )
],aes(y=Estimate,x=Gene,fill=rn))+  theme_classic()+
  #geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~CellLine)+coord_flip() +
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Top10CellDeathGeneSpecificMainTreatmentAndSynergyByEnd.pdf", dpi=320, height=6, width=6)

ggplot(finFULLEffs[  ][],aes(y=Estimate,x=Gene,fill=rn#,alpha=pvalue<0.05
                             ))+  theme_classic()+
     #geom_boxplot()+
     geom_bar(position="dodge", stat="identity")+
     facet_wrap(~CellLine)+coord_flip() +
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificMainTreatmentAndSynergyByEnd.pdf", dpi=320, height=6, width=6)

ggplot(finFULLEffs[  ][],aes(y=Estimate,x=Gene,fill=rn#,alpha=pvalue<0.05
))+  theme_classic()+
  #geom_boxplot()+
  geom_bar(aes(alpha=pvalue<0.05),position="dodge", stat="identity")+
  facet_wrap(~CellLine)+coord_flip() +
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))

apgenetop<-unique(finFULLEffs[c(1:(10*6) ,(nrow(finFULLEffs)-(c(1:(10*6) )-1)) )]$Gene)
finFULLEffs[pvalue<0.05]

plotdd0<-dd0[Gene%in% apgenetop][Day ==max(Day)][Treatment!="DMSO"]
plotdd0$Gene<- factor(plotdd0$Gene,levels=unique(finFULLEffs$Gene))
ggplot(plotdd0,aes(y=NormSensExpression.log2cpmB,x=Gene,
                                                 col=Treatment))+  theme_classic()+
  geom_boxplot()+
  geom_point(aes(group=Treatment),position=position_dodge(width=1), stat="identity")+
  facet_wrap(~State)+coord_flip() 



ccfinFULLEffs<-data.table(finFULL[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%ccgenes]%>%group_by(Gene)%>%
                          mutate(maxeff=max(Estimate),
                                 medeff=median(Estimate),
                                 meaneff=mean(Estimate),
                                 meansynergeff=mean((rn=="CombTreated")*Estimate)
                          ))[
                            order(meaneff)]
setnames(ccfinFULLEffs,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
ccfinFULLEffs[,minpval:=min(pvalue),by="Gene"]
ccfinFULLEffs$Gene<- factor(ccfinFULLEffs$Gene,levels=rev(unique(ccfinFULLEffs$Gene)))
ccfinFULLEffs[, CellLine :="Resistant"]
ccfinFULLEffs[Contrast=="Sen_Comb_vs_DMSO", CellLine :="Sensitive"]
ccfinFULLEffs$CellLine<- factor(ccfinFULLEffs$CellLine,levels=c("Sensitive","Resistant"))
ccfinFULLEffs$rn<- factor(ccfinFULLEffs$rn,levels=c("RiboTreated","AfatTreated","CombTreated"))
setnames(ccfinFULLEffs,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))

ccfinFULLEffs
ggplot(ccfinFULLEffs[c(1:(50*6) ,(nrow(ccfinFULLEffs)-(c(1:(50*6) )-1)) )
],aes(y=Estimate,x=Gene,fill=rn, alpha=pvalue<0.05))+  theme_classic()+
  #geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~CellLine)+coord_flip() +
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))

ggplot(ccfinFULLEffs[ minpval<0.05 ][],aes(y=Estimate,x=Gene,fill=rn,alpha=pvalue<0.05))+
  theme_classic(base_size=5)+
  #geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~CellLine)+coord_flip() +
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))


ccgenetop<-unique(ccfinFULLEffs[c(1:(15*6) ,(nrow(ccfinFULLEffs)-(c(1:(15*6) )-1)) )]$Gene)

ccplotdd0<-dd0[Gene%in% ccgenetop][Day ==max(Day)][Treatment!="DMSO"]
ccplotdd0$Gene<- factor(ccplotdd0$Gene,levels=unique(ccfinFULLEffs$Gene))
ggplot(ccplotdd0,aes(y=NormSensExpression.log2cpmB,x=Gene,
                   col=Treatment))+  theme_classic()+
  geom_boxplot()+
  geom_point(aes(group=Treatment),position=position_dodge(width=1), stat="identity")+
  facet_wrap(~State)+coord_flip() 




ggplot(finFULL[Contrast!="DeltaState_Comb_vs_DMSO"][grepl("CASP",Gene)],aes(y=Estimate,x=Gene,fill=rn))+
  theme_classic()+
  geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~Contrast)+coord_flip() 



apstatsR<-fincomp[Contrast=="RiboR_Comb_vs_DMSO"][Gene%in%apgenes]
apstatsR[,apFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
apfincompR <- apstatsR[pvalue<0.05][order(-(Estimate))]
apfincompR[, EstID:=1:nrow(apfincompR)]
apfincompR[, EstDir:=-1+2*(Estimate>0)]

apstatsS<-fincomp[Contrast=="Sen_Comb_vs_DMSO"][Gene%in%apgenes]
apstatsS[,apFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
apfincompS <- apstatsS[pvalue<0.05][order(-(Estimate))]
apfincompS[, EstID:=1:nrow(apfincompS)]
apfincompS[, EstDir:=-1+2*(Estimate>0)]


ccstatsR<-fincomp[Contrast=="RiboR_Comb_vs_DMSO"][Gene%in%ccgenes]
ccstatsR[,ccFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
ccfincompR <- ccstatsR[pvalue<0.05][order(-(Estimate))]
ccfincompR <- ccstatsR[ccFDR<0.05][order(-(Estimate))]
ccfincompR[, EstID:=1:nrow(ccfincompR)]
ccfincompR[, EstDir:=-1+2*(Estimate>0)]

ccstatsS<-fincomp[Contrast=="Sen_Comb_vs_DMSO"][Gene%in%ccgenes]
ccstatsS[,ccFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
ccfincompS <- ccstatsS[pvalue<0.05][order(-(Estimate))]
ccfincompS <- ccstatsS[ccFDR<0.05][order(-(Estimate))]
ccfincompS[, EstID:=1:nrow(ccfincompS)]
ccfincompS[, EstDir:=-1+2*(Estimate>0)]





ccRplotGenes <- ccfincompR$Gene
ccSplotGenes <- ccfincompS$Gene

ccimpacted<-unique(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%c(ccRplotGenes,ccSplotGenes)]$gene_symbol)
ccimpactedEffects<-fincomp[Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%ccimpacted][order(-Estimate)]
ccimpactedEffects$Gene<- factor(ccimpactedEffects$Gene, levels=unique(ccimpactedEffects$Gene))
ccmuimpactedEffects <-data.table(ccimpactedEffects%>%group_by(Gene)%>%dplyr::summarise(Estimate=mean(Estimate)))[order(-Estimate)]
ccmuimpactedEffects$Gene<- factor(ccmuimpactedEffects$Gene, levels=unique(ccmuimpactedEffects$Gene))
ccmuimpactedEffects[,rankabsEst:=rank(-abs(Estimate))]
ggplot( ccimpactedEffects, 
        aes(y= Gene, x= Estimate,col=Contrast))+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_point()+ labs(x="Synergistic effect")+
  scale_color_npg(name="Resistance", labels=c("Resistant","Sensitive"))

ggplot( ccmuimpactedEffects, 
        aes(y= Gene, x= Estimate))+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_point()+ labs(x="Synergistic effect")
ggplot( ccmuimpactedEffects[rankabsEst<=20], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_void(base_size=14)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+
  geom_text(aes(label=Gene ))+
  theme(axis.ticks=element_blank(),axis.text.y=element_blank())+
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEnd.pdf", dpi=320, height=6, width=6)


#labs(x="Synergistic effect")
ggplot( ccmuimpactedEffects[rankabsEst<=40], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)
ggplot( ccmuimpactedEffects[c(1:20, (nrow(ccmuimpactedEffects)-19):nrow(ccmuimpactedEffects))], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)

#lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[c(1:20, (nrow(muimpactedEffects)-19):nrow(muimpactedEffects))]$Gene ))
#lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[c(1:10, (nrow(muimpactedEffects)-19):nrow(muimpactedEffects))]$Gene ))
#lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[]$Gene ))
cclutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=ccmuimpactedEffects[c(1:2, (nrow(ccmuimpactedEffects)-37):nrow(ccmuimpactedEffects))]$Gene ))

cclutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=ccmuimpactedEffects[rankabsEst<=40]$Gene ))
cclutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=ccmuimpactedEffects[rankabsEst<=30]$Gene ))
cclutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=ccmuimpactedEffects[]$Gene ))

ccoutsyn<- rbindlist(lapply(1:nrow(cclutab), function(i){
  tmp<- dd0[Gene%in% cclutab[i]$Gene][Day ==max(Day)][Lineage==cclutab[i]$Lineage][timeunderdrug%in%c(0,6,24)]
  m1 <- lmer(NormSensExpression.log2cpm~ -1 +  RiboTreated + AfatTreated + CombTreated +
               (1|Hour) , data=tmp)
  cbind(Gene=cclutab[i]$Gene,Lineage=cclutab[i]$Lineage,Contrast="DeltaState_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )[rn=="CombTreated"]
}))

ccoutsyn$Lineage <- factor(ccoutsyn$Lineage, levels=c("CAMA1 Sen" , "CAMA1 RiboR" , "MCF7 Sen" , "MCF7 RiboR") )
ggplot( ccoutsyn,
        aes(y= Gene, x= Lineage, fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=3)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSI.pdf", dpi=320, height=12, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSI.png", dpi=320, height=12, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSI.tiff", dpi=320, height=12, width=6)

corecc<-c("TOP2A","BIRC5","CDC25C","CDK1","CCNB1","CCNB2","AURKA","MYC","CDC25A","CDK2","CCNE1","CDK4","CCND1","E2F1","CDKN1A"
          #"ORC1","STMN1","MCM7","TPX2","PLK1",
          )
ggplot( ccoutsyn[Gene%in%corecc],
        aes(y= Gene, x= Lineage, fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergy \n effect",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.pdf", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.png", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.tiff", dpi=320, height=6, width=6)



RplotGenes <- apfincompR$Gene
SplotGenes <- apfincompS$Gene

impacted<-unique(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%c(RplotGenes,SplotGenes)]$gene_symbol)
impacted<-unique(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%apgenes]$gene_symbol)
impactedEffects<-fincomp[Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%impacted][order(-Estimate)]
impactedEffects$Gene<- factor(impactedEffects$Gene, levels=unique(impactedEffects$Gene))
muimpactedEffects <-data.table(impactedEffects%>%group_by(Gene)%>%dplyr::summarise(Estimate=mean(Estimate)))[order(-Estimate)]
muimpactedEffects$Gene<- factor(muimpactedEffects$Gene, levels=rev(unique(muimpactedEffects$Gene)))
muimpactedEffects[,rankabsEst:=rank(-abs(Estimate))]
ggplot( impactedEffects, 
        aes(y= Gene, x= Estimate,col=Contrast))+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_point()+ labs(x="Synergistic effect")+
  scale_color_npg(name="Resistance", labels=c("Resistant","Sensitive"))

ggplot( muimpactedEffects, 
        aes(y= Gene, x= Estimate))+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_point()+ labs(x="Synergistic effect")
ggplot( muimpactedEffects[rankabsEst<=20], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_void(base_size=14)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+
  geom_text(aes(label=Gene ))+
  theme(axis.ticks=element_blank(),axis.text.y=element_blank())+
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEnd.pdf", dpi=320, height=6, width=6)


#labs(x="Synergistic effect")
ggplot( muimpactedEffects[rankabsEst<=30], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)

#lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[c(1:20, (nrow(muimpactedEffects)-19):nrow(muimpactedEffects))]$Gene ))
#lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[c(1:10, (nrow(muimpactedEffects)-19):nrow(muimpactedEffects))]$Gene ))
lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[]$Gene ))
#lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[c(1:2, (nrow(muimpactedEffects)-27):nrow(muimpactedEffects))]$Gene ))
#lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[rankabsEst<=30]$Gene ))
outsyn<- rbindlist(lapply(1:nrow(lutab), function(i){
  tmp<- dd0[Gene%in% lutab[i]$Gene][Day ==max(Day)][Lineage==lutab[i]$Lineage][timeunderdrug%in%c(0,6,24)]
  m1 <- lmer(NormSensExpression.log2cpm~ -1 +  RiboTreated + AfatTreated + CombTreated +
               (1|Hour) , data=tmp)
  cbind(Gene=lutab[i]$Gene,Lineage=lutab[i]$Lineage,Contrast="DeltaState_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )[rn=="CombTreated"]
}))

outsyn$Lineage <- factor(outsyn$Lineage, levels=c("CAMA1 Sen" , "CAMA1 RiboR" , "MCF7 Sen" , "MCF7 RiboR") )
ggplot( outsyn,
        aes(y= Gene, x= Lineage, fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=3)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSI.pdf", dpi=320, height=12, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSI.png", dpi=320, height=12, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSI.tiff", dpi=320, height=12, width=6)

coreap<-c("CASP3","CASP6","CASP7","CASP8", "CASP9" ,"CASP10","BIRC3","TNFSF10","MAP3K14","TRADD","BIK","BAK1","BAX","BID","BCL2")
ggplot( outsyn[Gene%in%coreap],
        aes(y= Gene, x= Lineage, fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergy \n effect",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.pdf", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.png", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.tiff", dpi=320, height=6, width=6)







finFULLEffs<-data.table(finFULL[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%c(corecc,coreap)]%>%group_by(Gene)%>%
                          mutate(                                 maxeff=max(Estimate),
                                                                  mineff=max(Estimate),
                                                                  medeff=median(Estimate),
                                                                  meaneff=mean(Estimate),
                                                                  meansynergeff=mean((rn=="CombTreated")*Estimate)
                          ))[
                            order(meaneff)]



finFULLEffs$Gene<- factor(finFULLEffs$Gene,levels=unique(finFULLEffs$Gene))
finFULLEffs$Gene<- factor(finFULLEffs$Gene,levels=rev(c(corecc,coreap)))
finFULLEffs[, CellLine :="Resistant"]
finFULLEffs[Contrast=="Sen_Comb_vs_DMSO", CellLine :="Sensitive"]
finFULLEffs$CellLine<- factor(finFULLEffs$CellLine,levels=c("Sensitive","Resistant"))
finFULLEffs$rn<- factor(finFULLEffs$rn,levels=c("RiboTreated","AfatTreated","CombTreated"))
setnames(finFULLEffs,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
finFULLEffs[,Theme:="Cell Death"]
finFULLEffs[Gene%in%corecc,Theme:="Cell Cycle"]
finFULLEffs[,Effect:=(Estimate )]
finFULLEffs[,Effect2:=sum(Estimate ), by=c("Gene","CellLine")]
finFULLEffs[rn=="CombTreated",Effect:=(Effect2 )]


ggplot(finFULLEffs[#c(1:(10*6) ,(nrow(finFULLEffs)-(c(1:(10*6) )-1)) )
],aes(y=Effect,x=Gene,fill=rn))+  theme_classic(base_size=12)+
  #geom_boxplot()+
  #geom_bar(position="dodge", stat="identity")+
  geom_col(aes(col=rn),width=0.25,    
           position=position_dodge(0.5))+
 # geom_point(size=3,aes(shape=Theme))+
  facet_grid(Theme~CellLine,scale="free_y")+coord_flip() +
  scale_fill_manual(name="Treatment", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Combination"))+
  scale_color_manual(name="Treatment", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Combination"))+
  labs(y="Treatment impact")+
  theme(aspect.ratio=1)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Top15CellCYCLEandCellDeathGeneSpecificMainTreatmentAndSynergyByEnd.pdf", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Top15CellCYCLEandCellDeathGeneSpecificMainTreatmentAndSynergyByEnd.png", dpi=320, height=6, width=6)



ggplot(finFULLEffs[  ][],aes(y=Estimate,x=Gene,fill=rn#,alpha=pvalue<0.05
))+  theme_classic()+
  #geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~CellLine)+coord_flip() +
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificMainTreatmentAndSynergyByEnd.pdf", dpi=320, height=6, width=6)











ggplot(subap[Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%outsyn$Gene],aes(y=Estimate,x=Gene,fill=rn))+geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~Contrast)+coord_flip() 
ggplot(subcc[Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%ccoutsyn$Gene],aes(y=Estimate,x=Gene,fill=rn))+geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~Contrast)+coord_flip() 

