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

# Specify source data location 
fileloc <- ("~/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse")
setwd(fileloc)

# load timecourse rnaseq data
# If first time using code, translate .xlsx file to .csv format for ease of reading 
#dd0mcf <-read.xlsx(xlsxFile = "sh10050_MCF7_genes_gam_short.term_inputData_synergy.xlsx")
#write.csv(dd0mcf,file= "sh10050_MCF7_genes_gam_short.term_inputData_synergy.csv")
#dd0cama<-read.xlsx(xlsxFile = "sh11141_CAMA1_genes_gam_short.term_inputData_synergy.xlsx")
#write.csv(dd0cama,file= "sh11141_CAMA1_genes_gam_short.term_inputData_synergy.csv")
dd0mcf <- read.csv(file = "sh10050_MCF7_genes_gam_short.term_inputData_synergy.csv")
dd0cama <- read.csv(file = "sh11141_CAMA1_genes_gam_short.term_inputData_synergy.csv")

# Join data across cell lines
dd0 <- data.table(rbind(dd0cama%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ),
                        dd0mcf%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )))
# order factors
dd0$Treatment <- factor(dd0$Treatment,c("DMSO", "Afat","Ribo", "Comb"))
dd0$State <- factor(dd0$State,c("Sen","RiboR"))
dd0$CellLine <- factor(dd0$CellLine,c("CAMA1", "MCF7"))

# Generate day number column
dd0[, Day:=(Hour-as.numeric(as.character(timeunderdrug)))/24 ]

# Scale each genes expression data for each cell line (line & resistance)
dd0[,scalelog2cpm:= scale(Expression.log2cpm), by=c("Gene","CellLine","State")]
dd0[,scalelog2cpmB:= scale(Expression.log2cpm), by=c("Gene","CellLine")]
dd0[,scalelog2cpmC:= scale(Expression.log2cpm), by=c("Gene")]

# Scale cell line expression of gene g relative to sensitive cells at time t under DMSO
dd0[,SensExpression.log2cpm := sum(Expression.log2cpm*(State=="Sen"&Treatment=="DMSO")), by=c("Gene","CellLine","Hour")]
dd0[,SensExpression.log2cpmB:= sum(Expression.log2cpm*(Treatment=="DMSO" )), by=c("Gene","CellLine","Hour","State")]
dd0[,NormSensExpression.log2cpm:=Expression.log2cpm -SensExpression.log2cpm, by=c("Gene","CellLine","Hour", "Treatment")]
dd0[,NormSensExpression.log2cpmB:=Expression.log2cpm -SensExpression.log2cpmB, by=c("Gene","CellLine","Hour", "Treatment")]

# Generate and order the cell lineage column
dd0[,Lineage:=paste0(CellLine," ",State)]
dd0$Lineage <- factor(dd0$Lineage  , levels=c("CAMA1 Sen" ,"MCF7 Sen", "CAMA1 RiboR"  ,  "MCF7 RiboR"  ))

# list all genes
selectgenes<- unique(dd0$Gene)

#load ssGSEA pathway gene set metadata
ssGSEAmeta <- data.table( msigdbr(species = "Homo sapiens") )
ssGSEAmeta[gs_cat%in%c("C2","H")][gs_subcat%in%c("","CP:BIOCARTA","CP:REACTOME","CP:KEGG" )]$gs_name%>%unique()
ssGSEAmeta[gs_cat%in%c("C2","H")][gs_subcat%in%c("CP:BIOCARTA")]$gs_name%>%unique()
ssGSEAmeta[,ngenes:=length(unique(gene_symbol)), by=gs_name]

# Compare pre-treatment gene expression of resistant and sensitive cells
dd0$Gene%>%unique()%>%length()
strtFULL<- rbindlist(mclapply(selectgenes, function(x){
  pred<-dd0[Gene%in% x][Day ==min(Day)][Hour==0][Treatment=="DMSO"]
  #ggplot(pred,  aes(x=Treatment, y=NormSensExpression.log2cpm,col=State))+theme_classic()+ geom_point()+facet_grid(~CellLine,scales="free")+geom_hline(linetype=3,yintercept = 0)
  cat(x)
  if(length(unique(pred$CellLine))>1){
    m1 <- lm(NormSensExpression.log2cpm~ -1 + State , data=pred)
    #summary(m1)
    outp1 <- cbind(Gene=x,Contrast="DeltaState_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )
    #outmat <- rbind(outp1[rn=="StateRiboR:CombTreated"],outpS[rn=="CombTreated"],outpR[rn=="CombTreated"])
    outmat <- rbind(outp1 )
    return(outmat)
  }else{ NULL }
} , mc.cores =detectCores()-1 ))
strtcomp <- strtFULL[rn=="StateRiboR"]#[rn%in%c("StateRiboR:CombTreated","CombTreated")]
setnames(strtcomp,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))

# Perform FDR correction and identify effect direction
strtcomp1 <- data.table(strtcomp)
strtcomp1[,FDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
strtcomp1[, EstID:=1:nrow(strtcomp1)]
strtcomp1[, EstDir:=-1+2*(Estimate>0)]

##List genes in focal myc and cell-ceycle-related gene sets
##ssGSEAmeta [gs_cat%in%c("H","C2")][gene_symbol=="PSMC3IP"][gs_subcat%in%c("CP:REACTOME","CP:BIOCARTA")]
#ccgenes <- unique(ssGSEAmeta[gs_name%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V2")]$gene_symbol)
#ccgenes<-unique(ssGSEAmeta[gs_name%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V2","HALLMARK_MYC_TARGETS_V1",
#                                        "BIOCARTA_CELLCYCLE_PATHWAY","BIOCARTA_G2M_PATHWAY")]$gene_symbol)
## extract gene expression data for genes in focal cell-ceycle-related gene sets
#cceffects <- strtcomp[Gene%in% ccgenes][order(-(Estimate))]
##write.csv(cceffects[pvalue<0.05], file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/tables/PretxResistantCellcycleGeneDE_HALLMARKorBIOCARTAcellcyclemyc.csv")

#List genes in focal cell-ceycle-related gene sets
ccgenes<-unique(ssGSEAmeta[gs_name%in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","BIOCARTA_G1_PATHWAY","BIOCARTA_G2M_PATHWAY")]$gene_symbol)
cceffects <- strtcomp[Gene%in% ccgenes][order(-(Estimate))]
#write.csv(cceffects[pvalue<0.05], file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/tables/PretxResistantCellcycleGeneDE_HALLMARKorBIOCARTA.csv")
ccgenes%>%length() 
# 341 cell cycle related genes (from 14860 total genes)

# look at output for specific genes
#strtcomp[Gene%in% c("CCNB1","CCNB2","CDK1","CDK2")][order(-(Estimate))]
#cceffects[pvalue<0.05][Gene%in% c("CCNB1","CCNB2","CDK1","CDK2")][order(-(Estimate))]

# Identify significantly DE genes between resistant and sensitive cells and rank by estimate
muimpactedEffects <- cceffects[pvalue<0.05][order(-(Estimate))]#[abs(Estimate)>log2(2)]
muimpactedEffects[,rankabsEst:=rank(-abs(Estimate))]
muimpactedEffects$Gene<- factor(muimpactedEffects$Gene, levels=rev(unique(muimpactedEffects$Gene)))

# 
ggplot( muimpactedEffects[rankabsEst<=15], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_void(base_size=14)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+
  geom_text(aes(label=Gene ))+
  theme(axis.ticks=element_blank(),axis.text.y=element_blank())+
  scale_fill_gradient2(name="Resistant expression \n (log2FC vs Sensitive)",low="blue", high="darkred", mid="white",midpoint=0)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/GeneSpecificDiffRvsSpreTx.pdf", dpi=320, height=6, width=6)



#Supervised analysis 1 of core cell cycle genes
tmp <- cceffects[Gene%in%c("CDKN1A","CDKN1B",
                         "CDC25A","CDK2","CDK1","CCNE1")]#,"CCNB1","CCNB2"
tmp$Gene<- factor(tmp$Gene, levels=rev(tmp$Gene))
ggplot( tmp, aes(y= Gene, x= "", fill=Estimate))+
  theme_void(base_size=14)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+
  geom_text(aes(label=Gene ))+
  theme(axis.ticks=element_blank(),axis.text.y=element_blank())+
  scale_fill_gradient2(name="Resistant expression \n (log2FC vs Sensitive)",low="blue", high="darkred", mid="white",midpoint=0)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CoreGeneSpecificDiffRvsSpreTx.pdf", dpi=320, height=6, width=6)

#Supervised analysis 2 of core cell cycle genes
coreset<-c("CCNA1","CCNA2", "CCNB1", "CCNB2", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2",
           "CDKN2B", "CDKN2A", "CDKN2C", "CDKN2D", "CDKN1A", "CDKN1B",
           "CDK1", "CDK2", "CDK4", "CDK6", "RB1",
           "E2F1", "E2F2", "E2F3", "CDC25A", "CDC25B", "CDC25C")
coreset<-c("CDC25A","CCNE1","MELK","TK1","RRM2","CKS1B","TUBB","CDC25B","CDK2","MCM4","CDK4",
           "CDC27","E2F3","E2F4","ORC2")
coremuimpactedEffects <- muimpactedEffects[Gene%in% coreset][order(-(Estimate))]
coremuimpactedEffects[,rankabsEst:=rank(-abs(Estimate))]
coremuimpactedEffects$Gene<- factor(coremuimpactedEffects$Gene, levels=rev(unique(coremuimpactedEffects$Gene)))
ggplot( coremuimpactedEffects, 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_void(base_size=14)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+
  geom_text(aes(label=Gene ))+
  theme(axis.ticks=element_blank(),axis.text.y=element_blank())+
  scale_fill_gradient2(name="Resistant expression \n (log2FC vs Sensitive)",low="blue", high="darkred", mid="white",midpoint=0)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CoreGeneSpecificDiffRvsSpreTx.pdf", dpi=320, height=6, width=6)

ggplot( coremuimpactedEffects[pvalue<0.05], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_void(base_size=14)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+
  geom_text(aes(label=Gene ))+
  theme(axis.ticks=element_blank(),axis.text.y=element_blank())+
  scale_fill_gradient2(name="Resistant expression \n (log2FC vs Sensitive)",low="blue", high="darkred", mid="white",midpoint=0)


# Extract gene expression of top DE genes pre-treatment
actualdata<-dd0[Gene%in%muimpactedEffects[rankabsEst<=15]$Gene][Day ==min(Day)][Hour==0][Treatment=="DMSO"]
actualdata<-dd0[Gene%in%coreset][Day ==min(Day)][Hour==0][Treatment=="DMSO"]
actualdata$Gene<-factor(actualdata$Gene, levels=rev(muimpactedEffects$Gene))
actualdata$Lineage<-factor(actualdata$Lineage, levels=c("CAMA1 Sen", "CAMA1 RiboR", "MCF7 Sen",  "MCF7 RiboR"))

# save output
resout<-data.table(actualdata%>%dplyr::select(Gene,Lineage,CellLine,State,Expression.cpm,scalelog2cpmC))[order(Gene,State,CellLine)]
#write.csv(resout, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/tables/PretxResistantCellcycleGeneExpression.csv")

# Visualize
ggplot(actualdata, 
       aes(y= Gene, x= Lineage, fill=scalelog2cpmC))+
  theme_classic(base_size=14)+
  theme(aspect.ratio=1.5)+labs(x="")+
  scale_x_discrete(labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive", "MCF-7 \n Resistant"))+
  geom_tile()+
  scale_fill_gradient2(name="Expression \nscaled(log2(CPM))",low="blue", high="darkred", mid="white",midpoint=0)+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Pre_tx_ExpressionGeneDE_1_2.png", dpi=320, height=6, width=17)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Pre_tx_ExpressionGeneDE_KKselected.png", dpi=320, height=6, width=8)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Pre_tx_ExpressionGeneDE_KKselected.pdf", dpi=320, height=6, width=8)


ggplot(actualdata, 
       aes(y= Gene, x= Lineage, fill=scalelog2cpmC))+
  theme_classic(base_size=14)+
  theme(aspect.ratio=1.5)+labs(x="")+
  scale_x_discrete(labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive", "MCF-7 \n Resistant"))+
  geom_tile()+
  scale_fill_gradient2(name="Expression \nscaled(log2(CPM))",low="blue", high="darkred", mid="white",midpoint=0)+
  theme(text=element_blank(), legend.text =element_blank() , legend.title =element_blank() )
  
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Pre_tx_ExpressionGeneDE_1_2.png", dpi=320, height=6, width=17)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/BLANK_Pre_tx_ExpressionGeneDE_KKselected.png", dpi=320, height=6, width=8)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/BLANK_Pre_tx_ExpressionGeneDE_KKselected.pdf", dpi=320, height=6, width=8)



ggplot( actualdata,
        aes(y= Gene, x= Lineage, fill=scalelog2cpmC))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Expression \n scaled(log2(x))",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Pre_tx_ExpressionGeneDE_1_2.pdf", dpi=320, height=6, width=17)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Pre_tx_ExpressionGeneDE_1_2.png", dpi=320, height=6, width=17)

ggplot( actualdata,
        aes(y= Gene, x= Lineage, fill=scalelog2cpmC))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Expression \n scaled(log2(x))",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))


#
coreactualdata<-dd0[Gene%in%coremuimpactedEffects$Gene][Day ==min(Day)][Hour==0][Treatment=="DMSO"]
coreactualdata$Gene<-factor(coreactualdata$Gene, levels=rev(coremuimpactedEffects$Gene))
coreactualdata$Lineage<-factor(coreactualdata$Lineage, levels=c("CAMA1 Sen", "CAMA1 RiboR", "MCF7 Sen",  "MCF7 RiboR"))

ggplot( coreactualdata,
        aes(y= Gene, x= Lineage, fill=scalelog2cpmC))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Expression \n scaled(log2(x))",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CorePre_tx_ExpressionGeneDE_1_2.pdf", dpi=320, height=6, width=17)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CorePre_tx_ExpressionGeneDE_1_2.png", dpi=320, height=6, width=17)



RplotGenes <- strtcomp1$Gene
RpathFindisregHall<- data.table(ssGSEAmeta[gs_cat%in%c( "H")&gs_subcat%in%c( "" )][gene_symbol%in%RplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
RpathFindisregBio<- data.table(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%RplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]

RpathFindisregHall[propgenes>0.05]
RpathFindisregBio[propgenes>0.05][ngene>10]

