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
dd0mcf <-read.csv(file = "sh10050_MCF7_genes_gam_short.term_inputData_synergy.csv")
dd0cama<-read.csv(file = "sh11141_CAMA1_genes_gam_short.term_inputData_synergy.csv")

# Join data across cell lines
dd0 <- data.table(rbind(dd0cama%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) ),dd0mcf%>%dplyr::select( intersect(names(dd0cama),names(dd0mcf)) )))

# order factors
dd0$Treatment <- factor(dd0$Treatment,c("DMSO", "Afat","Ribo", "Comb"))
dd0$State <- factor(dd0$State,c("Sen","RiboR"))
dd0$CellLine <- factor(dd0$CellLine,c("CAMA1", "MCF7"))

# Generate day number column
dd0[, Day:=(Hour-as.numeric(as.character(timeunderdrug)))/24 ]

# Scale each genes expression data for each cell line (line & resistance)
dd0[,scalelog2cpm:= scale(Expression.log2cpm), by=c("Gene","CellLine","State")]
dd0[,scalelog2cpmB:= scale(Expression.log2cpm), by=c("Gene","CellLine")]

# Scale cell line expression of gene g relative to sensitive cells at time t under DMSO
dd0[,SensExpression.log2cpm := sum(Expression.log2cpm*(State=="Sen"&Treatment=="DMSO")), by=c("Gene","CellLine","Hour")]
dd0[,SensExpression.log2cpmB:= sum(Expression.log2cpm*(Treatment=="DMSO" )), by=c("Gene","CellLine","Hour","State")]
dd0[,NormSensExpression.log2cpm:=Expression.log2cpm -SensExpression.log2cpm, by=c("Gene","CellLine","Hour", "Treatment")]
dd0[,NormSensExpression.log2cpmB:=Expression.log2cpm -SensExpression.log2cpmB, by=c("Gene","CellLine","Hour", "Treatment")]

# Define color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(4)

# Generate and order the cell lineage column
dd0[,Lineage:=paste0(CellLine," ",State)]
dd0$Lineage <- factor(dd0$Lineage  , levels=c("CAMA1 Sen" ,"MCF7 Sen", "CAMA1 RiboR"  ,  "MCF7 RiboR"  ))

# List all genes
selectgenes<- unique(dd0$Gene)

# Look specifically at CDK cell cycle genes
plotthis<-dd0[Gene%in%paste0("CDK",1:19)]
plotthis$Gene<- factor(plotthis$Gene,
                       levels=c(unique(plotthis$Gene)[unique(plotthis$Gene)%in%paste0("CDK",1:9)],
                                unique(plotthis$Gene)[unique(plotthis$Gene)%in%paste0("CDK",10:19)]))
plotthis[,DMSOscalelog2cpm:= scalelog2cpm*(Treatment=="DMSO"), by=c("Gene","CellLine","Hour","State")]
plotthis[,Normscalelog2cpm:=scalelog2cpm-DMSOscalelog2cpm, by=c("Gene","CellLine","Hour","State")]

# Load ssGSEA pathwaygeneset metadata
ssGSEAmeta <- data.table( msigdbr(species = "Homo sapiens") )
ssGSEAmeta[gs_cat%in%c("C2","H")][gs_subcat%in%c("","CP:BIOCARTA","CP:REACTOME","CP:KEGG" )]$gs_name%>%unique()
ssGSEAmeta[gs_cat%in%c("C2","H")][gs_subcat%in%c("CP:BIOCARTA")]$gs_name%>%unique()
ssGSEAmeta[,ngenes:=length(unique(gene_symbol)), by=gs_name]

# Perform Linear mixed effects model endpoint comparison of differential expression between treatments and synergy and how this varies between resistant and sensitive lineages
finFULL <- rbindlist(mclapply(selectgenes, function(x){
  pred<-dd0[Gene %in% x][Day == max(Day)]
  #ggplot(pred,  aes(x=Treatment, y=NormSensExpression.log2cpm,col=State))+theme_classic()+ geom_point()+facet_grid(~CellLine,scales="free")+geom_hline(linetype=3,yintercept = 0)
  cat(x)
  if(length(unique(pred$CellLine))>1){
    m1 <- lmer(NormSensExpression.log2cpm~ -1 + State +  RiboTreated + AfatTreated + CombTreated +
                 State:RiboTreated + State:AfatTreated + State:CombTreated  + (1|Hour) , data=pred)
    mS <- lmer(NormSensExpression.log2cpm~ -1 +  RiboTreated + AfatTreated + CombTreated + (1|Hour) ,
               data=pred[State!="RiboR"])
    mR <- lmer(NormSensExpression.log2cpm~ 1 +  RiboTreated + AfatTreated + CombTreated  + (1|Hour),
               data=pred[State=="RiboR"])
    #summary(mS);summary(mR)
    outp1 <- cbind(Gene=x,Contrast="DeltaState_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )
    outpS <- cbind(Gene=x,Contrast= "Sen_Comb_vs_DMSO" , data.table( coef(summary(mS)) , keep.rownames=T ) )
    outpR <- cbind(Gene=x,Contrast="RiboR_Comb_vs_DMSO", data.table( coef(summary(mR)) , keep.rownames=T ) )
    #outmat <- rbind(outp1[rn=="StateRiboR:CombTreated"],outpS[rn=="CombTreated"],outpR[rn=="CombTreated"])
    outmat <- rbind(outp1 ,outpS ,outpR )
    return(outmat)
  }else{ NULL }
}, mc.cores =detectCores()-1 ))

finFULL2<- rbindlist(mclapply(selectgenes, function(x){
  pred<-dd0[Gene%in% x][Day ==max(Day)];  cat(x)
  if(length(unique(pred$CellLine))>1){
    m1 <- lmer(NormSensExpression.log2cpmB~ -1 +  RiboTreated + AfatTreated + CombTreated +
                 (0+RiboTreated|State) +
                 (0+AfatTreated|State) +
                 (0+CombTreated|State) +
                 (1|Hour) , data=pred)
    outmat <- cbind(Gene=x,Contrast="All_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )
    return(outmat)
  }else{ NULL }
}, mc.cores =detectCores()-1 ))

fincomp <- finFULL[rn%in%c("StateRiboR:CombTreated","CombTreated")]
fincomp2 <- finFULL2[rn%in%c("CombTreated")]
setnames(fincomp,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
setnames(fincomp2,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))

fincomp[,FDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
fincomp2[,FDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]

# synergistically DEGs between treatments in Resistant cells
fincompR <- fincomp[Contrast=="RiboR_Comb_vs_DMSO"][FDR<0.05][order(-(Estimate))]#[abs(Estimate)>log2(1.5)]
fincompR[, EstID:=1:nrow(fincompR)]
fincompR[, EstDir:=-1+2*(Estimate>0)]
# synergistically DEGs between treatments in Sensitive cells
fincompS <- fincomp[Contrast=="Sen_Comb_vs_DMSO"][FDR<0.05][order(-(Estimate))]#[abs(Estimate)>log2(1.5)]
fincompS[, EstID:=1:nrow(fincompS)]
fincompS[, EstDir:=-1+2*(Estimate>0)]

#View top results
fincompR[1:50]
fincompR[nrow(fincompR):(nrow(fincompR)-(50-1))]
fincompS[1:50]
fincompS[nrow(fincompR):(nrow(fincompR)-(50-1))]

# synergistically DEGs with differing effect size between Resistant and Sensitive cells
fincomp1 <- fincomp[Contrast=="DeltaState_Comb_vs_DMSO"][pvalue<0.05][order(-(Estimate))]#[abs(Estimate)>log2(1.5)]
fincomp1[, EstID:=1:nrow(fincomp1)]
fincomp1[, EstDir:=-1+2*(Estimate>0)]

# Identify synergistically differentially expressed genes in Resistant cells
RplotGenes <- fincompR$Gene
# Find what fraction of the gene sets in each biocarta, kegg, reactome or hallmark pathways are synergisticallydifferentially expressed genes in Resistant cells
RpathFindisregHall<- data.table(ssGSEAmeta[gs_cat%in%c( "H")&gs_subcat%in%c( "" )][gene_symbol%in%RplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
RpathFindisregBio<- data.table(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%RplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
RpathFindisregKegg<- data.table(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c( "CP:KEGG" )][gene_symbol%in%RplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
RpathFindisregReact<- data.table(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c( "CP:REACTOME" )][gene_symbol%in%RplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
RpathFindisregall<- data.table(ssGSEAmeta[gs_cat%in%c( "H","C2")&gs_subcat%in%c( "","CP:BIOCARTA" #,"CP:KEGG"
)][gene_symbol%in%RplotGenes]%>%group_by(gs_name,gs_subcat,gs_cat)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]

# Identify synergistically differentially expressed genes in Sensitive cells
SplotGenes <- fincompS$Gene
# Find what fraction of the gene sets in each biocarta, kegg, reactome or hallmark pathways are synergisticallydifferentially expressed genes in Sensitive cells
SpathFindisregHall<- data.table(ssGSEAmeta[gs_cat%in%c( "H")&gs_subcat%in%c( "" )][gene_symbol%in%SplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
SpathFindisregBio<- data.table(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%SplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
SpathFindisregKegg<- data.table(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c( "CP:KEGG" )][gene_symbol%in%SplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
SpathFindisregReact<- data.table(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c( "CP:REACTOME" )][gene_symbol%in%SplotGenes]%>%group_by(gs_name)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]
SpathFindisregall<- data.table(ssGSEAmeta[gs_cat%in%c( "H","C2")&gs_subcat%in%c( "","CP:BIOCARTA" #,"CP:KEGG"
)][gene_symbol%in%SplotGenes]%>%group_by(gs_name,gs_subcat,gs_cat)%>%dplyr::summarise(ngene=max(ngenes),propgenes=length(unique(gene_symbol))/max(ngenes)))[order(-propgenes)]


# review the top gene sets and find intersection between resistant and sensitive cells
SpathFindisregall[1:10]$gs_name
RpathFindisregall[1:10]$gs_name
RpathFindisregHall[1:10]$gs_name
SpathFindisregHall[1:10]$gs_name
SpathFindisregBio[1:20]$gs_name
RpathFindisregBio[1:20]$gs_name
RpathFindisregReact[1:20]$gs_name
RpathFindisregKegg[1:20]$gs_name
SpathFindisregReact[1:20]$gs_name
SpathFindisregKegg[1:20]$gs_name
intersect( RpathFindisregHall[1:10]$gs_name,SpathFindisregHall[1:10]$gs_name)
intersect( RpathFindisregBio[1:20]$gs_name,SpathFindisregBio[1:20]$gs_name)
intersect( RpathFindisregKegg[1:20]$gs_name,SpathFindisregKegg[1:20]$gs_name)


#plotGenesssGSEA <- intersect(ssGSEAmeta[gs_name%in%c(x)]$gene_symbol,fincomp$Gene)
# List a selected subset of genes from each gene set of interest to plot in heatmaps
lu_gen_lst <- data.table( Gene_Set=rep(c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","BIOCARTA_DEATH_PATHWAY","BIOCARTA_MITOCHONDRIA_PATHWAY"),each=5),
                          Gene= c("CDC25A","DEPDC1","CDK2","KIF4A","RRM2",
                                  "AURKA","CCNB2","CDK1","NEK2","TOP2A",
                                  "CASP3","CASP7","CASP8","CASP9","CASP10",
                                  "AIFM1","BAK1","BAX","BIK","DIABLO") )

# select genes to plot
subst <- merge(fincomp[Contrast=="RiboR_Comb_vs_DMSO"][order(Estimate)],lu_gen_lst, by="Gene")#[Gene%in%c(RAD51)]
PLOTthis<-merge( dd0[Gene%in% subst$Gene][Treatment!="DMSO"],lu_gen_lst, by="Gene")#[State=="RiboR"]

# Adjust annotations/labelling
PLOTthis$Gene <- factor(PLOTthis$Gene , levels= rev(subst$Gene ))
PLOTthis[  ,Gene_SetLab:="Hallmark \n E2F Targets"]
PLOTthis[Gene_Set=="HALLMARK_G2M_CHECKPOINT",Gene_SetLab:="Hallmark \n G2M Checkpoint"]
PLOTthis[Gene_Set=="HALLMARK_MITOTIC_SPINDLE",Gene_SetLab:="Hallmark \n Mitotic Spindle"]
PLOTthis[Gene_Set=="HALLMARK_MYC_TARGETS_V1",Gene_SetLab:="Hallmark \n MYC Targets V1"]
PLOTthis[Gene_Set=="BIOCARTA_ATRBRCA_PATHWAY",Gene_SetLab:="Biocarta \n ATR-BRCA Pathway"]
PLOTthis[Gene_Set=="BIOCARTA_DEATH_PATHWAY",Gene_SetLab:="Biocarta \n Death Pathway"]
PLOTthis[Gene_Set=="BIOCARTA_MITOCHONDRIA_PATHWAY",Gene_SetLab:="Biocarta \n Mitochondria Pathway"]
PLOTthis$Gene_SetLab<- factor(PLOTthis$Gene_SetLab,
                              levels=c("Hallmark \n E2F Targets",
                                       #"Hallmark \n Mitotic Spindle",
                                       "Hallmark \n G2M Checkpoint",
                                       #"Hallmark \n MYC Targets V1" ,
                                       #"Biocarta \n ATR-BRCA Pathway",
                                       "Biocarta \n Death Pathway",
                                       "Biocarta \n Mitochondria Pathway"))
PLOTthis[   ,LineageLab:="CAMA-1 \n Sensitive"]
PLOTthis[CellLine=="MCF7"& State=="Sen",LineageLab:="MCF-7 \n Sensitive"]
PLOTthis[CellLine=="CAMA1"& State=="RiboR",LineageLab:="CAMA-1 \n Resistant"]
PLOTthis[CellLine=="MCF7"& State=="RiboR",LineageLab:="MCF-7 \n Resistant"]
PLOTthis$LineageLab<- factor(PLOTthis$LineageLab, levels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant" ,"MCF-7 \n Sensitive","MCF-7 \n Resistant"  ))

# normalization 
PLOTthis[,ScaledNormSensExpression.log2cpm:=scale(NormSensExpression.log2cpm), by= c("Gene_SetLab","Gene")]
PLOTthis[,CelltypeScaledNormSensExpression.log2cpm:=scale(NormSensExpression.log2cpmB), by= c("Gene_SetLab","Gene","Lineage")]

# Heatmaps of expression over time by treatment
p1<-ggplot(PLOTthis[Treatment=="Ribo"],
           aes(x=as.factor(Hour), y=Gene, fill=CelltypeScaledNormSensExpression.log2cpm))+#NormSensExpression.log2cpm ))+
  theme_classic(base_size=26)+
  geom_tile()+
  facet_grid(Gene_SetLab~LineageLab,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  #scale_fill_gradient2(name="Expression \n log2FC \n vs \n sensitive cells \n under DMSO" , 
  #                     low="blue", high="darkred", mid="white",midpoint=0)+#,  limits=c(-4,4))+
  scale_fill_gradient2(name="Expression \n scaled log2FC \n vs DMSO" , 
                       low="blue", high="darkred", mid="white",midpoint=0, limits=c(-3.5,3.5))+
  theme(aspect.ratio=1)+#labs(title=gsub("_"," ",x), x="Hour")+
  theme(axis.text.x = element_text(size=15.5,angle = 45, vjust = 0.5))+
  labs(x="Time (Hours)")+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
p1
MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig5c ReactivationUnderRibo_SubsetSignifvsSensDMSO_v4col.png"), height=18.5, width=18.5, dpi=600)
#ggsave(file=paste0(MSfigureoutput,"Fig5c ReactivationUnderRibo_SubsetSignifvsSensDMSO_v4col.pdf"), height=18.5, width=18.5, dpi=600)

#ggsave(p1,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/ReactivationUnderRibo_SubsetSignifvsSensDMSO_v4col.png", dpi=320, height=18, width=18.5)
#ggsave(p1,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/ReactivationUnderRibo_SubsetSignifvsSensDMSO_v4col.pdf", dpi=320, height=18, width=18.5)


p1<-ggplot(PLOTthis[Treatment=="Comb"],
           aes(x=as.factor(Hour), y=Gene, fill=CelltypeScaledNormSensExpression.log2cpm))+#NormSensExpression.log2cpm ))+
  theme_classic(base_size=26)+
  geom_tile()+
  facet_grid(Gene_SetLab~LineageLab,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  #scale_fill_gradient2(name="Expression \n log2FC \n vs \n sensitive cells \n under DMSO" , 
  #                     low="blue", high="darkred", mid="white",midpoint=0, limits=c(-4,4))+
  scale_fill_gradient2(name="Expression \n scaled log2FC \n vs DMSO" , 
                       low="blue", high="darkred", mid="white",midpoint=0, limits=c(-3.5,3.5))+
  theme(aspect.ratio=1)+#labs(title=gsub("_"," ",x), x="Hour")+
  theme(axis.text.x = element_text(size=15.5,angle = 45, vjust = 0.5))+
  labs(x="Time (Hours)")+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

p1

#ggsave(file=paste0(MSfigureoutput,"Fig5c ReactivationUnderComb_SubsetSignifvsSensDMSO_v4col.png"), height=18.5, width=18.5, dpi=600)
#ggsave(file=paste0(MSfigureoutput,"Fig5c ReactivationUnderComb_SubsetSignifvsSensDMSO_v4col.pdf"), height=18.5, width=18.5, dpi=600)

#ggsave(p1,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/ReactivationUnderComb_SubsetSignifvsSensDMSO_v4col.png", dpi=320, height=18, width=18.5)
#ggsave(p1,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/ReactivationUnderComb_SubsetSignifvsSensDMSO_v4col.pdf", dpi=320, height=18, width=18.5)

p1<-ggplot(PLOTthis[Treatment=="Afat"],
           aes(x=as.factor(Hour), y=Gene, fill=CelltypeScaledNormSensExpression.log2cpm))+#NormSensExpression.log2cpm ))+
  theme_classic(base_size=26)+
  geom_tile()+
  facet_grid(Gene_SetLab~LineageLab,scales="free")+
  #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
  #scale_fill_gradient2(name="Expression \n log2FC \n vs \n sensitive cells \n under DMSO" , 
  #                     low="blue", high="darkred", mid="white",midpoint=0, limits=c(-4,4))+
  scale_fill_gradient2(name="Expression \n scaled log2FC \n vs DMSO" , 
                       low="blue", high="darkred", mid="white",midpoint=0, limits=c(-3.5,3.5))+
  theme(aspect.ratio=1)+#labs(title=gsub("_"," ",x), x="Hour")+
  theme(axis.text.x = element_text(size=15.5,angle = 45, vjust = 0.5))+
  labs(x="Time (Hours)")+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

p1
#ggsave(file=paste0(MSfigureoutput,"Fig5c ReactivationUnderAfat_SubsetSignifvsSensDMSO_v4col.png"), height=18.5, width=18.5, dpi=600)
#ggsave(file=paste0(MSfigureoutput,"Fig5c ReactivationUnderAfat_SubsetSignifvsSensDMSO_v4col.pdf"), height=18.5, width=18.5, dpi=600)

#ggsave(p1,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/ReactivationUnderAfat_SubsetSignifvsSensDMSO_v4col.png", dpi=320, height=18, width=18.5)
#ggsave(p1,file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/ReactivationUnderAfat_SubsetSignifvsSensDMSO_v4col.pdf", dpi=320, height=18, width=18.5)


# generate output file 
#fincomp[Contrast=="RiboR_Comb_vs_DMSO"][order(Estimate)][Gene%in%unique(ssGSEAmeta[gs_name%in% unique(lu_gen_lst$Gene_Set)]$gene_symbol)][pvalue<0.05]
#fincomp[Contrast=="RiboR_Comb_vs_DMSO"][order(Estimate)][Gene%in%unique(lu_gen_lst$Gene) ]$pvalue%>%range()
#fincomp[Contrast=="Sen_Comb_vs_DMSO"][order(Estimate)][Gene%in%unique(lu_gen_lst$Gene) ]$pvalue%>%range()
output<- finFULL[grepl("Sen_",Contrast)|grepl("RiboR_",Contrast)][Gene%in%unique(lu_gen_lst$Gene) ][order(Gene)]
setnames(output,old=c("Pr(>|t|)","t value","Std. Error"),new=c("pvalue","tvalue","Std.Error"))
#output[Gene=="AIFM1"][Contrast=="Sen_Comb_vs_DMSO"]
output[ , isminpval:= abs(tvalue)==max (abs(tvalue)), by=c("Gene")]
output[isminpval==T]$Gene%>%unique()
output[isminpval==T][order(abs(tvalue))]$pvalue%>%range()



# Load gene set : gen correlation assessment
ddGSmcf <- fread(file="sh10050_MCF7_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
ddGScama <- fread(file="sh11141_CAMA1_H.C2.C5.C6_gam_short.term_inputData_synergy.csv")
ddGS <- data.table(rbind(ddGScama%>%dplyr::select( intersect(names(ddGScama),names(ddGSmcf)) ),ddGSmcf%>%dplyr::select( intersect(names(ddGScama),names(ddGSmcf)) )))

# List focal gene sets
GSlist<-c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT","BIOCARTA_MITOCHONDRIA_PATHWAY",  "BIOCARTA_DEATH_PATHWAY")
ddGSforgenecor <- ddGS[Gene_Set%in%GSlist]

# Scaling and normalization
ddGSforgenecor[,scaleSSGSEA:= scale(Enrichment_Score), by=c("Gene_Set","CellLine","State")]
ddGSforgenecor[,scaleSSGSEAB:= scale(Enrichment_Score), by=c("Gene_Set","CellLine")]
ddGSforgenecor[,SensExpression.SSGSEA := sum(Enrichment_Score*(State=="Sen"&Treatment=="DMSO")), by=c("Gene_Set","CellLine","Hour")]
ddGSforgenecor[,SensExpression.SSGSEAB:= sum(Enrichment_Score*(Treatment=="DMSO" )), by=c("Gene_Set","CellLine","Hour","State")]
ddGSforgenecor[,NormSensExpression.SSGSEA:=Enrichment_Score -SensExpression.SSGSEA, by=c("Gene_Set","CellLine","Hour", "Treatment")]
ddGSforgenecor[,NormSensExpression.SSGSEAB:=Enrichment_Score -SensExpression.SSGSEAB, by=c("Gene_Set","CellLine","Hour", "Treatment")]

#order treatments and lineages
ddGSforgenecor$Treatment <- factor(ddGSforgenecor$Treatment,c("DMSO", "Afat","Ribo", "Comb")) 
ddGSforgenecor[,Lineage:=paste0(CellLine," ",State)]
ddGSforgenecor$Lineage <- factor(ddGSforgenecor$Lineage  , levels=c("CAMA1 Sen" ,"MCF7 Sen", "CAMA1 RiboR"  ,  "MCF7 RiboR"  ))

# plot pathway changes over time
ggplot(ddGSforgenecor,
       aes(x=Hour/24, y=NormSensExpression.SSGSEA, col= Treatment, group=interaction(State,Gene_Set,Treatment) ))+
  theme_classic(base_size=12)+
  geom_hline(linetype="dashed", aes(yintercept=0))+
  geom_point(size=3)+ geom_line(linewidth=1.2)+
  facet_grid( Gene_Set~Lineage, scales="free_y")+
  # scale_color_manual(name="Treatment",values=cols[-1])+# , labels=c("Ribociclib","Afatinib","Combination"))
  #  scale_fill_manual(name="Treatment",values=cols[-1])+
  scale_shape_discrete(name="Cell line", labels=c("CAMA-1", "MCF7"))+
  theme(aspect.ratio=1)+#+#labs(title=gsub("_"," ",x), x="Hour")+
  labs(y= "ssGSEA pathway activity \n (vs  sensitive cell under DMSO)", x="Day")+
  scale_x_continuous(breaks=c(0,6,12,18) )+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

#ddGSforgenecor%>%dplyr::select(Gene_Set , NormSensExpression.SSGSEA, Cell:day,Lineage)



# For each of 4 focal pathways, assess gene set's correlation with ssGSEA pathway scores
##focal pathway 1 
# Merge pathway and gene data
#dd0[Gene%in%unique(ssGSEAmeta[gs_name%in% GSlist[1]]$gene_symbol)]
#ssGSEAmeta[gs_name%in% GSlist[1]]
GSplusGene<-merge(
  ddGSforgenecor[Gene_Set==GSlist[1]]%>%dplyr::select(Gene_Set , NormSensExpression.SSGSEA, Cell,day,Lineage, Treatment),
  dd0[Treatment!="DMSO"][Gene%in% c("CDK2","E2F3","E2F4",  unique(ssGSEAmeta[gs_name%in% GSlist[1]]$gene_symbol))] %>%select(Gene,Expression.log2cpm,Cell), 
  by=c("Cell")
)
#Test linear association of gene expression and pathway activity 
GS1res<-rbindlist(lapply(1:length(unique(GSplusGene$Gene)), function(i){
  data.table(Gene=unique(GSplusGene$Gene)[i],
             Gene_Set=unique(GSplusGene$Gene_Set),
             coef(summary(
               lm(NormSensExpression.SSGSEA~Expression.log2cpm*Treatment+Lineage ,
                  data= GSplusGene[Gene==unique(GSplusGene$Gene)[i]]))), keep.rownames = T)
}))[grepl("Expression.log2cpm",rn)]
setnames(GS1res,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
GS1res[ , isminpval:= abs(tvalue)==max (abs(tvalue)), by=c("Gene")]
GS1res[isminpval==T]
GS1res[isminpval==T][Gene %in% lu_gen_lst[Gene_Set==GSlist[1]]$Gene]

##focal pathway 1 
GSplusGene2<-merge(
  ddGSforgenecor[Gene_Set==GSlist[2]]%>%dplyr::select(Gene_Set , NormSensExpression.SSGSEA, Cell,day,Lineage, Treatment),
  dd0[Treatment!="DMSO"][Gene%in% c( unique(ssGSEAmeta[gs_name%in% GSlist[2]]$gene_symbol))] %>%select(Gene,Expression.log2cpm,Cell), 
  by=c("Cell")
)
GS2res<-rbindlist(lapply(1:length(unique(GSplusGene2$Gene)), function(i){
  data.table(Gene=unique(GSplusGene2$Gene)[i],
             Gene_Set=unique(GSplusGene2$Gene_Set),
             coef(summary(
               lm(NormSensExpression.SSGSEA~Expression.log2cpm*Treatment+Lineage ,
                  data= GSplusGene2[Gene==unique(GSplusGene2$Gene)[i]]))), keep.rownames = T)
}))[grepl("Expression.log2cpm",rn)]
setnames(GS2res,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
GS2res[ , isminpval:= abs(tvalue)==max (abs(tvalue)), by=c("Gene")]
GS2res[isminpval==T]
GS2res[isminpval==T][Gene %in% lu_gen_lst[Gene_Set==GSlist[2]]$Gene]

##focal pathway 1 
GSplusGene3<-merge(
  ddGSforgenecor[Gene_Set==GSlist[3]]%>%dplyr::select(Gene_Set , NormSensExpression.SSGSEA, Cell,day,Lineage, Treatment),
  dd0[Treatment!="DMSO"][Gene%in% c( unique(ssGSEAmeta[gs_name%in% GSlist[3]]$gene_symbol))] %>%select(Gene,Expression.log2cpm,Cell), 
  by=c("Cell")
)
GS3res<-rbindlist(lapply(1:length(unique(GSplusGene3$Gene)), function(i){
  data.table(Gene=unique(GSplusGene3$Gene)[i],
             Gene_Set=unique(GSplusGene3$Gene_Set),
             coef(summary(
               lm(NormSensExpression.SSGSEA~Expression.log2cpm*Treatment+Lineage ,
                  data= GSplusGene3[Gene==unique(GSplusGene3$Gene)[i]]))), keep.rownames = T)
}))[grepl("Expression.log2cpm",rn)]#[rn!="(Intercept)"];# GS3res[,rn:=NULL]
setnames(GS3res,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
GS3res[ , isminpval:= abs(tvalue)==max (abs(tvalue)), by=c("Gene")]
GS3res[isminpval==T]
GS3res[isminpval==T][Gene %in% lu_gen_lst[Gene_Set==GSlist[3]]$Gene]

##focal pathway 1 
GSplusGene4<-merge(
  ddGSforgenecor[Gene_Set==GSlist[4]]%>%dplyr::select(Gene_Set , NormSensExpression.SSGSEA, Cell,day,Lineage, Treatment),
  dd0[Treatment!="DMSO"][Gene%in% c( unique(ssGSEAmeta[gs_name%in% GSlist[4]]$gene_symbol))] %>%select(Gene,Expression.log2cpm,Cell), 
  by=c("Cell")
)
GS4res<-rbindlist(lapply(1:length(unique(GSplusGene4$Gene)), function(i){
  data.table(Gene=unique(GSplusGene4$Gene)[i],
             Gene_Set=unique(GSplusGene4$Gene_Set),
             coef(summary(
               lm(NormSensExpression.SSGSEA~Expression.log2cpm*Treatment+Lineage ,
                  data= GSplusGene4[Gene==unique(GSplusGene4$Gene)[i]]))), keep.rownames = T)
}))[grepl("Expression.log2cpm",rn)]#[rn!="(Intercept)"]; GS4res[,rn:=NULL]
setnames(GS4res,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
GS4res[ , isminpval:= abs(tvalue)==max (abs(tvalue)), by=c("Gene")]
GS4res[isminpval==T]
GS4res[isminpval==T][Gene %in% lu_gen_lst[Gene_Set==GSlist[4]]$Gene]

# Join results from across pathways
restab<-rbind( GS1res[isminpval==T][Gene %in% lu_gen_lst[Gene_Set==GSlist[1]]$Gene],
GS2res[isminpval==T][Gene %in% lu_gen_lst[Gene_Set==GSlist[2]]$Gene],
GS3res[isminpval==T][Gene %in% lu_gen_lst[Gene_Set==GSlist[3]]$Gene],
GS4res[isminpval==T][Gene %in% lu_gen_lst[Gene_Set==GSlist[4]]$Gene])
restab[,rn:=NULL];restab[,isminpval:=NULL];

# generate output file
restab$pvalue%>%max()
#write.csv( restab, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Table stats for fig 5c heatmap cell cycle and apoptosis gene correlation with gene sets under treatment raw.csv")
#ggplot( GSplusGene4[Gene=="CASP8"], aes(NormSensExpression.SSGSEA,Expression.log2cpm ,col=Treatment))+geom_point()
#ggplot( GSplusGene[Gene=="CCNE1"], aes(NormSensExpression.SSGSEA,Expression.log2cpm ,col=Treatment))+geom_point()


# Visualize gene expression over time using heatmaps
intersect(SpathFindisregBio[propgenes>=0.1][]$gs_name,RpathFindisregBio[propgenes>=0.1][]$gs_name)
intersect(SpathFindisregHall[propgenes>=0.1][]$gs_name,RpathFindisregHall[propgenes>=0.1][]$gs_name)
sloc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/GeneExpressionTimeCourse/HeatmapPathwayGeneTimecourse RelativeToSenDMSO/"
pltgenesheatNorm <- function(x,save=F){
  plotGenesssGSEA <- intersect(ssGSEAmeta[gs_name%in%c(x)]$gene_symbol,fincomp$Gene)
  PLOTthis<-dd0[Gene%in% plotGenesssGSEA][Treatment!="DMSO"]#[State=="RiboR"]
  PLOTthis$Gene <- factor(PLOTthis$Gene , levels= rev(fincomp[Contrast=="RiboR_Comb_vs_DMSO"][order(Estimate)][Gene%in%plotGenesssGSEA]$Gene ))
  PLOTthis[,CelltypeScaledNormSensExpression.log2cpm:=scale(NormSensExpression.log2cpmB), by= c("Gene","Lineage")]
  p1<-ggplot(PLOTthis,
             aes(x=as.factor(Hour), y=Gene, fill=CelltypeScaledNormSensExpression.log2cpm ))+
    theme_classic(base_size=16)+
    geom_tile()+
    facet_grid(Treatment~Lineage,scales="free")+
    #facet_grid(.~interaction(Lineage,Treatment),scales="free")+
    scale_fill_gradient2(name="Expression \n scaled log2FC \n vs DMSO" , low="blue", high="darkred", mid="white")+
    theme(aspect.ratio=3)+labs(title=gsub("_"," ",x), x="Hour")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
    scale_y_discrete(guide=guide_axis(n.dodge=3)) 
  p1
  if(save==T){ ggsave(p1,file=paste0(sloc, file="FinalCombovsDMSOnormalized_",x,".pdf"), dpi=320, height=18, width=18)
  }else{p1}
}

# Examine cell-cycle pathways showing reactivation
pltgenesheatNorm("BIOCARTA_MCM_PATHWAY" ,save=F)
pltgenesheatNorm("BIOCARTA_CDC25_PATHWAY" ,save=F)
pltgenesheatNorm("BIOCARTA_SRCRPTP_PATHWAY" ,save=F)
pltgenesheatNorm("BIOCARTA_RB_PATHWAY" ,save=F)
pltgenesheatNorm("BIOCARTA_EFP_PATHWAY" ,save=F)
pltgenesheatNorm("BIOCARTA_G2_PATHWAY" ,save=F)
pltgenesheatNorm("BIOCARTA_PTC1_PATHWAY" ,save=F)
pltgenesheatNorm("BIOCARTA_RANMS_PATHWAY" ,save=F)
syngslist<-c("BIOCARTA_CDC25_PATHWAY" ,"BIOCARTA_SRCRPTP_PATHWAY" ,"BIOCARTA_RB_PATHWAY" ,
             "BIOCARTA_EFP_PATHWAY" ,"BIOCARTA_G2_PATHWAY" ,"BIOCARTA_PTC1_PATHWAY" ,"BIOCARTA_RANMS_PATHWAY",
             "HALLMARK_E2F_TARGETS" ,    "HALLMARK_G2M_CHECKPOINT" , "HALLMARK_MITOTIC_SPINDLE")
unique(ssGSEAmeta[gs_name%in%syngslist ]$gene_symbol)


# synergistically differentially expressed genes in Resistant cells
fincomp[Contrast=="RiboR_Comb_vs_DMSO"][FDR<0.1][order(Estimate)][Gene%in%sort(unique(ssGSEAmeta[gs_name%in%syngslist ]$gene_symbol))]



# Identify genes in biocarta pathways that are synergistically differentially (impacted) expressed in either Resistant or Sensitive cells
ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%c(RplotGenes,SplotGenes)]
impacted<-unique(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%c(RplotGenes,SplotGenes)]$gene_symbol)
#impacted<-unique(ssGSEAmeta[gs_cat%in%c( "H","C2")&gs_subcat%in%c( "" , "CP:BIOCARTA")][gene_symbol%in%c(RplotGenes,SplotGenes)]$gene_symbol)

# Subset synergy effect results from GAM trend data for these genes
impactedEffects<-fincomp[Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%impacted][order(-Estimate)]
impactedEffects$Gene<- factor(impactedEffects$Gene, levels=unique(impactedEffects$Gene))

# Summarize mean synergy effect size estimate across cell lines and rank genes by effect size
muimpactedEffects <-data.table(impactedEffects%>%group_by(Gene)%>%dplyr::summarise(Estimate=mean(Estimate)))[order(-Estimate)]
muimpactedEffects$Gene<- factor(muimpactedEffects$Gene, levels=unique(muimpactedEffects$Gene))
muimpactedEffects[,rankabsEst:=rank(-abs(Estimate))]
ggplot( impactedEffects, 
        aes(y= Gene, x= Estimate,col=Contrast))+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size=3)+
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
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/GeneSpecificTreatmentSynergyByEnd.pdf", dpi=320, height=6, width=6)


#labs(x="Synergistic effect")
ggplot( muimpactedEffects[rankabsEst<=50], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)
ggplot( muimpactedEffects[c(1:20, (nrow(muimpactedEffects)-19):nrow(muimpactedEffects))], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)

# Select 30 genes to focus on in more detail (top synergistic treatment targets)
lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[c(1:2, (nrow(muimpactedEffects)-27):nrow(muimpactedEffects))]$Gene ))
# Use LMER to quantify synergy within each cell line, using all samples on the final treatment day (end point) and accounting for repeated measures with a random intercept for each hour 
outsyn<- rbindlist(lapply(1:nrow(lutab), function(i){
  tmp<- dd0[Gene%in% lutab[i]$Gene][Day ==max(Day)][Lineage==lutab[i]$Lineage][timeunderdrug%in%c(0,6,24)]
  m1 <- lmer(NormSensExpression.log2cpm~ -1 +  RiboTreated + AfatTreated + CombTreated +
               (1|Hour) , data=tmp)
  cbind(Gene=lutab[i]$Gene,Lineage=lutab[i]$Lineage,Contrast="DeltaState_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )[rn=="CombTreated"]
}))

# Produce heatmap of synergistic treatment effect on cell cycle genes in each cell line
outsyn$Lineage <- factor(outsyn$Lineage, levels=c("CAMA1 Sen" , "CAMA1 RiboR" , "MCF7 Sen" , "MCF7 RiboR") )
ggplot( outsyn,
        aes(y= Gene, x= Lineage, fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/GeneSpecificTreatmentSynergyByEndPerCellLine.pdf", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/GeneSpecificTreatmentSynergyByEndPerCellLine.png", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/GeneSpecificTreatmentSynergyByEndPerCellLine.tiff", dpi=320, height=6, width=6)


#unique(ssGSEAmeta[gs_name%in%syngslist][gene_symbol%in%c(RplotGenes,SplotGenes)]$gene_symbol)




### focus on cell cycle and apoptosis pathways
ssGSEAmeta[gs_name %in%c("BIOCARTA_DEATH_PATHWAY","BIOCARTA_MITOCHONDRIA_PATHWAY")]
ssGSEAmeta[gs_name %in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")]

#ccgenes <- (ssGSEAmeta[gs_name %in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT"
#                                     )]$gene_symbol)%>%unique()%>%sort()
# list genes in cell cycle-related pathways
ccgenes <-data.table(unique(ssGSEAmeta[gs_name %in%c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT",
                                                     "HALLMARK_MITOTIC_SPINDLE","BIOCARTA_RACCYCD_PATHWAY",
                                                     "BIOCARTA_MCM_PATHWAY",
                                                     "BIOCARTA_CELLCYCLE_PATHWAY",
                                                     "BIOCARTA_G1_PATHWAY","BIOCARTA_G2_PATHWAY"
)]%>%dplyr::select(gs_name,gene_symbol))%>%group_by(gene_symbol)%>%summarise(n=n()))$gene_symbol%>%sort()

# list genes in cell-death-related pathways
apgenes <- (ssGSEAmeta[gs_name %in%c("BIOCARTA_DEATH_PATHWAY","BIOCARTA_MITOCHONDRIA_PATHWAY"#,
                                     #"HALLMARK_APOPTOSIS"#, 
                                     #"HALLMARK_P53_PATHWAY"
)][]$gene_symbol)%>%unique()%>%sort()

# Extract statistics for cell cycle and cell death related genes
subcc<-finFULL[Gene%in%c(ccgenes)]
setnames(subcc,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
subap<-finFULL[Gene%in%c(apgenes)]
setnames(subap,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))

#ggplot(subcc[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][pvalue<0.05],aes(y=Estimate,x=rn))+geom_boxplot()+geom_point()+
#  facet_wrap(~Contrast)+coord_flip() #+geom_violin()
#ggplot(subap[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][pvalue<0.05],aes(y=Estimate,x=rn))+geom_boxplot()+
 # geom_point()+  facet_wrap(~Contrast)+coord_flip() 
#fincomp2
## Cell death
# Summarise the mean median max and min synergy effects across cell lines for each cell death related gene
finFULLEffs<-data.table(finFULL[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%apgenes]%>%group_by(Gene)%>%
                          mutate(                                 maxeff=max(Estimate),
                                                                  mineff=min(Estimate),
                                                                  medeff=median(Estimate),
                                                                  meaneff=mean(Estimate),
                                                                  meansynergeff=mean((rn=="CombTreated")*Estimate)
                          ))[
                            order(meaneff)]

# Order genes, and annotate resistance state and treatment
finFULLEffs$Gene<- factor(finFULLEffs$Gene,levels=unique(finFULLEffs$Gene))
finFULLEffs[, CellLine :="Resistant"]
finFULLEffs[Contrast=="Sen_Comb_vs_DMSO", CellLine :="Sensitive"]
finFULLEffs$CellLine<- factor(finFULLEffs$CellLine,levels=c("Sensitive","Resistant"))
finFULLEffs$rn<- factor(finFULLEffs$rn,levels=c("RiboTreated","AfatTreated","CombTreated"))
setnames(finFULLEffs,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))

# Bar plot of treatment effect sizes
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
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~CellLine)+coord_flip() +
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificMainTreatmentAndSynergyByEnd.pdf", dpi=320, height=6, width=6)

# plot gene expression changes in cell death genes post treatment with ribociclib,afatinib or combination
apgenetop <- unique(finFULLEffs[c(1:(10*6) ,(nrow(finFULLEffs)-(c(1:(10*6) )-1)) )]$Gene)
plotdd0<-dd0[Gene%in% apgenetop][Day ==max(Day)][Treatment!="DMSO"]
plotdd0$Gene<- factor(plotdd0$Gene,levels=unique(finFULLEffs$Gene))
ggplot(plotdd0,aes(y=NormSensExpression.log2cpmB,x=Gene,
                   col=Treatment))+  theme_classic()+
  geom_boxplot()+
  geom_point(aes(group=Treatment),position=position_dodge(width=1), stat="identity")+
  facet_wrap(~State)+coord_flip() 


## Cell cycle
# Summarise the mean median max and min synergy effects across cell lines for each cell cycle related gene
ccfinFULLEffs<-data.table(finFULL[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%ccgenes]%>%group_by(Gene)%>%
                            mutate(maxeff=max(Estimate),
                                   mineff=min(Estimate),
                                   medeff=median(Estimate),
                                   meaneff=mean(Estimate),
                                   meansynergeff=mean((rn=="CombTreated")*Estimate)
                            ))[
                              order(meaneff)]
setnames(ccfinFULLEffs,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))

# Order genes, and annotate resistance state and treatment
ccfinFULLEffs[,minpval:=min(pvalue),by="Gene"]
ccfinFULLEffs$Gene<- factor(ccfinFULLEffs$Gene,levels=rev(unique(ccfinFULLEffs$Gene)))
ccfinFULLEffs[, CellLine :="Resistant"]
ccfinFULLEffs[Contrast=="Sen_Comb_vs_DMSO", CellLine :="Sensitive"]
ccfinFULLEffs$CellLine<- factor(ccfinFULLEffs$CellLine,levels=c("Sensitive","Resistant"))
ccfinFULLEffs$rn<- factor(ccfinFULLEffs$rn,levels=c("RiboTreated","AfatTreated","CombTreated"))
setnames(ccfinFULLEffs,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))

# Bar plot of treatment effect sizes
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

# plot gene expression changes in cell cycle genes post treatment with ribociclib,afatinib or combination
ccgenetop<-unique(ccfinFULLEffs[c(1:(15*6) ,(nrow(ccfinFULLEffs)-(c(1:(15*6) )-1)) )]$Gene)
ccplotdd0<-dd0[Gene%in% ccgenetop][Day ==max(Day)][Treatment!="DMSO"]
ccplotdd0$Gene<- factor(ccplotdd0$Gene,levels=unique(ccfinFULLEffs$Gene))
ggplot(ccplotdd0,aes(y=NormSensExpression.log2cpmB,x=Gene,
                     col=Treatment))+  theme_classic()+
  geom_boxplot()+
  geom_point(aes(group=Treatment),position=position_dodge(width=1), stat="identity")+
  facet_wrap(~State)+coord_flip() 



# caspase gene specific plot
ggplot(finFULL[Contrast!="DeltaState_Comb_vs_DMSO"][grepl("CASP",Gene)],aes(y=Estimate,x=Gene,fill=rn))+
  theme_classic()+
  geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~Contrast)+coord_flip() 


##subset results for cell death genes in resistant cells
apstatsR<-fincomp[Contrast=="RiboR_Comb_vs_DMSO"][Gene%in%apgenes]
apstatsR[,apFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
apfincompR <- apstatsR[pvalue<0.05][order(-(Estimate))]
apfincompR[, EstID:=1:nrow(apfincompR)]
apfincompR[, EstDir:=-1+2*(Estimate>0)]

##subset results for cell death genes in sensitive cells
apstatsS<-fincomp[Contrast=="Sen_Comb_vs_DMSO"][Gene%in%apgenes]
apstatsS[,apFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
apfincompS <- apstatsS[pvalue<0.05][order(-(Estimate))]
apfincompS[, EstID:=1:nrow(apfincompS)]
apfincompS[, EstDir:=-1+2*(Estimate>0)]

##subset results for cell cycle genes in resistant cells
ccstatsR<-fincomp[Contrast=="RiboR_Comb_vs_DMSO"][Gene%in%ccgenes]
ccstatsR[,ccFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
ccfincompR <- ccstatsR[pvalue<0.05][order(-(Estimate))]
ccfincompR <- ccstatsR[ccFDR<0.05][order(-(Estimate))]
ccfincompR[, EstID:=1:nrow(ccfincompR)]
ccfincompR[, EstDir:=-1+2*(Estimate>0)]

##subset results for cell cycle genes in sensitive cells
ccstatsS<-fincomp[Contrast=="Sen_Comb_vs_DMSO"][Gene%in%ccgenes]
ccstatsS[,ccFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
ccfincompS <- ccstatsS[pvalue<0.05][order(-(Estimate))]
ccfincompS <- ccstatsS[ccFDR<0.05][order(-(Estimate))]
ccfincompS[, EstID:=1:nrow(ccfincompS)]
ccfincompS[, EstDir:=-1+2*(Estimate>0)]

##subset results for cell cycle genes across sensitive and resistant cells
ccstatsRS<-fincomp2[Contrast=="All_Comb_vs_DMSO"][Gene%in%ccgenes]
ccstatsRS[,ccFDR:=p.adjust(pvalue, method="fdr"),by=c("rn","Contrast")]
ccfincompRS <- ccstatsRS[pvalue<0.05][order(-(Estimate))]
ccfincompRS <- ccstatsRS[ccFDR<0.05][order(-(Estimate))]
ccfincompRS[, EstID:=1:nrow(ccfincompRS)]
ccfincompRS[, EstDir:=-1+2*(Estimate>0)]
ccfincompRS




#List genes synergisticall impacted in either resistant, sensitive or both cell states
ccRplotGenes <- ccfincompR$Gene
ccSplotGenes <- ccfincompS$Gene
ccRSplotGenes <- ccfincompRS$Gene
ccimpacted<-unique(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%c(ccRplotGenes,ccSplotGenes)]$gene_symbol)
#ccimpacted<-unique(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%c(ccRSplotGenes)]$gene_symbol)

# subset statistics for these synergistically impacted cell cycle  genes
ccimpactedEffects<-fincomp[Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%ccimpacted][order(-Estimate)]
ccimpactedEffects$Gene<- factor(ccimpactedEffects$Gene, levels=unique(ccimpactedEffects$Gene))
ccmuimpactedEffects <-data.table(ccimpactedEffects%>%group_by(Gene)%>%dplyr::summarise(Estimate=mean(Estimate)))[order(-Estimate)]
ccmuimpactedEffects$Gene<- factor(ccmuimpactedEffects$Gene, levels=unique(ccmuimpactedEffects$Gene))
ccmuimpactedEffects[,rankabsEst:=rank(-abs(Estimate))]

# Plot effect sizes as dotpot and heatmap
ggplot( ccimpactedEffects, 
        aes(y= Gene, x= Estimate,col=Contrast))+geom_vline(xintercept = 0,linetype="dashed")+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_point()+ labs(x="Synergistic effect")+
  scale_color_npg(name="Resistance", labels=c("Resistant","Sensitive"))

ggplot( ccmuimpactedEffects[rankabsEst<=20], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_void(base_size=14)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+
  geom_text(aes(label=Gene ))+
  theme(axis.ticks=element_blank(),axis.text.y=element_blank())+
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEnd.pdf", dpi=320, height=6, width=6)

# collect results
cc_resout <- ccimpactedEffects%>%mutate(Effect="SynegryByEndpoint", Process="CellCycle")
#write.csv(cc_resout)

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


# Get gene expression data for the identified cell cycle genes
cclutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=ccmuimpactedEffects[]$Gene ))

#Linear mixed effects model testing synergistic treatment effects on normalized log2CPM expression
ccoutsyn<- rbindlist(lapply(1:nrow(cclutab), function(i){
  tmp<- dd0[Gene%in% cclutab[i]$Gene][Day ==max(Day)][Lineage==cclutab[i]$Lineage][timeunderdrug%in%c(0,6,24)]
  m1 <- lmer(NormSensExpression.log2cpm~ -1 +  RiboTreated + AfatTreated + CombTreated +
               (1|Hour) , data=tmp)
  cbind(Gene=cclutab[i]$Gene,Lineage=cclutab[i]$Lineage,Contrast="DeltaState_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )[rn=="CombTreated"]
}))
#setnames(ccoutsyn,old=c("Pr(>|t|)","t value"),new=c("pvalue","tvalue"))
setnames(ccoutsyn,old=c("Pr(>|t|)","Std. Error","t value"),new=c("pvalue","Std.Error","tvalue"))
ccoutsyn$Lineage <- factor(ccoutsyn$Lineage, levels=c("CAMA1 Sen" , "CAMA1 RiboR" , "MCF7 Sen" , "MCF7 RiboR") )
cc_resout2 <- ccoutsyn%>%mutate(Effect="SynegryByEndpoint", Process="CellCycle")
#write.csv(cc_resout2)

# Visualize synergy effects on each cell line (estimated by LME model)
ggplot( ccoutsyn,
        aes(y= Gene, x= Lineage, fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=3)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))


#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSI.pdf", dpi=320, height=12, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSI.png", dpi=320, height=12, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSI.tiff", dpi=320, height=12, width=6)

# Show results for core cell cycle related gene set
corecc<- c("TOP2A","BIRC5","CDC25C","CDK1","CCNB1","CCNB2","AURKA","MYC","CDC25A","CDK2","CCNE1","CDK4","CCND1","E2F1","CDKN1A"
          #"ORC1","STMN1","MCM7","TPX2","PLK1",
)

ggplot( ccoutsyn[Gene%in%corecc],
        aes(y= Gene, x= Lineage, fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergy \n effect",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig4c CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.png"), height=6, width=6, dpi=600)
#ggsave(file=paste0(MSfigureoutput,"Fig4c CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.pdf"), height=6, width=6, dpi=600)

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.pdf", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.png", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellCycleGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.tiff", dpi=320, height=6, width=6)


##Repeat process for apoptosis genes
RplotGenes <- apfincompR$Gene
SplotGenes <- apfincompS$Gene


#Identify impacted genes and gen statistics
#impacted<-unique(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%c(RplotGenes,SplotGenes)]$gene_symbol)
impacted<-unique(ssGSEAmeta[gs_cat%in%c( "C2")&gs_subcat%in%c(  "CP:BIOCARTA")][gene_symbol%in%apgenes]$gene_symbol)
impactedEffects<-fincomp[Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%impacted][order(-Estimate)]
impactedEffects$Gene<- factor(impactedEffects$Gene, levels=unique(impactedEffects$Gene))
muimpactedEffects <-data.table(impactedEffects%>%group_by(Gene)%>%dplyr::summarise(Estimate=mean(Estimate)))[order(-Estimate)]
muimpactedEffects$Gene<- factor(muimpactedEffects$Gene, levels=rev(unique(muimpactedEffects$Gene)))
muimpactedEffects[,rankabsEst:=rank(-abs(Estimate))]

ap_resout <- impactedEffects%>%mutate(Effect="SynegryByEndpoint", Process="CellDeath")
#write.csv(cc_resout)

#Plot GAM coefficients
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


ggplot( muimpactedEffects[rankabsEst<=30], 
        aes(y= Gene, x= "", fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=2)+labs(x="")+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)

#Linear mixed effects model testing synergistic treatment effects on normalized log2CPM expression
lutab<- data.table(expand.grid(Lineage=unique(dd0$Lineage),Gene=muimpactedEffects[]$Gene ))
outsyn<- rbindlist(lapply(1:nrow(lutab), function(i){
  tmp<- dd0[Gene%in% lutab[i]$Gene][Day ==max(Day)][Lineage==lutab[i]$Lineage][timeunderdrug%in%c(0,6,24)]
  m1 <- lmer(NormSensExpression.log2cpm~ -1 +  RiboTreated + AfatTreated + CombTreated +
               (1|Hour) , data=tmp)
  cbind(Gene=lutab[i]$Gene,Lineage=lutab[i]$Lineage,Contrast="DeltaState_Comb_vs_DMSO", data.table( coef(summary(m1)) , keep.rownames=T ) )[rn=="CombTreated"]
}))
outsyn$Lineage <- factor(outsyn$Lineage, levels=c("CAMA1 Sen" , "CAMA1 RiboR" , "MCF7 Sen" , "MCF7 RiboR") )
setnames(outsyn,old=c("Pr(>|t|)","Std. Error","t value"),new=c("pvalue","Std.Error","tvalue"))

#Gather results
ap_resout2 <- outsyn%>%mutate(Effect="SynegryByEndpoint", Process="CellDeath")
#write.csv(cc_resout2)

#Visualize cell line specific response (obtained from LME model)
ggplot( outsyn[pvalue<0.05],
        aes(y= Gene, x= Lineage, fill=Estimate))+
  theme_classic(base_size=12)+
  theme(aspect.ratio=3)+
  geom_tile()+ 
  scale_fill_gradient2(name="Synergistic \n effect",low="blue", high="darkred", mid="white",midpoint=0)+
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

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
  scale_x_discrete(name="", labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive",  "MCF-7 \n Resistant"))+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig4c CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.png"), height=6, width=6, dpi=600)
#ggsave(file=paste0(MSfigureoutput,"Fig4c CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.pdf"), height=6, width=6, dpi=600)

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.pdf", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.png", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificTreatmentSynergyByEndPerCellLineSELECTED.tiff", dpi=320, height=6, width=6)


# resoutputboth1 <- rbind(
#   cc_resout%>%dplyr::select(Effect,Process,Gene,Estimate, df,tvalue,pvalue ),
#   ap_resout%>%dplyr::select(Effect,Process,Gene,Estimate, df,tvalue,pvalue )
# )
# write.csv(resoutputboth1[pvalue<0.05],file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/tables/
#.csv")

resoutputboth <- rbind(
  cc_resout2%>%dplyr::select(Effect,Process,Lineage,Gene,Estimate,Std.Error, df,tvalue,pvalue ),
  ap_resout2%>%dplyr::select(Effect,Process,Lineage,Gene,Estimate,Std.Error, df,tvalue,pvalue )
)

# write.csv(resoutputboth[pvalue<0.05],file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/tables/CellCycleCellDeathSynergyByEOT.csv")
actualdata <- dd0[Gene%in%c(coreap,corecc)]
resout<-data.table(actualdata%>%dplyr::select(Gene,Lineage,CellLine,State,Treatment, Hour,Expression.cpm,NormSensExpression.log2cpmB))[order(Gene,State,CellLine)]
#write.csv(resout, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/tables/SynergyEOTCellcycleCellDeathGeneExpression.csv")

# Look at gene specific trajectories over time
ggplot(resout[Hour==max(Hour)], 
       aes(col= Treatment, x= Lineage, y=log2(Expression.cpm)))+
  theme_classic(base_size=14)+
  theme(aspect.ratio=1.5)+labs(x="")+
  scale_x_discrete(labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive", "MCF-7 \n Resistant"))+
  geom_point()+
  facet_wrap(~Gene,scales="free_y")#scale_fill_gradient2(name="Expression \n scaled(log2(x))",low="blue", high="darkred", mid="white",midpoint=0)

ggplot(resout[Hour==max(Hour)], 
       aes(col= Treatment, x= Lineage, y=NormSensExpression.log2cpmB))+
  theme_classic(base_size=14)+
  theme(aspect.ratio=1.5)+labs(x="")+
  scale_x_discrete(labels=c("CAMA-1 \n Sensitive","CAMA-1 \n Resistant","MCF-7 \n Sensitive", "MCF-7 \n Resistant"))+
  geom_point()+
  facet_wrap(~Gene,scales="free_y")#scale_fill_gradient2(name="Expression \n scaled(log2(x))",low="blue", high="darkred", mid="white",midpoint=0)




## Gather statistics across cell cycle and cell death genes to combine into one plot 
# show effect sizes of each treatment on each gene
finFULLEffs<-data.table(finFULL[rn!="(Intercept)"][Contrast!="DeltaState_Comb_vs_DMSO"][Gene%in%c(corecc,coreap)]%>%group_by(Gene)%>%
                          mutate(                                 maxeff=max(Estimate),
                                                                  mineff=max(Estimate),
                                                                  medeff=median(Estimate),
                                                                  meaneff=mean(Estimate),
                                                                  meansynergeff=mean((rn=="CombTreated")*Estimate)
                          ))[
                            order(meaneff)]

#Annotation and ordering
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

#Visualize treatment effects as barplot (facet by theme and resistance state)
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
  theme(aspect.ratio=1)+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))


MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig4d CellDeathGeneSpecificMainTreatmentAndCombinationEffectsByEndwithoutpoints.png"), width=6, height=6, dpi=320)
#ggsave(file=paste0(MSfigureoutput,"Fig4d CellDeathGeneSpecificMainTreatmentAndCombinationEffectsByEndwithoutpoints.pdf"), width=6, height=6, dpi=320)


#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Top15CellCYCLEandCellDeathGeneSpecificMainTreatmentAndSynergyByEnd.pdf", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/Top15CellCYCLEandCellDeathGeneSpecificMainTreatmentAndSynergyByEnd.png", dpi=320, height=6, width=6)



ggplot(finFULLEffs[  ][],aes(y=Estimate,x=Gene,fill=rn#,alpha=pvalue<0.05
))+  theme_classic()+
  #geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~CellLine)+coord_flip() +
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Synergy"))+
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificMainTreatmentAndSynergyByEnd.pdf", dpi=320, height=6, width=6)



# create a version overlaying the estimated treatment effects and the observed data
plotdd0<-dd0[Gene%in% rev(c(corecc[],coreap))][Day ==max(Day)][Treatment!="DMSO"]
plotdd0$Gene<- factor(plotdd0$Gene,levels=rev(c(corecc,coreap)))
plotdd0[,Theme:="Cell Death"]
plotdd0[Gene%in%corecc,Theme:="Cell Cycle"]
plotdd0[, CellLineRes :="Resistant"]
plotdd0[State=="Sen", CellLineRes :="Sensitive"]

ggplot(plotdd0,aes(y=NormSensExpression.log2cpmB,x=Gene,
                   col=Treatment,group=interaction(Lineage,Treatment,Gene)))+  theme_classic()+
  geom_boxplot(aes(fill=CellLine,group=interaction(Lineage,Treatment,Gene)))+
  geom_point(aes(group=interaction(Lineage,Treatment,Gene)),position=position_dodge(width=.75), stat="identity")+
  facet_grid(Theme~CellLineRes,scales = "free_y")+coord_flip() +
  scale_fill_manual(name="Cell \n line",values=c("grey","white"))+
  scale_color_manual(name="Treatment \n effect", 
                     values=gg_color_hue(4)[-1],
                     labels=c("Ribociclib","Afatinib", "Combination"))+
  labs(y="Post-tx expression \n (log2(CPM) normalized to sensitive cells under DMSO)")+
  theme(aspect.ratio=1)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificMainTreatmentAndComboExpressionByEnd.pdf", dpi=320, height=6, width=6)


finFULLEffs[,Treatment:="Comb"]
finFULLEffs[rn=="AfatTreated",Treatment:="Afat"]
finFULLEffs[rn=="RiboTreated",Treatment:="Ribo"]
finFULLEffs$Treatment<- factor(finFULLEffs$Treatment, levels=c("DMSO", "Afat", "Ribo", "Comb"))
finFULLEffs[,CellLineRes:=CellLine]

# Modified plotting aesthetics
ggplot(plotdd0,aes(y=NormSensExpression.log2cpmB,x=Gene,
                   col=Treatment,fill=Treatment))+  theme_classic(base_size=22)+
  geom_bar(data=finFULLEffs[  ][],aes(y=Effect,x=Gene,fill=Treatment,col=Treatment),position="dodge",width=1, stat="identity")+
  #geom_boxplot(aes(shape=CellLine,group=interaction(Lineage,Treatment,Gene)))+
  geom_point(aes(group=interaction(Treatment)), #interaction(Lineage,Treatment,Gene,Hour)),
             position=position_dodge(width=.85), col=1, pch=21, alpha=0.4,size=1.2,stat="identity" )+
  facet_grid(Theme~CellLineRes,scales = "free_y")+coord_flip() +
  scale_shape_manual(name="Cell \n line",values=c("grey","white"))+
  scale_color_manual(name="Treatment \n effect", 
                     values=gg_color_hue(4)[-1],
                     labels=c("Ribociclib","Afatinib", "Combination"))+
  labs(y="Post-tx expression \n (log2(CPM) on day 18 compared to DMSO)")+
  theme(aspect.ratio=1.618)+
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Combination")) +
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
#ggsave(file=paste0(MSfigureoutput,"Fig4d CellDeathGeneSpecificMainTreatmentAndCombinationEffectsByEnd.png"), width=11, height=10.5, dpi=320)
#ggsave(file=paste0(MSfigureoutput,"Fig4d CellDeathGeneSpecificMainTreatmentAndCombinationEffectsByEnd.pdf"), width=11, height=10.5, dpi=320)


#Save output
#write.csv(finFULLEffs
tabl4d<-finFULLEffs
setnames(tabl4d, old=c("Std. Error"), new=c("Std.Error"))
tabl4d_out<-tabl4d%>%dplyr::select(Gene,Theme,Treatment, CellLineRes,Effect,Estimate, Std.Error, df,     tvalue ,      pvalue)
#write.csv( tabl4d_out, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Table stats for fig 4d bar and scatterplot cell cycle and apoptosis gene expression under treatment raw2.csv")



# Modified plotting aesthetics
plotdd0$CellLineRes <- factor( plotdd0$CellLineRes, levels=c("Sensitive" ,"Resistant"))
ggplot(plotdd0,aes(y=NormSensExpression.log2cpmB,x=Gene,
                   col=Treatment,fill=Treatment))+  theme_classic(base_size=22)+
  geom_bar(data=finFULLEffs[  ][],aes(y=Effect,x=Gene,fill=Treatment,col=Treatment),position="dodge",width=1, stat="identity")+
  #geom_boxplot(aes(shape=CellLine,group=interaction(Lineage,Treatment,Gene)))+
  geom_point(aes(shape=CellLine,group=interaction(Treatment)), #interaction(Lineage,Treatment,Gene,Hour)),
             position=position_dodge(width=.85), col=1,  alpha=0.6,size=1.02,stat="identity" )+#pch=21,
  facet_grid(Theme~CellLineRes,scales = "free_y")+coord_flip() +
  scale_shape_manual(name="Cell \n line",labels=c("CAMA-1", "MCF-7"),values= c(21,24))+
  scale_color_manual(name="Treatment \n effect", 
                     values=gg_color_hue(4)[-1],
                     labels=c("Ribociclib","Afatinib", "Combination"))+
  labs(y="Post-tx expression \n (log2(CPM) on day 18 compared to DMSO)")+
  theme(aspect.ratio=1)+
  scale_fill_manual(name="Treatment \n effect", 
                    values=gg_color_hue(4)[-1],
                    labels=c("Ribociclib","Afatinib", "Combination")) +
  theme( axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))+
  guides(shape = guide_legend(override.aes = list(size=2,shape = c(21,24)) ) )



MSfigureoutput <- "/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Manuscript Figures/"
ggsave(file=paste0(MSfigureoutput,"Fig4d CellDeathGeneSpecificMainTreatmentAndCombinationEffectsByEndSquare.png"), width=11, height=11, dpi=600)
ggsave(file=paste0(MSfigureoutput,"Fig4d CellDeathGeneSpecificMainTreatmentAndCombinationEffectsByEndSquare.pdf"), width=11, height=11, dpi=600)

#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificMainTreatmentAndCombinationEffectsByEnd.pdf", dpi=320, height=6, width=6)
#ggsave(file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/R Analysis/ComboRiboAfatinib/Fig3/CellDeathGeneSpecificMainTreatmentAndCombinationEffectsByEnd.png", dpi=320, height=6, width=6)


