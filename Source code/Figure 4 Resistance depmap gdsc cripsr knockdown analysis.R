#CRISPR-Cas9 gene effect scores represent baseline gene dependency measured in untreated cell lines and were integrated with independent drug sensitivity profiles from GDSC to assess associations with resistance phenotypes.
#CRISPR screens show baseline gene dependencies
# Drug response data (GDSC1/2) shows phenotypic sensitivity
# By comparing dependencies between ribociclib-sensitive vs resistant lines, we identify:
#  genes that become critical for survival in the resistant state
#  i.e. rewired dependencies, even without drug present

rm(list=ls())
wd<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/U01_gene_express_timecourse/depmap_data/"

# Load necessary packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(msigdbr)
library(readr)
library(data.table)
require(glmnet)
require(pheatmap)
require(ggsci)

require(lme4)
require(lmerTest)

# Specify candidate genes of interest: derived from temporal analysis of RNAseq trajectories
genes_of_interest <- c(
  "MELK","TK1","CKS1B","TUBB","MCM4","CDK4","CDC27","E2F4","ORC2",
  "CDK1","CDK2","CCNB1","CCNB2", "CCNE1","CDC25A","CDC25B", "CDC25C", "E2F1"  ,"E2F3"  ,"TOP2A","CDKN1A","AUKA","NEK2","KIF4A","DEPDC1",
  "RRM2","BIRC5","CDK6",#"MYC",
  "BAX", "BAK1","BIK", "DIABLO", "AIFM1", "CASP3", "CASP7", "CASP8", "CASP9","CASP10",
  "EGFR", "ERBB2", "ERBB3","ERBB4"
)

# Construct dataframe linking gene to pathway
gene_groups<-data.table(Gene=genes_of_interest)
gene_groups[, pathway:="Cell cycle"]
gene_groups[Gene%in%c("BAX", "BAK1","BIK", "DIABLO", "AIFM1", "CASP3", "CASP7", "CASP8", "CASP9","CASP10") , pathway:="Apoptosis"]
gene_groups[Gene%in%c("EGFR", "ERBB2", "ERBB3","ERBB4") , pathway:="ERBB"]

# Load msigDB
msigdb<-data.table(msigdbr())[gs_cat %in% c("H","C2")][grepl("HALLMARK",gs_name)|grepl("BIOCARTA",gs_name)|grepl("REACTOME",gs_name)|grepl("KEGG",gs_name)]

# Analyse Pan cancer data or restrict to BRCA
pancan<- T#F

# Function to clean column names once, globally
clean_gene_names <- function(x) {
  gsub(" \\(.*\\)$", "", x)
}

# Load CRISPR gene effect data
gene_effect <- data.table(fread(paste0(wd,"CRISPRGeneEffect.csv") ) )
gene_effect[1:4,1:14]
setnames(gene_effect, old="V1", new="ModelID")
colnames(gene_effect) <- clean_gene_names(colnames(gene_effect))

# CRISPR dependency probabilities
gene_dependency <- fread(paste0(wd,"CRISPRGeneDependency.csv") )
setnames(gene_dependency, old="V1", new="ModelID")
colnames(gene_dependency) <- clean_gene_names(colnames(gene_dependency))

# Model metadata
model <- fread(paste0(wd, "Model.csv") )

# Load cell line metadata
metadata <- fread(paste0(wd,"ModelCondition.csv")  )
screen_map <- fread(paste0(wd, "CRISPRScreenMap.csv") )

# Contains cell line name, lineage, ER status, etc.
if(pancan){
  # Define pan cancer cell lines
  er_breast_models <- model[] %>%
    dplyr::select(ModelID, CellLineName, OncotreeSubtype)
  breast_models <- model[] %>%
    dplyr::select(ModelID, CellLineName, OncotreeSubtype)
}else{
  # Define ER+ breast cancer cell lines
  er_breast_models <- model[OncotreeLineage == "Breast"][
    (!grepl("HER+",ModelSubtypeFeatures)&grepl("ER+",ModelSubtypeFeatures))|
      grepl("luminal",ModelSubtypeFeatures)
  ] %>%
    dplyr::select(ModelID, CellLineName, OncotreeSubtype)
  breast_models <- model[OncotreeLineage == "Breast"][] %>%
    dplyr::select(ModelID, CellLineName, OncotreeSubtype)
}

# Subset CRISPR data to pan cancer or breast cancer cell lines
gene_effect_breast <- gene_effect[ ModelID %in% breast_models$ModelID]  #ModelID %in% er_breast_models$ModelID] 
gene_dependency_breast <- gene_dependency [ ModelID %in% breast_models$ModelID]  #ModelID %in% er_breast_models$ModelID] 

# List all gene names
all_genes <- names(gene_effect_breast)[-1]

# Subset data for candidate genes
gene_effect_breast_sub <- gene_effect_breast #%>% dplyr::select(ModelID, any_of(genes_of_interest))
gene_dependency_breast_sub <- gene_dependency_breast# %>%  dplyr::select(ModelID, any_of(genes_of_interest))


# Downloaded Filtered GDSC1/2 for BRCA cell lines 
# drug effect anova results : https://www.cancerrxgene.org/downloads/anova
# Files downloaded from CancerRxGene portal (GDSC v1/v2), BRCA-filtered
gdsc1 <-  fread(paste0(wd,"GDSC_BRCA_ANOVA_Tue Feb  3 23_32_54 2026.csv"))
gdsc2 <-  fread(paste0(wd,"GDSC_BRCA_ANOVA_Tue Feb  3 23_32_59 2026.csv"))
gdsc0 <- rbind(gdsc1,gdsc2)
cdk_drug_list <-   names(table(gdsc0[`Drug target`=="CDK4, CDK6"]$`Drug name`))
erbb_drug_list <-  names(table(gdsc0[`Target Pathway`=="EGFR signaling"]$`Drug name`))

# AUC values per cell line: https://www.cancerrxgene.org/downloads/drug_data?tissue=BRCA
if(pancan){
  gdsc_cellcycle <-  fread(paste0(wd,"GDSC_PANCANCER_IC_Wed Feb  4 11_28_05 2026.csv"))
  gdsc_erbb <-  fread(paste0(wd,"GDSC_PANCANCER_IC_Wed Feb  4 11_28_11 2026.csv"))
  gdsc <- rbind(gdsc_cellcycle,gdsc_erbb)
}else{
  gdsc <-  fread(paste0(wd,"GDSC_BRCA_IC_Tue Feb  3 23_42_37 2026.csv"))
}
gdsc_res <- gdsc[`Drug Name`%in%c("Ribociclib","Afatinib")] %>% 
  dplyr::select(c("Drug Name":"Cell Line Name", "TCGA Classification", "Tissue", "Tissue Sub-type", "AUC", "Dataset Version" ))
setnames( gdsc_res, old= c("Cell Line Name"), new= c("CellLineName")) 

# Merge everything
plot_data <- merge(
  data.table( 
  gene_effect_breast_sub %>%
    left_join(
      merge(breast_models,gdsc_res, by="CellLineName"),
      by = "ModelID") %>%
    pivot_longer(
      cols = any_of(all_genes),#any_of(genes_of_interest),
      names_to = "Gene",
      values_to = "CERES"
    ) )[is.finite(AUC)] ,
  gene_groups, by="Gene", all.x=T)[!is.na(CERES)]


# Tidy the data
plot_data[is.na(pathway),pathway:="Other"]
plot_data[,FocalGene:=TRUE]
plot_data[is.na(pathway),FocalGene:=FALSE]

# remove some columns that are no longer needed: 
plot_data[, ModelID:=NULL]
plot_data[, `Drug ID`:=NULL]
plot_data[, Tissue:=NULL]
plot_data[, `Tissue Sub-type`:=NULL]
plot_data[, `Dataset Version`:=NULL]
plot_data[, TCGAClassification:=`TCGA Classification`]
plot_data[, `TCGA Classification`:=NULL]

# classify cell lines as resistant or sensitive to each drug (Resistant= logit(AUC) > mean )
plot_data[, Sensitivity := ifelse(boot::logit(AUC) > mean(boot::logit(AUC)), "Resistant", "Sensitive"), by="Drug Name"]#plot_data[, Sensitivity := ifelse(AUC > median(AUC), "Resistant", "Sensitive"), by="Drug Name"]
plot_data$Sensitivity <- factor(plot_data$Sensitivity, levels=c("Sensitive","Resistant" ))
plot_data[, DrugName:=`Drug Name`]
plot_data[, logitAUC:=boot::logit(AUC)]


##### Pathway essentiality bootstrap analysis : Q) are these genes more essential in breast cancer than average.
# For each communication pathway, day and cell type signal senders and receiver, randomize the tumor response annotations to remove the response related network structure
# Repeat many times to determine the distribtution of communication differences that may be expected by chance under the null model that communication strength do not differ between resistant and sensitive tumors.
# number of repeats
nRandomizations <- 1000

# Focal gene set mean essentiality
summaryEssentiality <- data.table( plot_data[TCGAClassification=="BRCA"&`Drug Name`=="Ribociclib"][Gene%in%c(genes_of_interest)] %>% 
                            group_by(pathway) %>% 
                            dplyr::summarise(meannegCERES= mean(-CERES) ) )

# Randomized gene set mean essentiality
MultisummaryRandGenesets <- rbindlist( #mc
  lapply(1:nRandomizations, function(i){
    set.seed(i)
    genes_not_interest<- all_genes[sample(which(!all_genes%in%genes_of_interest),length(genes_of_interest))]
    stats_data<- plot_data[TCGAClassification=="BRCA"&`Drug Name`=="Ribociclib"][Gene%in%c(genes_of_interest,genes_not_interest)]
    stats_data$pathway <- factor(stats_data$pathway  , levels=c("Other","Apoptosis" ,"Cell cycle", "ERBB") )
    mod1 <- lm(-CERES~pathway, data=stats_data)
    result<- cbind(sim_id= i,data.table( coef( summary(mod1) ) , keep.rownames = T))
    setnames(result, old=c("rn","Std. Error","t value","Pr(>|t|)"), new=c("pathway","Std.Error","tvalue","pvalue"))
    return(result)
  })) [pathway!="(Intercept)"] 
MultisummaryRandGenesets[, pathway:= gsub("pathway","",pathway)]     

# Compute mean and sd of gene essentiality estimates by pathway
rand_meandiff <- MultisummaryRandGenesets[, .(rand.meandiff = mean(Estimate)), by=pathway]
rand_sddiff   <- MultisummaryRandGenesets[, .(rand.sddiff = sd(Estimate)), by=pathway]

# Calculate sest statistics and p-values
z_scores <- merge(summaryEssentiality, rand_meandiff, by="pathway")
z_scores[, mean_nullrandomized:= meannegCERES - rand.meandiff]
z_scores <- merge(z_scores, rand_sddiff, by="pathway")
z_scores[, obs.z := (rand.meandiff)/rand.sddiff ]
z_scores[, p_z := 2 * pnorm(-abs(obs.z))]

# two-sided, testing whether pathway genes are more or less essential than bootstrap randomized gene sets
# Pathway-level essentiality was assessed using a bootstrap randomazation framework in which observed gene essentiality (-CERES scores) were compared against 1,000 randomly sampled gene sets of equal size. Significance was quantified using standardized z-scores.
ggplot(MultisummaryRandGenesets[pathway!="(Intercept)"], aes(x=pathway, fill=pathway, y=Estimate)) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(outlier.shape = NA) + # fill="lightblue"
  stat_boxplot(geom ='errorbar') +
  geom_jitter(width=0.04,alpha=0.1)+
  theme_classic(base_size=12) +
  labs(y="Gene set essentiality \n (versus randomized null gene sets)", 
       x="Pathway", 
       #title="Pathway essentiality vs random gene sets"
       ) +
  theme(aspect.ratio=1,
        #axis.text.x = element_text(angle=90, hjust=1), 
        legend.position = "none")+
  scale_fill_npg()+
  theme(text = element_text(color = "black"))

ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/Gene_essentiality_v2.pdf", width=4, height=4,dpi=320)
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/Gene_essentiality_v2.png", width=4, height=4,dpi=320)
# biological interpretation: these pathways are unuasally essential at baseline
unique(plot_data[TCGAClassification=="BRCA"&`Drug Name`=="Ribociclib"]$CellLineName)%>%length()


##### Prediction of ribociclib drug resistanc (AUC) from baseline gene essentiality
#Tests predictive relevance, not just association, to show that focal genes are informative
wide <-  data.table(
    plot_data[`Drug Name`=="Ribociclib"][Gene %in% c(all_genes[])]%>%
      group_by(Gene)%>%
      mutate(scale_CERES= ( CERES-mean(CERES, na.rm=T))/ sd(CERES, na.rm=T))%>%
      select(-c(CERES,pathway))%>%
      spread(Gene,scale_CERES) 
    ) 

# Specify inputs and outputs
X <- as.matrix( wide %>% select(any_of(all_genes))  )
X[is.na(X)] <- 0
X_focal <- X[,colnames(X)%in%genes_of_interest]
Y <- boot::logit(wide[`Drug Name`=="Ribociclib"]$AUC)

# Fit Lasso models using all genes vs focal gene set 
set.seed(12345)
cvfit <- cv.glmnet( y=Y, x= X, nfolds=10,alpha = 1, scale=F,type.measure ="mse",keep =T)
cvfit_focal <- cv.glmnet( y=Y, x= X_focal, nfolds=10,alpha = 1, scale=F,type.measure ="mse",keep =T)
cvfit_focal_1se <- cv.glmnet(y=Y, x=X_focal, alpha=1, nfolds=10,scale=F,type.measure ="mse",keep =T)
cvfit_focal_en <- cv.glmnet(y=Y, x=X_focal, alpha=0.5, nfolds=10,scale=F,type.measure ="mse",keep =T)

# idetify significant predictors
coef_lasso <- coef(cvfit, s = "lambda.min")
coef_lasso_focal <- coef(cvfit_focal, s = "lambda.min")
selected_genes <- setdiff(
  rownames(coef_lasso)[which(coef_lasso != 0)], 
  "(Intercept)"
  )
selected_genes_focal <- setdiff(
  rownames(coef_lasso_focal)[which(coef_lasso_focal != 0)], 
  "(Intercept)"
)

# Cross-validated prediction performance
pred <- predict(cvfit, X, s = "lambda.min")[,1]
pred_focal <- predict(cvfit_focal, X_focal, s = "lambda.min")[,1]
pred_cv <- cvfit$fit.preval[, cvfit$lambda == cvfit$lambda.min]
pred_cv_focal <- cvfit_focal$fit.preval[, cvfit_focal$lambda == cvfit_focal$lambda.min]


pred_cv_1se <- cvfit_focal_1se$fit.preval[, cvfit_focal_1se$lambda == cvfit_focal_1se$lambda.1se]
pred_cv_focal_en <- cvfit_focal_en$fit.preval[, cvfit_focal_en$lambda == cvfit_focal_en$lambda.min]
cor(pred_cv_1se, Y)
cor(pred_cv_focal_en, Y)

ggplot(data.table( pred=pred ,
                   pred_focal=pred_focal,
                   logitAUC=as.vector(Y)) , 
       aes(x=boot::inv.logit(pred), y=boot::inv.logit( logitAUC)) )+
  geom_abline(linetype="dashed")+
  theme_classic()+
  geom_point(alpha=0.5)+
  theme(aspect.ratio=1)+#coord_flip()+
  labs(y="Observed AUC", x="Predicted AUC") 

ggplot(data.table( pred= pred ,
                   pred_focal= pred_focal,
                   logitAUC= as.vector(Y)) , 
       aes(x= boot::inv.logit(pred), y= boot::inv.logit( logitAUC)) )+
  geom_abline(linetype= "dashed")+
  theme_classic()+
  geom_point(alpha= 0.5)+
  geom_point(aes( x= boot::inv.logit(pred_focal) ),alpha=0.5, col="red")+
  theme(aspect.ratio=1)+#coord_flip()+
  labs(y="Observed AUC", x="Predicted AUC") 

cor(pred, Y);cor(pred_focal, Y)
cor(pred_cv, Y);cor(pred_cv_focal, Y)

rmse_res<-data.table(
  GeneSet= c("All genes", "Focal genes"), #c("All \n genes", "Focal \n genes"),
  rmse= c( sqrt(mean((pred - Y)^2))   , 
           sqrt(mean((pred_focal - Y)^2))  ), # RMSE
  cor= c( cor(pred_cv, Y)  , 
          cor(pred_cv_focal, Y) ) # RMSE
)

# Contrast the predictive accuracy of the full and focal gene set models  
ggplot(rmse_res , 
       aes(x= GeneSet, y= cor, fill=GeneSet ) )+
  theme_classic()+theme(aspect.ratio=1, legend.position="none")+
  geom_bar(stat="identity")+
  #geom_text(aes(label=sprintf("RMSE: %.2f", rmse)), vjust=1.5) +
  labs(y="Resistance prediction accuracy \n (observed vs predicted AUC correlation)", x="Lasso gene set predictors") +
  theme(text = element_text(color = "black"))
# Focal genes retain predictive signal with reduced variance, consistent with a parsimonious resistance signature.
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/PredictivePerformance_V2.pdf", width=4, height=4,dpi=320)
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/PredictivePerformance_V2.png", width=4, height=4,dpi=320)
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/PredictivePerformance2_V2.pdf", width=3.5, height=4,dpi=320)
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/PredictivePerformance2_V2.png", width=3.5, height=4,dpi=320)

## re-fit the models and determine the stable predictors of ribociclib resistance
iterativ_fits <- rbindlist( lapply(1:100, function(i){
  set.seed(i)
  cvfit_focal <- cv.glmnet( y=Y, x= X_focal, nfolds=10, alpha = 1, scale=F, type.measure ="mse", keep =T)
  coef_lasso_focal <- coef(cvfit_focal, s = "lambda.min")
  selected_genes_focal <- rownames(coef_lasso_focal)[which(coef_lasso_focal != 0)]
  selected_genes_focal <- setdiff(selected_genes_focal, "(Intercept)")
  return( data.table(i=i,
                     Gene= rownames(coef_lasso_focal), 
                     coef=as.vector(coef_lasso_focal)  )   )
}))
iterativ_fits[, maxabseff:= max(abs(coef)), by=Gene]
iterativ_fits[, zero_not_in_quantile:=!( quantile(coef , 0.25)<= 0 & quantile(coef , 0.75) >= 0), by=Gene]
iterativ_fitsplot <- merge(gene_groups,  iterativ_fits[Gene!="(Intercept)"][maxabseff>0][zero_not_in_quantile==T], by="Gene")
ggplot(  iterativ_fitsplot,
        aes(y=coef, x=Gene))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_boxplot(aes(fill=pathway))+
  theme_classic(base_size=12)+
  geom_point(alpha=0.5, size=0.4)+
  coord_flip()+ 
  theme(aspect.ratio=1)+
  scale_fill_npg()+
  theme(text = element_text(color = "black"))+
  labs(fill="Pathway",y="Lasso resistance prediction coefficients \n (across cross-validation iterations)")
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/LassoCoefficients.pdf", width=5, height=4,dpi=320)
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/LassoCoefficients.png", width=5, height=4,dpi=320)



##### Compare gene essentiallity in Resistant vs sensitive cell lines
# Visualization to show resistant-specific dependencies
sensit_lm_mod <- lm(-CERES ~ -1+Gene+Sensitivity:Gene, data= plot_data[Gene%in%genes_of_interest][DrugName=="Ribociclib"])
signifRvsS <- data.table(coef( summary(sensit_lm_mod)), keep.rownames=T)[grep("SensitivityResistant",rn)][`Pr(>|t|)`<0.05]
signifRvsS[, Gene := gsub(":SensitivityResistant","",rn)]
signifRvsS[, Gene := gsub("Gene","",Gene)]
deltaRdata <- data.table( plot_data[pathway != "Other" ] %>%
  group_by(Sensitivity, Gene, pathway, DrugName) %>%
  dplyr::summarise(CERES=median(CERES)) %>%
  spread( Sensitivity, CERES)%>%
  mutate(deltaR=Resistant-Sensitive) )
deltaRdata$Gene <- factor(deltaRdata$Gene,
                          levels=  deltaRdata[DrugName=="Ribociclib"][order(deltaR)]$Gene)
deltaRdata$DrugName <- factor(deltaRdata$DrugName, 
                              levels=  rev(sort(unique( deltaRdata[order(deltaR)]$DrugName)) ))
#
ggplot(data= deltaRdata,#[pathway=="Cell cycle"],
       aes(y=-deltaR, x=Gene, col=DrugName))+
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point()+
  #facet_wrap(~pathway )+
  coord_flip()+
  theme( aspect.ratio=1)+
  labs(y = "Differential median gene essentiality\n(resistant − sensitive cell lines)",color = "Drug")

 
##### Does baseline dependency on this gene significantly associate with drug resistance?
mod_resist_Ribo<- lm(I(-CERES)~-1+Gene+Gene:logitAUC , 
          data=plot_data[`Drug Name`=="Ribociclib"][Gene%in%c(genes_of_interest)]
          )
mod_resist_Afatin<- lm(I(-CERES)~-1+Gene+Gene:logitAUC , 
              data=plot_data[`Drug Name`=="Afatinib"][Gene%in%c(genes_of_interest)]
)

result_tab_Ribo <- cbind(DrugName="Ribociclib",as.data.table(coef(summary(mod_resist_Ribo)), keep.rownames = T))
result_tab_Afat <- cbind(DrugName="Afatinib",as.data.table(coef(summary(mod_resist_Afatin)), keep.rownames = T))
result_tab <- rbind(result_tab_Ribo,result_tab_Afat)
result_tab[, effect:="gene_essent"]
result_tab[grepl("AUC", rn), effect:="resistance"]
result_tab[, Gene:=gsub("Gene","", gsub(":I\\(boot::logit\\(AUC\\)\\)", "", rn) )]
result_tab[, Gene:=gsub("Gene","", gsub(":\\logitAUC", "", rn) )]
setnames(result_tab, old=c("Std. Error",  "t value", "Pr(>|t|)" ), new=c("Std.Error",  "tvalue", "pvalue"))
result_tab[effect=="resistance"][order(Estimate)]
result_tab[,drug_effect:=paste0(DrugName,effect)]

wide_result_tab <-data.table(result_tab%>%select(DrugName,Gene,effect,Estimate)%>%spread(effect,Estimate))[order(resistance)]
result_tab_reord<- result_tab
result_tab_reord$Gene <- factor(result_tab_reord$Gene, levels=wide_result_tab[DrugName=="Ribociclib"]$Gene)
result_tab_reord$DrugName <- factor(result_tab_reord$DrugName, levels=c("Ribociclib","Afatinib"))

# significance annotations
result_tab_reord[, sig := ""]
result_tab_reord[effect == "resistance" & pvalue < 0.05,  sig := "*"]
result_tab_reord[effect == "resistance" & pvalue < 0.01,  sig := "**"]
result_tab_reord[effect == "resistance" & pvalue < 0.001, sig := "***"]

lim <- max(abs(result_tab_reord[effect=="resistance"]$Estimate), na.rm = TRUE)
ggplot(result_tab_reord[effect=="resistance"], 
       aes(x = DrugName, y = Gene, fill=Estimate ))+#, group=Gene, col=Gene)) +
  geom_tile()+
  geom_text(
    aes(label = sig),
    color = "white",
    size = 4
  ) +
  theme_classic(base_size = 12) +
  labs(  x = "Drug",  fill="Gene \n essentiality \n change \n with \n resistance"  ) +
  theme(
    aspect.ratio=1#,
    #legend.position = "none"
  )+
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    limits = c(-lim, lim)
  )

#Asterisks denote significance of the gene × resistance interaction term in linear models of baseline gene dependency: * p < 0.05, ** p < 0.01, *** p < 0.001.
#Significant gene–resistance associations identified by linear modeling are indicated directly on the heatmap


#Plot gene essentiality (CERES)
scatter_dat <- plot_data[Gene%in%unique(result_tab[effect=="resistance"][pvalue<0.05]$Gene)]
scatter_order<- unique(result_tab[DrugName=="Ribociclib"][effect=="resistance"][pvalue<0.05][order(Estimate)]$Gene)
scatter_order<-unique(c(scatter_order, result_tab[effect=="resistance"][pvalue<0.05][order(Estimate)]$Gene))
scatter_dat$Gene <- factor(scatter_dat$Gene, levels=scatter_order )
ggplot(scatter_dat, aes(x = AUC, y = -CERES, group=Gene, col=Gene)) +
  geom_point( size = 0.5, alpha = 0.8) +
  geom_smooth(se=T,method="lm", col="black", fill="black", alpha=0.4)+
  #geom_smooth(se=F,method="gam", formula=y~s(x,k=4),method.args = list(family = "betar"))+
  facet_wrap(DrugName~ Gene,nrow=2, scales = "free") +
  theme_classic(base_size = 12) +
  labs(
    title = "Gene essentiality in resistant vs sensitive cancer cell lines",
    y = "Gene essentiality (-CERES score)", #"CERES score (lower = more essential)",
    x = "Cell line drug resistance (AUC)"
  ) +
  scale_x_continuous(
    trans = scales::logit_trans(),
    breaks = (c(0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.975,0.99,0.995)),
    #labels = scales::percent_format(accuracy = 1)
  ) +
  theme(
    aspect.ratio=1,legend.position = "none",
  )


##### Analysis of Ribo-resistant lines
# Q) Does ribociclib-resistance change the relationship between gene essentiality and afatinib response?
# Define the ribo-resistant cell lines
rib_res_lines <- unique(plot_data[DrugName=="Ribociclib"][Sensitivity=="Resistant"]$CellLineName)
cell_lines_used_for_prediction <- unique(plot_data[`Drug Name`=="Ribociclib"]$CellLineName)
# Select afatinib response data and focal genes
rib_resist_dd<- plot_data[DrugName=="Afatinib"][Gene%in% genes_of_interest][CellLineName%in%cell_lines_used_for_prediction]
# annotate cell lines as resistant or sensitive
rib_resist_dd[,ribo_R:=CellLineName%in%rib_res_lines]
rib_resist_dd[, scalelogitAUC:=scale(logitAUC)]

# Fit linear model
rewir_res_raw <- rbindlist( lapply( 1:length(unique(rib_resist_dd$Gene)), function(i){
  #Does ribociclib resistance change the relationship between baseline gene-g dependency and afatinib response?
  g<-unique(rib_resist_dd$Gene)[i]
  lm_rewire_g<- lm( I(-CERES)~ scalelogitAUC * ribo_R , data=rib_resist_dd[Gene==g])
  afatesent_stats <- data.table(coef(summary(lm_rewire_g)),keep.rownames = T)
  afatesent_stats[, Gene:=g]
  return(afatesent_stats)
} ) )
rewir_res <- rewir_res_raw[rn!="(Intercept)"]

setnames(rewir_res_raw, old=c("rn","Std. Error", "t value", "Pr(>|t|)"), new=c("Parameter","Std.Error", "t-value", "p-value"))
rewir_res_raw[,ParameterName:= "Baseline"]
rewir_res_raw[Parameter=="ribo_RTRUE",ParameterName:= "Ribociclib resistance"]
rewir_res_raw[Parameter=="scalelogitAUC",ParameterName:= "Afatinib resistance"]
rewir_res_raw[Parameter=="scalelogitAUC:ribo_RTRUE",ParameterName:= "Afatinib resistance modulated by Ribociclib resistance"]
rewir_res_raw <- rewir_res_raw%>%select(c("Gene","ParameterName", "Estimate", "Std.Error", "t-value", "p-value"))
write.csv(rewir_res_raw, file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/Lab U01/Table S5 drug resistance effect on gene essentiality.csv")
# Keep only main effects and interaction terms
plot_dataheat <- rewir_res[rn %in% c("ribo_RTRUE","scalelogitAUC", "scalelogitAUC:ribo_RTRUE")]

# Add a more readable column for effect type
plot_dataheat[, effect := "Afatinib\nresistance"]
plot_dataheat[rn == "ribo_RTRUE", effect :="Ribociclib\nresistance"]
plot_dataheat[rn == "scalelogitAUC:ribo_RTRUE", effect :="Afatinib resistance \n modulated by \n  Ribociclib resistance"]
plot_dataheat[, effect := factor(effect, levels=c( "Ribociclib\nresistance","Afatinib\nresistance", "Afatinib resistance \n modulated by \n  Ribociclib resistance"))]

#plot_data[, effect := ifelse(rn == "ribo_RTRUE", "Ribociclib\nresistance", "Afatinib resistance \n modulated by \n  Ribociclib resistance")]
#plot_data[, effect := factor(effect, levels=c( "Ribociclib\nresistance", "Afatinib resistance \n modulated by \n  Ribociclib resistance"))]
# Optional: mark significance
plot_dataheat[, sig := ""]
plot_dataheat[`Pr(>|t|)` < 0.05, sig := "*"]
plot_dataheat[`Pr(>|t|)`  < 0.01, sig := "**"]
plot_dataheat[`Pr(>|t|)`  < 0.001, sig := "***"]

# Ensure genes are ordered by effect size (for nicer plotting)
#rev(plot_data[rn=="ribo_RTRUE"][order(Estimate)]$Gene)
plot_dataheat[, Gene := factor(Gene, levels=(plot_dataheat[rn=="ribo_RTRUE"][order(Estimate)]$Gene) )]#     unique(Gene[order(Estimate)]))]
ggplot(plot_dataheat, aes(x=effect, y=Gene, fill=Estimate)) +
  geom_tile(color="white") +
  #geom_text(aes(label=sig), color="slategrey", size=4) +
  geom_text(aes(label=sig), color="black", size=4) +
  scale_fill_gradient2(low="darkblue", mid="white", high="darkred", midpoint=0) +
  theme_classic(base_size=12) +
  labs(
    x="Drug-resistance effect",
    y="Gene",
    fill="Gene \n essentiality \n change"#Effect on\n gene\nessentiality"
  ) +
  theme(
    axis.text.x = element_text(angle=90, vjust=0.5,hjust=1),
    aspect.ratio=1/5
  )+coord_flip()+
  theme(text = element_text(color = "black"))
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/LinearModelCoefficients_V2.pdf", width=12, height=4,dpi=320)
ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/LinearModelCoefficients_V2.png", width=12, height=4,dpi=320)


plot_dataheat[, label_effect := paste("Afatinib", effect)]

ggplot(plot_dataheat, aes(x=label_effect, y=Gene, fill=Estimate)) +
  geom_tile(color="white") +
  geom_text(aes(label=sig), color="black", size=4) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_classic(base_size=12) +
  labs(
    x="Drug and effect type",
    y="Gene",
    fill="Effect on\ngene essentiality"
  ) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1),
    aspect.ratio=1
  )


rewir_res[rn%in%c("scalelogitAUC:ribo_RTRUE","ribo_RTRUE")][Estimate>0][`Pr(>|t|)`< 0.05]
rewir_res[rn=="scalelogitAUC:ribo_RTRUE"][  Gene%in%c("EGFR","ERBB2","ERBB3","ERBB4") ]

glist <- rewir_res[rn%in%c("scalelogitAUC:ribo_RTRUE","ribo_RTRUE")][Estimate>0][`Pr(>|t|)`< 0.05][Gene!="CDKN1A"]$Gene


ggplot(rib_resist_dd[ ][  Gene%in%glist ],
       aes(y=-CERES, x=AUC, col=ribo_R))+#, group=(CellLineName%in%rib_res_lines)) )+
  theme_classic(base_size = 12) +
  geom_point()+
  facet_wrap(Gene~., scales = "free_y", nrow=2)+#(CellLineName%in%rib_res_lines), scales = "free_x" )+
  geom_smooth(method="lm")+
  labs(
    #title = "Gene essentiality in afatinib resistant vs sensitive cancer cell lines",
    y = "Gene essentiality (-CERES score)", #"CERES score (lower = more essential)",
    x = "Afatinib resistance (AUC)"
  ) +
  scale_x_continuous(
    trans = scales::logit_trans(),
    breaks = (c(0.2,0.4, 0.6,0.8,0.9,0.95,0.975,0.99)),
    #labels = scales::percent_format(accuracy = 1)
  ) +
  theme(
    aspect.ratio=1,legend.position = "none",
    #strip.background = element_rect(fill = "grey90")
  )
