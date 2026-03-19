rm(list=ls())
require(data.table)
require(dplyr)
require(tidyr)
require(ggplot2)
require(synergyfinder)

dd <- fread( file="/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/Lab U01/ribociclib_plus_afatinib_synergy_analysis_data.csv")
dd[,Intensity:=as.numeric(Intensity)]
dd_mean <- data.table( dd %>% group_by(Lineage, Afatinib_dose, Ribociclib_dose, Resistance, CellLine)%>%
                         summarise(Intensity=mean(Intensity)))

ggplot(dd_mean, aes(y=Intensity,x=Afatinib_dose, col=Ribociclib_dose, group=Ribociclib_dose))+
  theme_classic(base_size=16)+
  geom_point()+geom_smooth(method="gam", formula=y~s(x, k=4))+
  facet_grid(Resistance~ CellLine)+
  labs(col="Ribociclib (nM)",
       x="Afatinib (nM)", y="Intensity")+
  theme(aspect.ratio=1)


# Define synergy metric: using Loewe additivity, which evaluates whether the combination produces greater effects than expected if the drugs were acting through overlapping or similar pathways
synergy_method <- "Loewe" #"Bliss" #"ZIP" #"HSA" # "Loewe"

# Define groups for which to assess synergy  
dd_levs <- unique(data.table(dd_mean)%>%select(CellLine,Resistance))

# Perform synergy analyses 
output_full_all <- rbindlist( lapply(1:nrow(dd_levs), function(i){
  metadd <- dd_levs[i]
  cellline_name <- metadd$CellLine
  resistance_state<- metadd$Resistance
  
  # Normalize to untreated control intensity
  # Determine control intensity per cell lineage
  dd_mean[, control_Intensity:=sum(Intensity*(Afatinib_dose == 0 & Ribociclib_dose == 0)), by="Lineage" ]
  # Convert to % Inhibition
  dd_mean[, Response:=100 * (1 - Intensity / control_Intensity) , by="Lineage" ]
  
  #Format Data for synergyfinder
  df_synergy <- dd_mean[Resistance==resistance_state&CellLine==cellline_name] %>%
    mutate(
      BlockID = interaction(CellLine, Resistance),
      DrugRow = "Afatinib",
      DrugCol = "Ribociclib",
      ConcRow = Afatinib_dose,
      ConcCol = Ribociclib_dose,
      Response = Response
    ) %>%
    select(BlockID, DrugRow, DrugCol, ConcRow, ConcCol, Response) %>%
    group_by(BlockID) %>%
    mutate(ConcUnit = "nM") # label units

  # Structrure data for synergy analysis
  data_sf <- ReshapeData(
    data = df_synergy,
    data_type = "inhibition"   # or "viability" depending on normalization
  )
  
  # Perform synergy analysis
  res <- CalculateSynergy(
    data_sf,
    method = synergy_method
  )

  # Merge results, stats and predictions
  stats <- merge(
    data.table(res$ drug_pairs),
    data.table(res$ response), by="block_id")
  setnames(stats, old="Loewe_synergy",new="Loewe_effect")
  
  output_full<-merge(
    stats,
    data.table(res$synergy_scores)
    , by=c("block_id","conc1", "conc2" ))
  output_full[,CellLine:=cellline_name]
  output_full[,Resistance:=resistance_state]
  return(output_full)
})
)

# View statistics
unique(output_full_all%>%select(Resistance,CellLine,block_id,Loewe_synergy_p_value, Loewe_effect))[order(Resistance)]
output_full_all[,Loewe_synergy_p_value:=as.numeric(Loewe_synergy_p_value)]

max_abs <- max(abs(output_full_all$Loewe_synergy), na.rm = TRUE)
ggplot(output_full_all, 
       aes(x = sqrt(conc1), y = sqrt(conc2), fill = Loewe_synergy)) +
  theme_classic(base_size = 21) +
  facet_wrap(Resistance ~ CellLine, scales="free") +
  
  # Overlay actual measured dose combinations
  geom_point(pch=21,col="black",size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        limits = c(-max_abs, max_abs) ) +
  labs(x = "Afatinib (nM)",
       y = "Ribociclib (nM)",
       fill = "Loewe synergy",
       color = "Loewe synergy") +
  theme(aspect.ratio = 1)+
  scale_x_continuous(breaks=sqrt(sort(unique(output_full_all$conc1))), labels=sort(unique(output_full_all$conc1)))+
  scale_y_continuous(breaks=sqrt(sort(unique(output_full_all$conc2))), labels=sort(unique(output_full_all$conc2)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )
#ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/Synergy_2DDoseResponseCurves_ribo_afatinib.png", width=8, height=8,dpi=320)




  
  