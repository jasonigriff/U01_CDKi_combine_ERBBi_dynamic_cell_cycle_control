rm(list=ls())
require(data.table)
require(dplyr)
require(ggplot2)
require(tidyr)
require(drc)
data_loc<-"/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/data/Lab U01/"
data.file<- "Dose Response Curves Raw Data.csv"

dd <- fread(file.path(data_loc,data.file))
# calculate % viability
dd[,mean_control:=sum(intensity*(dose==0))/sum(dose==0) , by=c("drug","lineage")]
dd[,viability_norm :=intensity / mean_control * 100]

dd[,log1pdose:=log10(1+dose)]
dd[,cell_line:=factor(cell_line)]
dd[,drug := factor(drug)]
dd[,resistance := factor(resistance)]
dd[, curve := interaction(drug, lineage)]

dd_fit <- dd[dose > 0]

comb_drm <- drm(viability_norm ~ dose,
                curveid = curve,
                data = dd_fit,
                fct = LL.4()
                )

ED(comb_drm, 50, interval = "delta")
compParm(comb_drm, "e")
plot(comb_drm, log = "x", col = 1:8)

newdata <- data.table(
  merge(
  expand.grid(
    dose = exp(seq(log(min(dd_fit$dose)),
                   log(max(dd_fit$dose)),
                   length = 100)),
    lineage = unique(dd_fit$lineage),
    drug = factor(unique(dd_fit$drug))
  ),
  unique( dd_fit%>%select(lineage,cell_line,resistance,drug,curve))
  )
)
preds<- data.table(predict(comb_drm, newdata=newdata,se.fit =T))
preds$Prediction
newdata[,pred := preds$Prediction]
newdata[,ucl := preds$Prediction+1.96*preds$SE]
newdata[,lcl := preds$Prediction-1.96*preds$SE]
newdata[lineage=="CAMA-1 resistant"][drug=="afatinib"]

ggplot(dd_fit,
       aes(x = dose, y = viability_norm,  group=(curve), color = drug)) +
  theme_classic(base_size=16)+
  geom_point() +
  geom_line(data = newdata,
            aes(y = pred)) +
  scale_x_log10() +
  facet_wrap(.~lineage) +
  theme(aspect.ratio=1)+
  labs(y="Viability (% relative to DMSO control)",
       x="Dose (nM)",
       color="Drug"
       )
#ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/DoseResponseCurves.png", width=7, height=7,dpi=320)

ggplot(dd_fit,
       aes(x = dose, y = viability_norm,  group=(curve), color = resistance)) +
  theme_classic(base_size=16)+
  geom_point() +
  geom_line(data = newdata,
            aes(y = pred)) +
  scale_x_log10() +
  facet_grid(drug~gsub("_", "-", cell_line )) +
  theme(aspect.ratio=1)+
  labs(y="Viability (% relative to DMSO control)",
       x="Dose (nM)",
       color="Resistance"
  )
#ggsave("/Users/jgriffiths/Jason Griffiths Dropbox/jason griffiths/Cancer_pheno_evo/images and presentations/Lab U01/DoseResponseCurvesResistanceConstrast.png", width=7, height=7,dpi=320)

compParm(comb_drm, "e")
summary(comb_drm)

test_resist_effect <- function(dat){
  
  # full model: separate curves per drug
  m_full <- drm(viability_norm ~ dose,
                curveid = resistance,
                data = dat,
                fct = LL.4())
  
  # reduced model: same curve for both drugs
  m_reduced <- drm(viability_norm ~ dose,
                   data = dat,
                   fct = LL.4())
  x<-anova(m_reduced, m_full)
  x<-cbind(curveresist=unique(dat$curveresist), as.data.table(x))
}
dd_fit[,curveresist:=interaction(drug,cell_line)]
results <- rbindlist(lapply(split(dd_fit, dd_fit$curveresist), test_resist_effect))
na.omit(results)

test_resist_effect2 <- function(dat){
  
  # full model: separate curves per drug
  m_full <- drm(viability_norm ~ dose,
                curveid = resistance,
                data = dat,
                fct = LL.4())
  cbind(curveresist=unique(dat$curveresist),data.table(summary(m_full)$coef,keep.rownames=T))
}
results2 <- rbindlist(lapply(split(dd_fit, dd_fit$curveresist), test_resist_effect2))
results2[curveresist=="afatinib.CAMA_1"]



