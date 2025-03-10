# U01_CDKi_combine_ERBBi_dynamic_cell_cycle_control
This repository contains code accompanying the project: "Combination CDK4/6 and ErbB inhibition induces cell cycle arrest and apoptosis in ribociclib-resistant estrogen receptor-positive (ER+) breast cancer". 

It provides code used to analyze the dynamic transcriptional RNA sequencing data of CDK4/6 inhibitor-resistant and -sensitive ER+ breast cell cell lines during individual or combined treatment with CDK4/6i (ribociclib) or ERBBi (afatinib).


Related publications (Please consult and cite accordingly when using this repository):
Griffiths, J. I., Chen, J., Cosgrove, P. A., O’Dea, A., Sharma, P., Ma, C., ... & Bild, A. H. (2021). Serial single-cell genomics reveals convergent subclonal evolution of resistance as patients with early-stage breast cancer progress on endocrine plus CDK4/6 therapy. Nature cancer, 2(6), 658-671.https://www.nature.com/articles/s43018-021-00215-7

Griffiths, J. I., Cosgrove, P. A., Medina, E. F., Nath, A., Chen, J., Adler, F. R., ... & Bild, A. H. (2025). Cellular interactions within the immune microenvironment underpins resistance to cell cycle inhibition in breast cancers. Nature Communications, 16(1), 2132.https://www.nature.com/articles/s41467-025-56279-x

# Environment set up
The following software is required:

*R version 1.9.0 or later with core packages:
  * GSVA (version 1.48.3)
  * msigdbr (version 7.5.1)
  * openxlsx (version 4.2.5.2)
  * readxl (version 1.4.3)
  * mgcv (version 1.9-1)
  * lme4 (version 1.1.35.1)
  * lmerTest (version 3.1.3)
  * graphite (version 1.48.0)
  * ggplot2 (version .5.1)
  * data.table (version 1.14.10)
  * tidyr (version 1.3.0)
  * dplyr (version 1.1.4)
  * ggsci (version 3.0.0)
  * scales (version 1.3.0)
  * colorspace (version 2.1.0)
# Data availablity


# RNA sequencing data processing 


# Pathway activity score calculation


# Generalized additive model (GAM) of dynamic phenotype change during treatment
Dynamic changes in pathway activity during treatment were characterized in CDK4/6i-sensitive and -resistant cells using generalized additive models (GAM). Models integrated RNA data from across ribociclib-resistant and -sensitive cell states, over time and across individual and combined treatments. Temporal patterns of pathway activity were partitionaed into four components: i) pre-treatment differences between CDK4/6i-resistant and -sensitive cells, ii) initial response to treatment upon drug addition, iii) pathway reactivation/deactivation during treatment, and iv) plasticity in treatment response upon subsequent drug additions. Treatment effects of monotherapies were additively combined, yielding an additive expectation of the combination treatment effect. Deviation from the additive expectation quantified synergistic or antagonistic effects of combination therapy.

# Cancer population growth rate and trajectory analyses
To assess the ability of mono- and combination-treatments to control cancer proliferation, we fluorescently labelled cell lines to track cancer population size over time in 3D spheroids. 

The individual and combined effects of CDK4/6i and ERBBoi on the cancer growth trajectory of each cell line over time was characterized using a Gaussian process (GP) model. The model quantified treatment impacts across time, with 95% confidence intervals highlighting timepoints where cancer population dynamics diverged between treatment conditions. 

The speed of growth or shrinkage of each replicate cancer population over the experiment was quantified using relative growth rate (rgr). Statistical significance of synergy effects on cancer growth rates was tested using a linear mixed effects model (lme). This partition main treatment effects, to define the additive expectation (under the BLISS synergy model) and included an interaction term to quantify synergy by measuring deviation of combination-treated cancer growth rates from this additive expectation. The random effects component accounted for cell line specific growth rate differences.


# Western blot protein quantification 
Immunoblotting was applied to test GAM model predicted CDK4/6i (ribociclib) and ERBBi (afatinib) induced protein level and phosphorylation changes. Replicate experiments were performed and western blot bands were quantified using imageJ.

