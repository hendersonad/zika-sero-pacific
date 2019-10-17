# - - - - - - - - - - - - - - - - - - - - - - - 
# Pacific seroprevalence: Main run analysis
# Authors: Alasdair Henderson & Adam Kucharski
# https://github.com/a-henderson91/zika-sero-pacific
# - - - - - - - - - - - - - - - - - - - - - - - 

# Set up libraries --------------------------------------------------------
library('dplyr')
library('tidyverse')
library('ggplot2')
library('gridExtra')
library('epitools')
library('irr')
library('mgcv')

rm(list=ls()) # Clear workspace
setwd("~/Documents/GitHub/zika-sero-pacific/")

# load functions for plotting and tables ----------------------------------
source("Rscripts/Functions/functions_zika_seroprevalence_tables_figures.R")

# Load data ---------------------------------------------------------------
fp_data <- data.frame(read.csv("data/dset1-FrPoly-data.csv",stringsAsFactors = F)) # French Polynesia serology summary
fj_data <- data.frame(read.csv("data/dset2-fiji-MIAassay.csv",stringsAsFactors = F)) # Fiji MIA data (n=189)
fj_age_data <- data.frame(read.csv("data/dset2b-fiji-age.csv",stringsAsFactors = F)) # Fiji age distribution of participants
fj_neut_data <- data.frame(read.csv("data/dset3-fiji-neutralizationassay.csv",stringsAsFactors = F))  # Fiji neut data (n=45 with all three years for ZIKV)
fj_neut_data_extra_13_15 <- data.frame(read.csv("data/dset3A-fiji-neutralizationassay.csv",stringsAsFactors = F)) # Fiji neut data (with additional 2013/15 samples to test rise)
fp_ind_data <- data.frame(read.csv("data/dset4-FrPoly-individual.csv",stringsAsFactors = F)) # General pop French Polynesia serology in 2015
fj_all_data <- data.frame(read.csv("data/dset5-fiji-mergedassaydata.csv",stringsAsFactors = F)) # Fiji merged MIA/neut dataset (n=189)

fp_pop_data <- read_csv("data/dset-fp-pop-pyramid.csv") # French Polynesia age data

fj_pcr_data <- read_csv("data/dset-fj-pcr-data-raw.csv") # Fiji PCR and epi data
fp_pcr_data <- read_csv("data/dset-fp-pcr-data-raw.csv") # French Polynesia PCR and epi data

fj_all_data[fj_all_data$raw_age<=16,"raw_age"] <- 16
fj_all_data[fj_all_data$raw_age>=17,"raw_age"] <- 17
write_csv(fj_all_data,"data/dset5-fiji-mergedassaydata.csv") 


# Tests for association ---------------------------------------------------
# tests for table 1
#French Poly. 2014-2015
matrix(c(18,31,154,546),byrow = T,nrow=2) %>% chisq.test() 
#French Poly. schoolchildren 2014-2018
matrix(c(312, 476-312, 291, 457-291),byrow = T,nrow=2) %>% chisq.test() 

#Fiji 2013-2015
matrix(c(12, 189-12, 45, 189-45),byrow = T,nrow=2) %>% chisq.test() 
#Fiji 2015-2017
matrix(c(45, 189-45, 23, 189-23),byrow = T,nrow=2) %>% chisq.test() 

#French Poly. adults only in general population (n=48 in 2014, n=672 in 2015)
binom.test(17, 48) ## 2014 - 35.4%
binom.test(143, 672) ## 2015 - 21.3%

# Mcnemar tests for seroconversion/reversion in Fiji dataset
mcnemar.test.fiji("DENV1")
mcnemar.test.fiji("DENV2")
mcnemar.test.fiji("DENV3")
mcnemar.test.fiji("DENV4")
mcnemar.test.fiji("ZIKV")

# McNemar tests by age group
# adults
matrix(c(7, 2, 23, 90), byrow = T,nrow=2) %>% mcnemar.test() 
# adults
matrix(c(10, 5, 4, 48), byrow = T,nrow=2) %>% mcnemar.test() 

## Gather all tests
tests <- tests_for_association()
##French Poly. schoolchildren 2014-2018 DENV1 to DENV4
tests[[30]] ## DENV1
tests[[32]] ## DENV2
tests[[34]] ## DENV3
tests[[36]] ## DENV4
##French Poly. adults 2014-2015 DENV1 to DENV4
tests[[20]] ## DENV1
tests[[22]] ## DENV2
tests[[24]] ## DENV3
tests[[26]] ## DENV4


# Figure 1 - seroprevalence pattern by age --------------------------------
fig1_plot <- plot_fig1("ZIKV")
write.csv(fig1_plot$table_for_plot, "outputs/data_for_fig1.csv")

# Figure 1 - Supplement 1-4 - DENV seroprevalence pattern by age
plot_fig1("DENV1")
plot_fig1("DENV2")
plot_fig1("DENV3")
plot_fig1("DENV4")

# Figure 2 - histograms of neutralization titres in Fiji (ZIKV/DEN3) --------------------------------
# Figure 2 - Supplement 1 - histograms of neutralization titres in Fiji (DENV1/DEN2)
plot_fig2(virus1="Zs",virus2="D3s")

# Figure 2 - Supplement 2 - histograms of neutralization titres in Fiji (DENV1/DEN2)
plot_fig2(virus1="D1s",virus2="D2s")

# Figure 3 - correspondence between neutralization assay and MIA ------
plot_comparison_MIA_NT()

# Figure 4 - ZIKV titre distributions  --------------------------------
plot_ZIKV_titres()

# Figure 4 - supplement 1 - individual participant titres -
plot.individual.2017()

# Figure 4 - supplement 2 - neutralization titre correlation between DENV1-4 and ZIKV - no prior infection
plot_supp_titres(denv_prior=F)
dev.copy(pdf, "outputs/supp_fig3_denv_naive.pdf", width=6, height=5)
dev.off()

# Figure 4 - supplement 3 - neutralization titre correlation between DENV1-4 and ZIKV - at least one DENV
plot_supp_titres(denv_prior=T)
dev.copy(pdf, "outputs/supp_fig3_denv_infection.pdf", width=6, height=5)
dev.off()


# Table 1 - summary of seroprevalence -------------------------------------
# Table 2 - Age adjusted seroprevalence ----------
## Can be run for ZIKV or DENV1-4
tab1 <- tab1_create("ZIKV")
tab1
write_csv(tab1, "outputs/Table_1.csv")

# Table 3 - seroconversion and seroreversion to ZIKV (Fiji) ----------
supp_tab1 <- data.frame(supp_tab1_create())
supp_tab1
write_csv(supp_tab1, file.path("outputs/supp_tab1.csv"))

# Table 4 - Age and DENV exposure profile (French Poly.) -------------
supp_tab2_tab3 <- supp_tab2_tab3_create()
supp_tab2 <- data.frame(supp_tab2_tab3$supp_tab2)
supp_tab2
write_csv(supp_tab2, file.path("outputs/supp_tab2.csv"))

# Table 5 - Bootstrap estimated DENV seroprevalence ------------------
supp_tab3 <- data.frame(supp_tab2_tab3$supp_tab3)
supp_tab3
write_csv(supp_tab3, file.path("outputs/supp_tab3.csv"))

# Table 6 - Change in neutralization titre to DENV3 and ZIKV  --------
supp_tab4 <- supp_tab4_create(virus1="Zs",virus2="D3s")
supp_tab4
write_csv(supp_tab4, file.path("outputs/supp_tab4.csv"))





