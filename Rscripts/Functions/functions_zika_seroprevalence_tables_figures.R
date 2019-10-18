# - - - - - - - - - - - - - - - - - - - - - - - 
# Pacific seroprevalence: Functions to create tables and figures
# Author: Alasdair Henderson
# github.com/a-henderson91/zika-sero-pacific.git
# - - - - - - - - - - - - - - - - - - - - - - - 


# Define virus colours
viruses_1 <- c("DENV1","DENV2","DENV3","DENV4","ZIKV")
viruses_epi <- c("DENV1","DENV2","DENV3","DENV4","ZIKV","Unknown")
lab_names <- c("DENV-1","DENV-2","DENV-3","DENV-4","ZIKV")

pick_sera <- c("D1s13","D1s15","D1s17","D2s13","D2s15","D2s17","D3s13","D3s15","D3s17","D4s13","D4s15","D4s17","Zs13","Zs15","Zs17")
pick_5 <- c("Zs13","Zs15","Zs17")

col_def <- list(col1=rgb(0.8,0.6,0),col2=rgb(0,0.8,0.8),col3=rgb(0.1,0.4,0.1),col4=rgb(1,0.4,1),col5=rgb(0.2,0,0.8),"grey")
col_def_F <- list(col1=rgb(0.8,0.6,0,0.5),col2=rgb(0,0.8,0.8,0.5),col3=rgb(0.1,0.4,0.1,0.5),col4=rgb(1,0.4,1,0.5),col5=rgb(0.2,0,0.8,0.5),"grey")

b_test <- function(x,n){
  if(length(x)==1){
    test_r <- binom.test(x,n)
    test_out <- 100*c(x/n,test_r$conf.int[1:2]) %>% signif(3)
    text_out <- paste0(test_out[1],"% (",test_out[2],"-",test_out[3],"%)")
  }else{
    text_out <- NULL
    for(ii in 1:length(x)){
      test_r <- binom.test(x[ii],n[ii])
      test_out <- 100*c(test_r$estimate,test_r$conf.int[1:2]) %>% signif(3)
      text_out <- c(text_out,paste0(test_out[1],"% (",test_out[2],"-",test_out[3],"%)"))
    }
  }
  
  text_out
}


# Table 1 - summary of seroprevalence -------------------------------------
tab1_create <- function(virus){
  ### FIJI 
  ## Survey dates and population
  table1_fj <- data.frame(Date=rep(NA,3), Country=rep(NA,3), Population=rep(NA,3))
  table1_fj$Date <- c("Oct-Nov 2013", "Nov 2015", "Jun-2017")
  table1_fj$Country <- rep("Central Division, Fiji", 3)
  table1_fj$Population <- rep("General", 3)
  ## Age range
  age_2013 <- fj_age_data$age[!is.na(fj_age_data$age)]
  age_2015 <- age_2013 + as.numeric((as.Date("2015-11-01") - as.Date("2013-11-01"))/365)
  age_2017 <- round(age_2013 + as.numeric((as.Date("2017-06-01") - as.Date("2013-11-01"))/365))
  table1_fj$AgeRange <- rbind(paste0(min(age_2013),"-",max(age_2013)," (",median(age_2013),")"),
                              paste0(min(age_2015),"-",max(age_2015)," (",median(age_2015),")"),
                              paste0(min(age_2017),"-",max(age_2017)," (",median(age_2017),")"))
  ## Seroprevalence estimates
  # Summarise number positive and sample size
  fj_positive <- fj_data %>% 
    dplyr::select(paste0(virus, 13),paste0(virus, 15),paste0(virus, 17)) %>%
    summarise_all(~sum(.)) %>%
    as.numeric()
  fj_n <- fj_data %>% summarise(n=n()) %>% as.numeric()
  table1_fj$total_tested_n <- rbind(
    paste0(fj_positive[1],"/",fj_n),
    paste0(fj_positive[2],"/",fj_n),
    paste0(fj_positive[3],"/",fj_n)
  )
  
  # Binomial test in 3 samples
  CI.calc13 <- binom.test(x=fj_positive[1], n=fj_n)        
  CI.calc15 <- binom.test(x=fj_positive[2], n=fj_n)        
  CI.calc17 <- binom.test(x=fj_positive[3], n=fj_n)        
  
  # present results in table format
  fj_seroprevalence <- rbind(
        paste0(m=signif(CI.calc13$estimate*100,2)," [", lci=signif(CI.calc13$conf.int[1]*100,2),"-",uci=signif(CI.calc13$conf.int[2]*100,2),"]"),
        paste0(m=signif(CI.calc15$estimate*100,2)," [", lci=signif(CI.calc15$conf.int[1]*100,2),"-",uci=signif(CI.calc15$conf.int[2]*100,2),"]"),
        paste0(m=signif(CI.calc17$estimate*100,2)," [", lci=signif(CI.calc17$conf.int[1]*100,2),"-",uci=signif(CI.calc17$conf.int[2]*100,2),"]"))
  table1_fj <- cbind(table1_fj, fj_seroprevalence)
  names(table1_fj)[6] <- "Seroprevalence"
  
  ### FRENCH POLYNESIA - results for blood donors previously published. ref - 
  ## Survey dates and population 
  table1_fp <- data.frame(Date=rep(NA,5), Country=rep(NA,5), Population=rep(NA,5))
  table1_fp$Date <- c("Jul 2011- Oct 2013","Feb-Mar 2014", "Sep-Nov 2015", "May-Jun 2014", "Jun 2018")
  table1_fp$Country <- rep("Society Islands, French Polynesia", 5)
  table1_fp$Population <- c("Blood donors", rep("General", 2), rep("School children", 2))
  table1_fp$AgeRange <- c("18-75 (36)", "13-77 (47)", "4-88 (43)", "6-16 (11)", "6-16 (11)")
  
  ## Seroprevalence estimates
  # Summarise number positive and sample size
  fp_positive <- fp_data %>% 
    dplyr::select("Age", "pop", "year", paste0(virus, "P"),paste0(virus, "N")) %>%
    group_by(year,pop) %>%
    summarise_at(.funs=~sum(.), .vars=c(paste0(virus, "P"), paste0(virus, "N"))) %>%
    mutate(sero=signif(get(paste0(virus, "P"))/(get(paste0(virus, "P"))+get(paste0(virus, "N"))),2)) 
  table1_fp$total_tested_n <- rbind(
    "5/593",
    paste0(fp_positive[1,3],"/",fp_positive[1,3]+fp_positive[1,4]),
    paste0(fp_positive[3,3],"/",fp_positive[3,3]+fp_positive[3,4]),
    paste0(fp_positive[2,3],"/",fp_positive[2,3]+fp_positive[2,4]),
    paste0(fp_positive[4,3],"/",fp_positive[4,3]+fp_positive[4,4])
  )
  
  # Age adjusted seroprevalence for French Polynesia
  fp_ind_data_over_x <- fp_ind_data %>% filter(Age>=17)
  fp_pop_data_over_X<- fp_pop_data #%>% filter(age>=20)
  fp_weights <- fp_pop_data_over_X$all/sum(fp_pop_data_over_X$all)
  virus_pick <- "ZIKV"
  
  # DENV1
  fp_serology_weighted <- NULL
  for(ii in 1:(length(fp_pop_data_over_X$all)-1)){
    fp_2014 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==14) %>% select(DENV1))$DENV1
    fp_2015 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==15) %>% select(DENV1))$DENV1
    fp_serology_weighted <- rbind(fp_serology_weighted,c(sum(fp_2014),length(fp_2014),sum(fp_2015),length(fp_2015)))
  }
  
  fp_serology_weighted <- data.frame(fp_serology_weighted);names(fp_serology_weighted) <- c("p2014","n2014","p2015","n2015")
  
  # Age adjustment
  nonzero2014 <- fp_serology_weighted$n2014>0; nonzero2015 <- fp_serology_weighted$n2015>0
  fp_weights2014 <- head(fp_weights,-1)[nonzero2014]/sum(head(fp_weights,-1)[nonzero2014]); fp_weights2015 <- head(fp_weights,-1)[nonzero2015]/sum(head(fp_weights,-1)[nonzero2015])
  sum(fp_weights2014*fp_serology_weighted$p2014[nonzero2014]/fp_serology_weighted$n2014[nonzero2014]); sum(fp_weights2015*fp_serology_weighted$p2015[nonzero2015]/fp_serology_weighted$n2015[nonzero2015]) 
  output14 <- round(100*ageadjust.direct(fp_serology_weighted$p2014[nonzero2014],fp_serology_weighted$n2014[nonzero2014],stdpop = fp_weights2014),1)
  output15 <- round(100*ageadjust.direct(fp_serology_weighted$p2015[nonzero2015],fp_serology_weighted$n2015[nonzero2015],stdpop = fp_weights2015),1)
  conf14 <- signif(as.numeric(100*binom.test(round(as.numeric(output14[1])*sum(fp_serology_weighted$n2014)/100), sum(fp_serology_weighted$n2014))$conf),3)
  conf15 <- signif(as.numeric(100*binom.test(round(as.numeric(output15[1])*sum(fp_serology_weighted$n2015)/100), sum(fp_serology_weighted$n2015))$conf),3)
  write_csv(data.frame(rbind(c(conf14,output14),c(conf15,output15))),"outputs/age_adjusted_DENV1.csv")
  
  # DENV2
  fp_serology_weighted <- NULL
  for(ii in 1:(length(fp_pop_data_over_X$all)-1)){
    fp_2014 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==14) %>% select(DENV2))$DENV2
    fp_2015 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==15) %>% select(DENV2))$DENV2
    fp_serology_weighted <- rbind(fp_serology_weighted,c(sum(fp_2014),length(fp_2014),sum(fp_2015),length(fp_2015)))
  }
  
  fp_serology_weighted <- data.frame(fp_serology_weighted);names(fp_serology_weighted) <- c("p2014","n2014","p2015","n2015")
  
  # Age adjustment
  nonzero2014 <- fp_serology_weighted$n2014>0; nonzero2015 <- fp_serology_weighted$n2015>0
  fp_weights2014 <- head(fp_weights,-1)[nonzero2014]/sum(head(fp_weights,-1)[nonzero2014]); fp_weights2015 <- head(fp_weights,-1)[nonzero2015]/sum(head(fp_weights,-1)[nonzero2015])
  sum(fp_weights2014*fp_serology_weighted$p2014[nonzero2014]/fp_serology_weighted$n2014[nonzero2014]); sum(fp_weights2015*fp_serology_weighted$p2015[nonzero2015]/fp_serology_weighted$n2015[nonzero2015]) 
  output14 <- round(100*ageadjust.direct(fp_serology_weighted$p2014[nonzero2014],fp_serology_weighted$n2014[nonzero2014],stdpop = fp_weights2014),1)
  output15 <- round(100*ageadjust.direct(fp_serology_weighted$p2015[nonzero2015],fp_serology_weighted$n2015[nonzero2015],stdpop = fp_weights2015),1)
  conf14 <- signif(as.numeric(100*binom.test(round(as.numeric(output14[1])*sum(fp_serology_weighted$n2014)/100), sum(fp_serology_weighted$n2014))$conf),3)
  conf15 <- signif(as.numeric(100*binom.test(round(as.numeric(output15[1])*sum(fp_serology_weighted$n2015)/100), sum(fp_serology_weighted$n2015))$conf),3)
  write_csv(data.frame(rbind(c(conf14,output14),c(conf15,output15))),"outputs/age_adjusted_DENV2.csv")
  
  # DENV3
  fp_serology_weighted <- NULL
  for(ii in 1:(length(fp_pop_data_over_X$all)-1)){
    fp_2014 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==14) %>% select(DENV3))$DENV3
    fp_2015 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==15) %>% select(DENV3))$DENV3
    fp_serology_weighted <- rbind(fp_serology_weighted,c(sum(fp_2014),length(fp_2014),sum(fp_2015),length(fp_2015)))
  }
  
  fp_serology_weighted <- data.frame(fp_serology_weighted);names(fp_serology_weighted) <- c("p2014","n2014","p2015","n2015")
  
  # Age adjustment
  nonzero2014 <- fp_serology_weighted$n2014>0; nonzero2015 <- fp_serology_weighted$n2015>0
  fp_weights2014 <- head(fp_weights,-1)[nonzero2014]/sum(head(fp_weights,-1)[nonzero2014]); fp_weights2015 <- head(fp_weights,-1)[nonzero2015]/sum(head(fp_weights,-1)[nonzero2015])
  sum(fp_weights2014*fp_serology_weighted$p2014[nonzero2014]/fp_serology_weighted$n2014[nonzero2014]); sum(fp_weights2015*fp_serology_weighted$p2015[nonzero2015]/fp_serology_weighted$n2015[nonzero2015]) 
  output14 <- round(100*ageadjust.direct(fp_serology_weighted$p2014[nonzero2014],fp_serology_weighted$n2014[nonzero2014],stdpop = fp_weights2014),1)
  output15 <- round(100*ageadjust.direct(fp_serology_weighted$p2015[nonzero2015],fp_serology_weighted$n2015[nonzero2015],stdpop = fp_weights2015),1)
  conf14 <- signif(as.numeric(100*binom.test(round(as.numeric(output14[1])*sum(fp_serology_weighted$n2014)/100), sum(fp_serology_weighted$n2014))$conf),3)
  conf15 <- signif(as.numeric(100*binom.test(round(as.numeric(output15[1])*sum(fp_serology_weighted$n2015)/100), sum(fp_serology_weighted$n2015))$conf),3)
  write_csv(data.frame(rbind(c(conf14,output14),c(conf15,output15))),"outputs/age_adjusted_DENV3.csv")
  
  # DENV4
  fp_serology_weighted <- NULL
  for(ii in 1:(length(fp_pop_data_over_X$all)-1)){
    fp_2014 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==14) %>% select(DENV4))$DENV4
    fp_2015 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==15) %>% select(DENV4))$DENV4
    fp_serology_weighted <- rbind(fp_serology_weighted,c(sum(fp_2014),length(fp_2014),sum(fp_2015),length(fp_2015)))
  }
  
  fp_serology_weighted <- data.frame(fp_serology_weighted);names(fp_serology_weighted) <- c("p2014","n2014","p2015","n2015")
  
  # Age adjustment
  nonzero2014 <- fp_serology_weighted$n2014>0; nonzero2015 <- fp_serology_weighted$n2015>0
  fp_weights2014 <- head(fp_weights,-1)[nonzero2014]/sum(head(fp_weights,-1)[nonzero2014]); fp_weights2015 <- head(fp_weights,-1)[nonzero2015]/sum(head(fp_weights,-1)[nonzero2015])
  sum(fp_weights2014*fp_serology_weighted$p2014[nonzero2014]/fp_serology_weighted$n2014[nonzero2014]); sum(fp_weights2015*fp_serology_weighted$p2015[nonzero2015]/fp_serology_weighted$n2015[nonzero2015]) 
  output14 <- round(100*ageadjust.direct(fp_serology_weighted$p2014[nonzero2014],fp_serology_weighted$n2014[nonzero2014],stdpop = fp_weights2014),1)
  output15 <- round(100*ageadjust.direct(fp_serology_weighted$p2015[nonzero2015],fp_serology_weighted$n2015[nonzero2015],stdpop = fp_weights2015),1)
  conf14 <- signif(as.numeric(100*binom.test(round(as.numeric(output14[1])*sum(fp_serology_weighted$n2014)/100), sum(fp_serology_weighted$n2014))$conf),3)
  conf15 <- signif(as.numeric(100*binom.test(round(as.numeric(output15[1])*sum(fp_serology_weighted$n2015)/100), sum(fp_serology_weighted$n2015))$conf),3) 
  write_csv(data.frame(rbind(c(conf14,output14),c(conf15,output15))),"outputs/age_adjusted_DENV4.csv")
  
  # ZIKV
  fp_serology_weighted <- NULL
  for(ii in 1:(length(fp_pop_data_over_X$all)-1)){
    fp_2014 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==14) %>% select(ZIKV))$ZIKV
    fp_2015 <- (fp_ind_data_over_x %>% filter(Age>= fp_pop_data_over_X$age[ii], Age < fp_pop_data_over_X$age[ii+1], year==15) %>% select(ZIKV))$ZIKV
    fp_serology_weighted <- rbind(fp_serology_weighted,c(sum(fp_2014),length(fp_2014),sum(fp_2015),length(fp_2015)))
  }
  
  fp_serology_weighted <- data.frame(fp_serology_weighted);names(fp_serology_weighted) <- c("p2014","n2014","p2015","n2015")
  
  # Age adjustment
  nonzero2014 <- fp_serology_weighted$n2014>0; nonzero2015 <- fp_serology_weighted$n2015>0
  fp_weights2014 <- head(fp_weights,-1)[nonzero2014]/sum(head(fp_weights,-1)[nonzero2014]); fp_weights2015 <- head(fp_weights,-1)[nonzero2015]/sum(head(fp_weights,-1)[nonzero2015])
  sum(fp_weights2014*fp_serology_weighted$p2014[nonzero2014]/fp_serology_weighted$n2014[nonzero2014]); sum(fp_weights2015*fp_serology_weighted$p2015[nonzero2015]/fp_serology_weighted$n2015[nonzero2015]) 
  output14 <- round(100*ageadjust.direct(fp_serology_weighted$p2014[nonzero2014],fp_serology_weighted$n2014[nonzero2014],stdpop = fp_weights2014),1)
  output15 <- round(100*ageadjust.direct(fp_serology_weighted$p2015[nonzero2015],fp_serology_weighted$n2015[nonzero2015],stdpop = fp_weights2015),1)
  conf14 <- signif(as.numeric(100*binom.test(round(as.numeric(output14[1])*sum(fp_serology_weighted$n2014)/100), sum(fp_serology_weighted$n2014))$conf),3)
  conf15 <- signif(as.numeric(100*binom.test(round(as.numeric(output15[1])*sum(fp_serology_weighted$n2015)/100), sum(fp_serology_weighted$n2015))$conf),3)
  write_csv(data.frame(rbind(c(conf14,output14),c(conf15,output15))),"outputs/age_adjusted_ZIKV.csv")

  # Descriptive table
  table_age <- cbind(head(fp_pop_data_over_X,-1),fp_serology_weighted$n2014,fp_serology_weighted$n2015)
  names(table_age) <- c("age","population","2014 survey","2015 survey")
  write_csv(table_age,"outputs/age_table.csv")
  
  # CIs were calculated taking into account the cluster sampling design1 and using the Fisher exact test. ref - https://wwwnc.cdc.gov/eid/article/23/4/16-1549_article
  fp_seroprevalence <- rbind(
    paste0("0.8 [0.3-2.0]"),
    paste0(signif(fp_positive[1,3]/(fp_positive[1,3]+fp_positive[1,4])*100,2)," [26-47]"),
    paste0(signif(fp_positive[3,3]/(fp_positive[3,3]+fp_positive[3,4])*100,2)," [16-28]"),
    paste0(signif(fp_positive[2,3]/(fp_positive[2,3]+fp_positive[2,4])*100,2)," [61-70]"),
    paste0(signif(fp_positive[4,3]/(fp_positive[4,3]+fp_positive[4,4])*100,2)," [58-69]")
  )
  table1_fp <- cbind(table1_fp,fp_seroprevalence)
  names(table1_fp)[6] <- "Seroprevalence"
  
  full.table.1 <- rbind(table1_fp, table1_fj)  
  full.table.1
}

# Figure 1 - seroprevalence pattern by age --------------------------------
plot_fig1 <- function(virus){
  ## FRENCH POLYNESIA
  fp_summary <- fp_data %>% select(year, Age, pop, paste0(virus, "P"),paste0(virus, "N")) %>%
    mutate(n=get(paste0(virus, "P"))+get(paste0(virus, "N"))) 
  ## Binomial calculations
  # children (14,18)
  ci_calc1 <- binom.test(x=as.numeric(fp_summary[5,4]), n=fp_summary[5,6])
  ci_calc2 <- binom.test(x=as.numeric(fp_summary[6,4]), n=fp_summary[6,6])
  # adults (14,15,13)
  ci_calc3 <- binom.test(x=as.numeric(fp_summary[2,4]), n=fp_summary[2,6])
  ci_calc4 <- binom.test(x=as.numeric(fp_summary[4,4]), n=fp_summary[4,6])
  ci_calc5 <- binom.test(x=as.numeric(fp_summary[7,4]), n=fp_summary[7,6])
  
  
  fp_sero_plot <- rbind.data.frame(
                        cbind(age=fp_summary[5,2], m=ci_calc1$estimate, lci=ci_calc1$conf.int[1], uci=ci_calc1$conf.int[2], year=2014.4),
                        cbind(age=fp_summary[6,2], m=ci_calc2$estimate, lci=ci_calc2$conf.int[1], uci=ci_calc2$conf.int[2], year=2018.4),
                        cbind(age=fp_summary[2,2], m=ci_calc3$estimate, lci=ci_calc3$conf.int[1], uci=ci_calc3$conf.int[2], year=2014.3),
                        cbind(age=fp_summary[4,2], m=ci_calc4$estimate, lci=ci_calc4$conf.int[1], uci=ci_calc4$conf.int[2], year=2015.9),
                        cbind(age=fp_summary[7,2], m=ci_calc5$estimate, lci=ci_calc5$conf.int[1], uci=ci_calc5$conf.int[2], year=2013.75))
  fp_sero_plot$Country <- "French Polynesia"
  fp_sero_plot$first.report <- c(2013.75)
  
  ## FIJI
  fj_summary <- fj_data %>% select(age, paste0(virus, 13),paste0(virus, 15),paste0(virus, 17)) %>%
    group_by(age) %>%
    summarise_at(c(paste0(virus, 13),paste0(virus, 15),paste0(virus, 17)),  ~sum(.)) 
  fj_denominator <- fj_data %>% select(age, paste0(virus, 13),paste0(virus, 15),paste0(virus, 17)) %>%
    group_by(age) %>%
    summarise(n())
  ## Binomial calculations
  nn_children <- dim(fj_data %>% group_by(age) %>% filter(age=="<=16"))[1]
  nn_adults <- dim(fj_data %>% group_by(age) %>% filter(age==">=17"))[1]
  
  # children (13,15,17)
  ci_calc1 <- binom.test(x=as.numeric(fj_summary[1,2]), n=nn_children)
  ci_calc2 <- binom.test(x=as.numeric(fj_summary[1,3]), n=nn_children)
  ci_calc3 <- binom.test(x=as.numeric(fj_summary[1,4]), n=nn_children)
  
  # adults (13,15,17)
  ci_calc4 <- binom.test(x=as.numeric(fj_summary[2,2]), n=nn_adults)
  ci_calc5 <- binom.test(x=as.numeric(fj_summary[2,3]), n=nn_adults)
  ci_calc6 <- binom.test(x=as.numeric(fj_summary[2,4]), n=nn_adults)
  
  fj_sero_plot <- rbind.data.frame(
                        cbind(age=fj_summary[1,1], m=ci_calc1$estimate, lci=ci_calc1$conf.int[1], uci=ci_calc1$conf.int[2], year=2013.9),
                        cbind(age=fj_summary[2,1], m=ci_calc4$estimate, lci=ci_calc4$conf.int[1], uci=ci_calc4$conf.int[2], year=2013.9),
                        cbind(age=fj_summary[1,1], m=ci_calc2$estimate, lci=ci_calc2$conf.int[1], uci=ci_calc2$conf.int[2], year=2015.9),
                        cbind(age=fj_summary[2,1], m=ci_calc5$estimate, lci=ci_calc5$conf.int[1], uci=ci_calc5$conf.int[2], year=2015.9),
                        cbind(age=fj_summary[1,1], m=ci_calc3$estimate, lci=ci_calc3$conf.int[1], uci=ci_calc3$conf.int[2], year=2017.4),
                        cbind(age=fj_summary[2,1], m=ci_calc6$estimate, lci=ci_calc6$conf.int[1], uci=ci_calc6$conf.int[2], year=2017.4))
  fj_sero_plot$Country <- "Fiji"
  fj_sero_plot$first.report <- c(2015.57)
  
  full_age_summary <- rbind(fj_sero_plot, fp_sero_plot)
  rownames(full_age_summary) <- NULL
  
  full_age_summary$Country <- factor(full_age_summary$Country)
  levels(full_age_summary$Country) <- c("Fiji", "French Polynesia")
  full_age_summary$year <- as.numeric(full_age_summary$year)
  full_age_summary$m <- as.numeric(full_age_summary$m)
  full_age_summary$lci <- as.numeric(full_age_summary$lci)
  full_age_summary$uci <- as.numeric(full_age_summary$uci)
  full_age_summary <- full_age_summary %>% arrange(Country, age)
  if(substr(virus,1,4)=="DENV"){first.report.line=0}else(first.report.line=3)
  dodge <- position_dodge(width=0.1)
  if(virus=="DENV1"){tt=paste0("A - ", virus)}else 
    if(virus=="DENV2"){tt=paste0("B- ", virus)}else 
      if(virus=="DENV3"){tt=paste0("C- ", virus)}else 
        if(virus=="DENV4"){tt=paste0("D- ", virus)}else{tt=""}
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Create new plot with epi curves

  par(mfrow=c(4,1),mar=c(3,3,1,3),las=0,mgp=c(2,0.7,0))
  col1b <- rgb(1,0,0)
  col1Z <- rgb(1,0.5,0)
  wk_dist <- 0.5*7/365
  
  if(virus=="ZIKV"){ymax=0.84}else{ymax=1.05}
  
  # - - - 
  # Fiji serological data
  
  fiji_c <- full_age_summary %>% filter(Country=="Fiji",age=="<=16") %>% arrange(year)
  fiji_a <- full_age_summary %>% filter(Country=="Fiji",age==">=17") %>% arrange(year)
  
  plot(fiji_c$year-wk_dist,fiji_c$m,bty="l",xlim=c(2013,2019),ylim=c(0,ymax),col=col1b,type="l",yaxs="i",xlab="",ylab="Seropositive by MIA",main=paste0("Fiji",if(virus!="ZIKV"){paste0(" - ",virus)}))
  grid(ny = NULL, nx = NA, col = "lightgray",lty=1)
  points(fiji_c$year-wk_dist,fiji_c$m,col=col1b,pch=19)
  for(ii in 1:length(fiji_c$m)){ lines(c(fiji_c[ii,"year"],fiji_c[ii,"year"])-wk_dist,c(fiji_c[ii,"lci"],fiji_c[ii,"uci"]),col=col1b)}
  lines(fiji_a$year+wk_dist,fiji_a$m,col=col1Z)
  points(fiji_a$year+wk_dist,fiji_a$m,col=col1Z,pch=19)
  for(ii in 1:length(fiji_a$m)){ lines(c(fiji_a[ii,"year"],fiji_a[ii,"year"])+wk_dist,c(fiji_a[ii,"lci"],fiji_a[ii,"uci"]),col=col1Z)}
  title(main=LETTERS[1],adj=0)
  
  lines(c(rep(fj_sero_plot$first.report[1],2)),c(0,0.8),lty=2)
  
  text(labels = "age <=16", x = max(fiji_c$year) + 0.1,y= tail(fiji_c$m,1) ,col=col1b,adj=0)
  text(labels = "age >16", x = max(fiji_a$year) + 0.1,y= tail(fiji_a$m,1) ,col=col1Z,adj=0)
  
  # - - - 
  # Fiji epi and PCR data
  
  fiji_pcr <- fj_pcr_data
  fiji_pcr <- fiji_pcr %>% mutate(date_x = Year + as.numeric(substr(WEEK,5,6))/52)

  plot(fiji_pcr$date_x,fiji_pcr$PF-1e3,xlim=c(2013,2019),ylim=c(0,12),col="white",xlab="",frame.plot=F,ylab="PCR positive samples",yaxs="i")

  for(jj in 1:length(fiji_pcr$date_x)){
    
    data_pick_jj <- fiji_pcr[jj,] %>% select(date_x,viruses_1); data_pick_jj[viruses_1][viruses_1[is.na(data_pick_jj[,viruses_1])]] <- 0
    tally_v <- 0
    
    for(ii in 1:5){
      data_pick_ii <-  data_pick_jj[viruses_1][ii]
      if(ii==1){ min_y <- -1 }else{ min_y <- tally_v} # remove axes points
      if(data_pick_ii>0){
        lines(c(data_pick_jj$date_x,data_pick_jj$date_x) ,c(min_y,tally_v+data_pick_ii),col=col_def[[ii]],lwd=1.5,lend='square')
      }
      tally_v <- tally_v + data_pick_ii
    }
  }

  par(new=TRUE)
  plot(fiji_pcr$date_x,fiji_pcr$AFR-1e3,bty="l",xlim=c(2013,2019),ylim=c(0,300),type="l",xlab="",col="white",ylab="",xaxt="n",yaxt="n",yaxs="i")
  title(main=LETTERS[2],adj=0)
  
  lines(fiji_pcr$date_x,fiji_pcr$PF,col="black")
  
  for(ii in 1:length(viruses_epi)){ text(labels = lab_names[ii], x = 2018.2,y= 310-20*ii ,col=col_def[[ii]],adj=0) }
  axis(4,col="black",col.axis="black")
  mtext("Prolonged fever", side=4, line=1.7,col="black",cex=0.7,srt=90) # Label for 2nd axis

  # - - - 
  # FP serological data
  
  fp_c <- full_age_summary %>% filter(Country=="French Polynesia",age=="<=16") %>% arrange(year)
  fp_a <- full_age_summary %>% filter(Country=="French Polynesia",age==">=17") %>% arrange(year)
  
  plot(fp_c$year,fp_c$m,bty="l",xlim=c(2013,2019),ylim=c(0,ymax),col=col1b,type="l",lty=2,yaxs="i",xlab="",ylab="Seropositive by MIA",main=paste0("French Polynesia",if(virus!="ZIKV"){paste0(" - ",virus)}))
  grid(ny = NULL, nx = NA, col = "lightgray",lty=1)
  points(fp_c$year,fp_c$m,col=col1b,pch=19)
  for(ii in 1:length(fp_c$m)){ lines(c(fp_c[ii,"year"],fp_c[ii,"year"]),c(fp_c[ii,"lci"],fp_c[ii,"uci"]),col=col1b)}
  lines(fp_a$year,fp_a$m,col=col1Z,lty=2)
  points(fp_a$year,fp_a$m,col=col1Z,pch=19)
  for(ii in 1:length(fp_a$m)){ lines(c(fp_a[ii,"year"],fp_a[ii,"year"]),c(fp_a[ii,"lci"],fp_a[ii,"uci"]),col=col1Z)}
  
  lines(c(rep(fp_sero_plot$first.report[1],2)),c(0,0.8),lty=2)
  title(main=LETTERS[3],adj=0)
  
  text(labels = "age <=16", x = max(fp_c$year) + 0.1,y= tail(fp_c$m,1)  +if(virus=="DENV2"){0.05},col=col1b,adj=0)
  text(labels = "age >16", x = max(fp_a$year) + 0.1,y= tail(fp_a$m,1),col=col1Z,adj=0)
  
  # - - - 
  # FP epi data
  fp_pcr <- fp_pcr_data
  fp_pcr <- fp_pcr %>% mutate(date_x = Year + as.numeric(substr(WEEK,5,6))/52)
  
  plot(fp_pcr$date_x,fp_pcr$PF-1e3,xlim=c(2013,2019),ylim=c(0,100),col="white",xlab="",frame.plot=F,ylab="PCR positive samples",yaxs="i")
  
  for(jj in 1:length(fp_pcr$date_x)){
    
    data_pick_jj <- fp_pcr[jj,] %>% select(date_x,viruses_1); data_pick_jj[viruses_1][viruses_1[is.na(data_pick_jj[,viruses_1])]] <- 0
    tally_v <- 0
    
    for(ii in c(1,3,5)){ # As only DENV-1, DENV-3 and ZIKV circulating
      data_pick_ii <-  data_pick_jj[viruses_1][ii]
      if(ii==1){ min_y <- -1 }else{ min_y <- tally_v} # remove axes points
      if(data_pick_ii>0){
        lines(c(data_pick_jj$date_x,data_pick_jj$date_x) ,c(min_y,tally_v+data_pick_ii),col=col_def[[ii]],lwd=1.5,lend='square')
      }
      tally_v <- tally_v + data_pick_ii
    }
    
  }
  
  par(new=TRUE)
  
  plot(fp_pcr$date_x,fp_pcr$AFR-1e3,bty="l",xlim=c(2013,2019),ylim=c(0,1050),type="l",xlab="",col="white",ylab="",xaxt="n",yaxt="n",yaxs="i")
  title(main=LETTERS[4],adj=0)
  lines(fp_pcr$date_x,fp_pcr$Dengue_Syndrome,col="black")
  lines(fp_pcr$date_x,fp_pcr$Zika_Syndrome,col="black",lty=2)
  
  for(ii in 1:length(viruses_epi)){ text(labels = lab_names[ii], x = 2018.2,y= 1000-75*ii ,col=col_def[[ii]],adj=0) }
  
  axis(4,col="black",col.axis="black")
  #par(las=0)
  mtext("Dengue/Zika cases", side=4, line=1.7,col="black",cex=0.7,srt=90) # Label for 2nd axis
  #par(las=1)

  
  
  dev.copy(pdf, paste0("outputs/Fig1",virus,".pdf"), width=6, height=8)
  dev.off()
  

  return(list(plot=plot, table_for_plot=full_age_summary))
}

# Figure 2 - histograms of neutralization titres in Fiji (ZIKV/DEN3) ------
# virus1="Zs";virus2="D3s"
# virus1="D1s";virus2="D2s"
plot_fig2 <- function(virus1, virus2){
  
  virus_list <- c("D1s","D2s","D3s","D4s","Zs")
  virus1_id <- match(virus1,virus_list)
  virus2_id <- match(virus2,virus_list)
  
  fj_neut_virus1 <- fj_neut_data %>%
    mutate_at(.funs=~pmax(.,0), .vars=3:17) %>%
    drop_na(paste0(virus1,"13"), paste0(virus1,"15"), paste0(virus1,"17")) %>%
    filter(get(paste0(virus1,"13"))<2 & (get(paste0(virus1,"15"))-get(paste0(virus1,"13")) )>=2) %>%
    mutate(rise_1=get(paste0(virus1,"15"))-get(paste0(virus1,"13")),
           drop_1=get(paste0(virus1,"17"))-get(paste0(virus1,"15")))
  fj_neut_virus2 <- fj_neut_data %>%
    mutate_at(.funs=~pmax(.,0), .vars=3:17) %>%
    drop_na(paste0(virus2,"13"), paste0(virus2,"15"), paste0(virus2,"17")) %>%
    filter(get(paste0(virus2,"13"))<2 & (get(paste0(virus2,"15"))-get(paste0(virus2,"13")))>=2) %>%
    mutate(rise_2=get(paste0(virus2,"15"))-get(paste0(virus2,"13")),
           drop_2=get(paste0(virus2,"17"))-get(paste0(virus2,"15")))
  
  # Calculate seroreversions
  
  fj_neut_virus1 %>% filter(Zs17<2)
    
  # Plot rise in titre
  par(mfrow = c(2,1), mar=c(3,3,1,1),mgp=c(1.5,0.5,0))
  ii=0
  entries = c("d1","d2","d3","d4","zd")
  ent.rise = c("rise1","rise2","rise3","rise4","rise")
  lab.names = c("DENV-1","DENV-2","DENV-3","DENV-4","ZIKV")
  bttdy = seq(-5.5,9.5,1)
  label1=lab.names[virus1_id]
  label2=lab.names[virus2_id]
  
  # Plot decline in titre
  hist(c(0),xlab="Change in neutralisation titre (2015-17)",col="white",border="white",ylim=c(0,15),breaks = bttdy,main=NULL,freq=T,xaxs="i",yaxs="i")
  shiftX = 0
  ypos = 0
    hist(fj_neut_virus1$drop_1,col=col_def_F[[virus1_id]],lwd=2,add=T,freq=T,breaks=bttdy)
    hist(fj_neut_virus2$drop_2,col=col_def_F[[virus2_id]],lwd=2,add=T,freq=T,breaks=bttdy)
    ypos = ypos+1
    text(x=3,y=15-1*ypos,labels=label1,adj=0,col=col_def[[virus1_id]])
    ypos = ypos+1
    text(x=3,y=15-1*ypos,labels=label2,adj=0,col=col_def[[virus2_id]])
    lines(c(0,0),c(0,100),lty=2)
    title(main=LETTERS[1],adj=0)
  
  hist(c(0),xlab="Change in neutralisation titre (2013-15)",col="white",border="white",ylim=c(0,15),breaks = bttdy,main=NULL,freq=T,xaxs="i",yaxs="i")
  shiftX = 0
  ypos = 0
    hist(fj_neut_virus1$rise_1,col=col_def_F[[virus1_id]],lwd=2,add=T,freq=T,breaks=bttdy)
    hist(fj_neut_virus2$rise_2,col=col_def_F[[virus2_id]],lwd=2,add=T,freq=T,breaks=bttdy)
    ypos = ypos+1
    text(x=3,y=15-1*ypos,labels=label1,adj=0,col=col_def[[virus1_id]])
    ypos = ypos+1
    text(x=3,y=15-1*ypos,labels=label2,adj=0,col=col_def[[virus2_id]])
    lines(c(0,0),c(0,100),lty=2)
    title(main=LETTERS[2],adj=0)
    
  dev.copy(pdf, paste0("outputs/Fig3",virus1,"_",virus2,".pdf"), width=6, height=8)
  dev.off()
  
  # Include people initially seropositive
  
  fj_neut_virus1 <- fj_neut_data %>%
    mutate_at(.funs=~pmax(.,0), .vars=3:17) %>%
    drop_na(paste0(virus1,"13"), paste0(virus1,"15"), paste0(virus1,"17")) %>%
    filter( (get(paste0(virus1,"15"))-get(paste0(virus1,"13")))>=2) %>%
    mutate(rise_1=get(paste0(virus1,"15"))-get(paste0(virus1,"13")),
           drop_1=get(paste0(virus1,"17"))-get(paste0(virus1,"15")))
  fj_neut_virus2 <- fj_neut_data %>%
    mutate_at(.funs=~pmax(.,0), .vars=3:17) %>%
    drop_na(paste0(virus2,"13"), paste0(virus2,"15"), paste0(virus2,"17")) %>%
    filter( (get(paste0(virus2,"15"))-get(paste0(virus2,"13")))>=2) %>%
    mutate(rise_2=get(paste0(virus2,"15"))-get(paste0(virus2,"13")),
           drop_2=get(paste0(virus2,"17"))-get(paste0(virus2,"15")))
  
  # Plot rise in titre
  par(mfrow = c(2,1), mar=c(3,3,1,1),mgp=c(1.5,0.5,0))
  ii=0
  entries = c("d1","d2","d3","d4","zd")
  ent.rise = c("rise1","rise2","rise3","rise4","rise")
  lab.names = c("DENV-1","DENV-2","DENV-3","DENV-4","ZIKV")
  bttdy = seq(-5.5,9.5,1)
  label1=lab.names[virus1_id]
  label2=lab.names[virus2_id]
  
  # Plot decline in titre
  hist(c(0),xlab="Change in neutralisation titre (2015-17)",col="white",border="white",ylim=c(0,15),breaks = bttdy,main=NULL,freq=T,xaxs="i",yaxs="i")
  shiftX = 0
  ypos = 0
  hist(fj_neut_virus1$drop_1,col=col_def_F[[virus1_id]],lwd=2,add=T,freq=T,breaks=bttdy)
  hist(fj_neut_virus2$drop_2,col=col_def_F[[virus2_id]],lwd=2,add=T,freq=T,breaks=bttdy)
  ypos = ypos+1
  text(x=3,y=15-1*ypos,labels=label1,adj=0,col=col_def[[virus1_id]])
  ypos = ypos+1
  text(x=3,y=15-1*ypos,labels=label2,adj=0,col=col_def[[virus2_id]])
  lines(c(0,0),c(0,100),lty=2)
  title(main=LETTERS[1],adj=0)
  
  hist(c(0),xlab="Change in neutralisation titre (2013-15)",col="white",border="white",ylim=c(0,15),breaks = bttdy,main=NULL,freq=T,xaxs="i",yaxs="i")
  shiftX = 0
  ypos = 0
  hist(fj_neut_virus1$rise_1,col=col_def_F[[virus1_id]],lwd=2,add=T,freq=T,breaks=bttdy)
  hist(fj_neut_virus2$rise_2,col=col_def_F[[virus2_id]],lwd=2,add=T,freq=T,breaks=bttdy)
  ypos = ypos+1
  text(x=3,y=15-1*ypos,labels=label1,adj=0,col=col_def[[virus1_id]])
  ypos = ypos+1
  text(x=3,y=15-1*ypos,labels=label2,adj=0,col=col_def[[virus2_id]])
  lines(c(0,0),c(0,100),lty=2)
  title(main=LETTERS[2],adj=0)
  
  dev.copy(pdf, paste0("outputs/Fig3_all_",virus1,"_",virus2,".pdf"), width=6, height=8)
  dev.off()
    
}


# Tests for association ---------------------------------------------------
mcnemar.test.fiji <- function(virus){
  fj_test <- fj_data %>% 
    mutate(PP=as.numeric(get(paste0(virus,"15"))+get(paste0(virus,"17"))==2),
           PN=as.numeric(get(paste0(virus,"15"))==1 & get(paste0(virus,"17"))==0),
           NP=as.numeric(get(paste0(virus,"15"))==0 & get(paste0(virus,"17"))==1),
           NN=as.numeric(get(paste0(virus,"15"))+get(paste0(virus,"17"))==0)) %>%
    group_by(age) %>%
    summarise(PP=sum(PP), PN=sum(PN), NP=sum(NP), NN=sum(NN))
  ## 2015-17 Adults DENV3
  fj_ad_data <- as.numeric(fj_test[2,2:5])
  fj_ad_data
  fj_1517_mcnemar_ad <- matrix(c(fj_ad_data[1],fj_ad_data[3],fj_ad_data[2],fj_ad_data[4]), byrow = T,nrow=2) %>% mcnemar.test() 
  ## 2015-17 Children DENV3
  fj_ch_data <- as.numeric(fj_test[1,2:5])
  fj_ch_data
  fj_1517_mcnemar_sch <- matrix(c(fj_ch_data[1],fj_ch_data[3],fj_ch_data[2],fj_ch_data[4]), byrow = T,nrow=2) %>% mcnemar.test() 
  ##
  total_fj_data <- fj_ad_data+fj_ch_data
  fj_1517_mcnemar <- matrix(c(total_fj_data[1],total_fj_data[3],total_fj_data[2],total_fj_data[4]), byrow = T,nrow=2) %>% mcnemar.test()
  
  return(list("child.p.val"=fj_1517_mcnemar_sch$p.value, "adult.p.val"=fj_1517_mcnemar_ad$p.value, "overall.p.val"=fj_1517_mcnemar$p.value))
}

## all tests for association
tests_for_association <- function(){
  ## McNemar test for Zika in Fiji
  fj_test <- fj_data %>% 
    mutate(PP=as.numeric(ZIKV15+ZIKV17==2),
           PN=as.numeric(ZIKV15==1 & ZIKV17==0),
           NP=as.numeric(ZIKV17==1 & ZIKV15==0),
           NN=as.numeric(ZIKV15+ZIKV17==0)) %>%
    group_by(age) %>%
    summarise(PP=sum(PP), PN=sum(PN), NP=sum(NP), NN=sum(NN))
  ## 2015-17 Adults
  fj_ad_data <- as.numeric(fj_test[2,2:5])
  fj_ad_data
  fj_1517_mcnemar_ad <- matrix(c(fj_ad_data[1],fj_ad_data[3],fj_ad_data[2],fj_ad_data[4]), byrow = T,nrow=2) %>% mcnemar.test() 
  ## 2015-17 Children
  fj_ch_data <- as.numeric(fj_test[1,2:5])
  fj_ch_data
  fj_1517_mcnemar_sch <- matrix(c(fj_ch_data[1],fj_ch_data[3],fj_ch_data[2],fj_ch_data[4]), byrow = T,nrow=2) %>% mcnemar.test() 
  ## 2015-17 total
  total_fj_data <- fj_ad_data+fj_ch_data
  fj_1517_mcnemar <- matrix(c(total_fj_data[1],total_fj_data[3],total_fj_data[2],total_fj_data[4]), byrow = T,nrow=2) %>% mcnemar.test()
  
  ## 2015-17 total - chi squared test
  n_ad <- sum(fj_ad_data)
  n_ch <- sum(fj_ch_data)
  n_tot <- n_ad+n_ch
  pos_15_ad <- fj_ad_data[1]+fj_ad_data[2]; pos_17_ad <- fj_ad_data[1]+fj_ad_data[3]
  pos_15_ch <- fj_ch_data[1]+fj_ch_data[2]; pos_17_ch <- fj_ch_data[1]+fj_ch_data[3]
  pos_15 <- total_fj_data[1]+total_fj_data[2]; pos_17 <- total_fj_data[1]+total_fj_data[3]
  matrix(c(pos_15_ad, n_ad-pos_15_ad, pos_17_ad, n_ad-pos_17_ad), byrow = T,nrow=2) %>% chisq.test()
  matrix(c(pos_15_ch, n_ch-pos_15_ch, pos_17_ch, n_ch-pos_17_ch), byrow = T,nrow=2) %>% chisq.test()
  fj_1517_chisq <- matrix(c(pos_15, n_tot-pos_15, pos_17, n_tot-pos_17), byrow = T,nrow=2) %>% chisq.test()
  
  ##FP
  fp_chisq <- matrix(c(18,31,154,546),byrow = T,nrow=2) %>% chisq.test() ## FP data from Mai table https://wwwnc.cdc.gov/eid/article/23/4/16-1549-t1
  ##adults
  fp_chisq_ad <- matrix(c(17,48-17,143,672-143),byrow = T,nrow=2) %>% chisq.test() ## FP data from Mai table https://wwwnc.cdc.gov/eid/article/23/4/16-1549-t1
  ##children
  fp_chisq_sch <- matrix(c(311,476-311,291,457-291),byrow = T,nrow=2) %>% chisq.test() ## FP data from Mai table https://wwwnc.cdc.gov/eid/article/23/4/16-1549-t1
  
  ## French Polynesia - general - DENV 1-4- adults only
  fp_data_for_tests <- fp_data %>%
    mutate(sample_size=DENV1P+DENV1N) %>%
    filter(Age==">=17") %>%
    dplyr::select(pop, year, DENV1P,DENV2P,DENV3P,DENV4P,ZIKVP, sample_size) %>%
    group_by(pop, year) %>%
    summarise_all(~sum(.)) %>%
    filter(pop=="general")
  fpG1 <- matrix(c(fp_data_for_tests$DENV1P[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$DENV1P[1], fp_data_for_tests$DENV1P[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$DENV1P[2]), byrow = T, nrow = 2) %>% chisq.test()
  fpG2 <- matrix(c(fp_data_for_tests$DENV2P[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$DENV2P[1], fp_data_for_tests$DENV2P[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$DENV2P[2]), byrow = T, nrow = 2) %>% chisq.test()
  fpG3 <- matrix(c(fp_data_for_tests$DENV3P[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$DENV3P[1], fp_data_for_tests$DENV3P[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$DENV3P[2]), byrow = T, nrow = 2) %>% chisq.test()
  fpG4 <- matrix(c(fp_data_for_tests$DENV4P[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$DENV4P[1], fp_data_for_tests$DENV4P[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$DENV4P[2]), byrow = T, nrow = 2) %>% chisq.test()
  fpG5 <- matrix(c(fp_data_for_tests$ZIKVP[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$ZIKVP[1], fp_data_for_tests$ZIKVP[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$ZIKVP[2]), byrow = T, nrow = 2) %>% chisq.test()
  
  ## French Polynesia - schoolchildren - DENV 1-4 
  fp_data_for_tests <- fp_data %>%
    mutate(sample_size=DENV1P+DENV1N) %>%
    dplyr::select(pop, year, DENV1P,DENV2P,DENV3P,DENV4P,ZIKVP, sample_size) %>%
    group_by(pop, year) %>%
    summarise_all(~sum(.)) %>%
    filter(pop=="schoolchildren")
  fp1 <- matrix(c(fp_data_for_tests$DENV1P[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$DENV1P[1], fp_data_for_tests$DENV1P[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$DENV1P[2]), byrow = T, nrow = 2) %>% chisq.test()
  fp2 <- matrix(c(fp_data_for_tests$DENV2P[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$DENV2P[1], fp_data_for_tests$DENV2P[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$DENV2P[2]), byrow = T, nrow = 2) %>% chisq.test()
  fp3 <- matrix(c(fp_data_for_tests$DENV3P[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$DENV3P[1], fp_data_for_tests$DENV3P[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$DENV3P[2]), byrow = T, nrow = 2) %>% chisq.test()
  fp4 <- matrix(c(fp_data_for_tests$DENV4P[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$DENV4P[1], fp_data_for_tests$DENV4P[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$DENV4P[2]), byrow = T, nrow = 2) %>% chisq.test()
  fp5 <- matrix(c(fp_data_for_tests$ZIKVP[1], fp_data_for_tests$sample_size[1]-fp_data_for_tests$ZIKVP[1], fp_data_for_tests$ZIKVP[2], fp_data_for_tests$sample_size[2]-fp_data_for_tests$ZIKVP[2]), byrow = T, nrow = 2) %>% chisq.test()
  
  ## t-tests for rise and drop in ZIKV/DENV3
  virus1 <- "Zs"
  virus2 <- "D3s"
  fj_neut_virus1 <- fj_neut_data %>%
    mutate_at(.funs=~pmax(.,0), .vars=3:17) %>%
    filter(!is.na(get(paste0(virus1,"13"))) & get(paste0(virus1,"13"))<2 & (get(paste0(virus1,"15"))-get(paste0(virus1,"13")))>=2) %>%
    mutate(rise_1=get(paste0(virus1,"15"))-get(paste0(virus1,"13")),
           drop_1=get(paste0(virus1,"17"))-get(paste0(virus1,"15")))
  fj_neut_virus2 <- fj_neut_data %>%
    mutate_at(.funs=~pmax(.,0), .vars=3:17) %>%
    filter(!is.na(get(paste0(virus2,"13"))) & get(paste0(virus2,"13"))<2 & (get(paste0(virus2,"15"))-get(paste0(virus2,"13")))>=2) %>%
    mutate(rise_2=get(paste0(virus2,"15"))-get(paste0(virus2,"13")),
           drop_2=get(paste0(virus2,"17"))-get(paste0(virus2,"15")))
  test_vals_Zrise <- fj_neut_virus1$rise_1
  test_vals_Zdrop <- fj_neut_virus1$drop_1
  test_vals_Drise <- fj_neut_virus2$rise_2
  test_vals_Ddrop <- fj_neut_virus2$drop_2
  length(test_vals_Zrise);length(test_vals_Zdrop);length(test_vals_Drise);length(test_vals_Ddrop)
  ## RISE 2013-15
  t_test2 = c(t.test(test_vals_Drise)$p.value,t.test(test_vals_Zrise)$p.value,
              t.test(test_vals_Drise)$estimate,t.test(test_vals_Zrise)$estimate)
  t_test2 = matrix(t_test2,ncol=2)
  rownames(t_test2) = c("d3_rise","z_rise") 
  colnames(t_test2) = c("p","mean") 
  t_test2
  ## DROP 2015-17
  t_test1 = c(t.test(test_vals_Ddrop)$p.value,t.test(test_vals_Zdrop)$p.value,
              t.test(test_vals_Ddrop)$estimate,t.test(test_vals_Zdrop)$estimate,
              t.test(test_vals_Ddrop)$conf[1],t.test(test_vals_Zdrop)$conf[1],
              t.test(test_vals_Ddrop)$conf[2],t.test(test_vals_Zdrop)$conf[2])
  
  t_test1 = matrix(t_test1,ncol=4)
  rownames(t_test1) = c("d3_wane","z_wane")
  colnames(t_test1) = c("p","mean", "LCI", "UCI") 
  t_test1
  return(list("FJ - Z - 2015-17 - McNemars test", fj_1517_mcnemar, "FJ - Z - 2015-17 - Chi-Squared test", fj_1517_chisq,
       "FJ - Z - 2015-17 - ADULTS - McNemars test", fj_1517_mcnemar_ad, "FJ - Z - 2015-17 - CHILDREN - McNemars test", fj_1517_mcnemar_sch,
       "FP - ZIKA - chisquared",fp_chisq,"FP - ZIKA - ADULTS - chisquared",fp_chisq_ad,"FP - ZIKA - SCHOOLCHILDREN - chisquared",fp_chisq_sch,
       "FP - D1 general",fpG1,"FP - D2 general",fpG2,"FP - D3 general",fpG3,"FP - D4 general",fpG4,"FP - Z general",fpG5,
       "FP - D1 schoolchildren",fp1,"FP - D2 schoolchildren",fp2,"FP - D3 schoolchildren",fp3,"FP - D4 schoolchildren",fp4,"FP - Z schoolchildren",fp5,
       "FJ paired neutralisation titre rise 2013-15", t_test1,"FJ paired neutralisation titre drop 2015-17", t_test2)
  )
}

# Figure 3 - ZIKV titre distributions -----------------------------

plot_ZIKV_titres <- function(){

  fj_triple_data <- fj_all_data[rowSums(!is.na(fj_all_data[,pick_5]))==3,] # select entries with triple ZIKV samples
  
  bdy <- seq(-0.5,8.5,1)
  
  par(mfrow = c(1,3),mar=c(3,3,1,1),mgp=c(1.5,0.5,0))
  
  titre_13 <- (fj_triple_data %>% select(Zs13) %>% drop_na()); titre_13[titre_13==-Inf,] <- 0
  hist(titre_13$Zs13,breaks = bdy,main="2013",col=col_def_F[[5]],xlab="ZIKV log titre",yaxs="i",ylim=c(0,35))
  lines(c(1.5,1.5),c(0,100),lty=2)
  
  titre_15 <- (fj_triple_data %>% select(Zs15) %>% drop_na()); titre_15[titre_15==-Inf,] <- 0
  hist(titre_15$Zs15,breaks = bdy,main="2015",col=col_def_F[[5]],xlab="ZIKV log titre",yaxs="i",ylim=c(0,35))
  lines(c(1.5,1.5),c(0,100),lty=2)
  
  titre_17<- (fj_triple_data %>% select(Zs17) %>% drop_na()); titre_17[titre_17==-Inf,] <- 0
  hist(titre_17$Zs17,breaks = bdy,main="2017",col=col_def_F[[5]],xlab="ZIKV log titre",yaxs="i",ylim=c(0,35))
  lines(c(1.5,1.5),c(0,100),lty=2)
  
  dev.copy(pdf, "outputs/supp_fig_hist_ZIKV_titres.pdf", width=10, height=3)
  dev.off()
  

  # - - - - - - - - - - - - - - - - - 
  # Compare positivity by MIA and neut
  
  titre_13 <- fj_triple_data %>% filter(Zs13>=2) %>% select(ZIKV13) %>% drop_na()
  titre_15 <- fj_triple_data %>% filter(Zs15>=2) %>% select(ZIKV15) %>% drop_na()
  titre_17 <- fj_triple_data %>% filter(Zs17>=2) %>% select(ZIKV17) %>% drop_na()
  merge_ZIKV <- c(titre_13$ZIKV13,titre_15$ZIKV15,titre_17$ZIKV17)
  sum(merge_ZIKV)/length(merge_ZIKV) # Sensitivity
  
  titre_13 <- fj_triple_data %>% filter(Zs13<2) %>% select(ZIKV13) %>% drop_na()
  titre_15 <- fj_triple_data %>% filter(Zs15<2) %>% select(ZIKV15) %>% drop_na()
  titre_17 <- fj_triple_data %>% filter(Zs17<2) %>% select(ZIKV17) %>% drop_na()
  merge_ZIKV <- c(titre_13$ZIKV13,titre_15$ZIKV15,titre_17$ZIKV17)
  1 - sum(merge_ZIKV)/length(merge_ZIKV) # Specificity

  
  
  
}

# # Supplementary estimate - Comparison of MIA and neut datasets -----------------------------

plot_comparison_MIA_NT <- function(){
  
  # - - - - - - - - - - - - - - - - - 
  # Compare positivity by neut for those positive by MIA
  
  fj_triple_data1 <- fj_all_data[rowSums(!is.na(fj_all_data[,pick_5]))==3,] # select entries with triple ZIKV samples
  
  neut_13 <- fj_triple_data1 %>% filter(ZIKV13==1) %>% select(Zs13) %>% drop_na()
  neut_15 <- fj_triple_data1 %>% filter(ZIKV15==1) %>% select(Zs15) %>% drop_na()
  neut_17 <- fj_triple_data1 %>% filter(ZIKV17==1) %>% select(Zs17) %>% drop_na()
  merge_titre <- as.numeric (c(neut_13$Zs13,neut_15$Zs15,neut_17$Zs17)>=2)
  positive_neut <- c(sum(merge_titre),length(merge_titre),sum(merge_titre)/length(merge_titre) ) # Sensitivity
  
  neut_13 <- fj_triple_data1 %>% filter(ZIKV13==0) %>% select(Zs13) %>% drop_na()
  neut_15 <- fj_triple_data1 %>% filter(ZIKV15==0) %>% select(Zs15) %>% drop_na()
  neut_17 <- fj_triple_data1 %>% filter(ZIKV17==0) %>% select(Zs17) %>% drop_na()
  merge_titre <- as.numeric (c(neut_13$Zs13,neut_15$Zs15,neut_17$Zs17)<2)
  negative_neut <- c(sum(merge_titre),length(merge_titre),sum(merge_titre)/length(merge_titre) ) # Sensitivity
  
  MIA_by_neut <- signif(c(positive_neut,negative_neut),3)
  
  # Calculate Cohen's kappa
  neut_results <- c(fj_triple_data1$ZIKV13,fj_triple_data1$ZIKV15,fj_triple_data1$ZIKV17)
  MIA_results <- as.numeric(c(fj_triple_data1$Zs13,fj_triple_data1$Zs15,fj_triple_data1$Zs17) >=2)
  k_out_neut_MIA <- kappa2(cbind(neut_results,MIA_results), "unweighted")  # Cohen's kappa

  
  # - - - - - - - - - - - - - - - - - 
  # Compare positivity by MIA for those positive by neut
  
  fj_triple_data <- fj_all_data[rowSums(!is.na(fj_all_data[,pick_5]))==3,] # select entries with triple ZIKV samples
  
  titre_13 <- fj_triple_data %>% filter(Zs13>=2) %>% select(ZIKV13) %>% drop_na()
  titre_15 <- fj_triple_data %>% filter(Zs15>=2) %>% select(ZIKV15) %>% drop_na()
  titre_17 <- fj_triple_data %>% filter(Zs17>=2) %>% select(ZIKV17) %>% drop_na()
  merge_titre <- c(titre_13$ZIKV13,titre_15$ZIKV15,titre_17$ZIKV17)
  positive_MIA <- c(sum(merge_titre),length(merge_titre),sum(merge_titre)/length(merge_titre) ) # Sensitivity
  
  titre_13 <- fj_triple_data %>% filter(Zs13<2) %>% select(ZIKV13) %>% drop_na()
  titre_15 <- fj_triple_data %>% filter(Zs15<2) %>% select(ZIKV15) %>% drop_na()
  titre_17 <- fj_triple_data %>% filter(Zs17<2) %>% select(ZIKV17) %>% drop_na()
  merge_titre <- c(titre_13$ZIKV13,titre_15$ZIKV15,titre_17$ZIKV17)
  negative_MIA <- c(length(merge_titre) - sum(merge_titre),length(merge_titre),  1 - sum(merge_titre)/length(merge_titre) ) # Specificity
  
  neut_by_MIA <- signif(c(positive_MIA,negative_MIA),3)
  
  # Summary table of positivity
  
  summary_tab <- rbind(MIA_by_neut,neut_by_MIA)
  row.names(summary_tab) <- c("MIA_by_neut","neut_by_MIA")
  summary_tab <- data.frame(summary_tab)

  names(summary_tab) <- c("n_pos_match","n_pos_total","percent_pos_match","n_neg_match","n_neg_total","percent_neg_match")
  
  print(summary_tab)
  print(k_out_neut_MIA)
  
  triple_samples_NT_2013 <- sum(fj_triple_data$Zs13>=2) # Positive by MIA
  triple_samples_pos_MIA <- sum(fj_triple_data1$ZIKV17==1) # Positive by MIA
  triple_samples_pos_NT <- sum(fj_triple_data1$Zs17>=2) # Positive by NT
  
  # - - - 
  # Analyse differences by mean DENV
  
  fj_data_by_denv <- fj_all_data[rowSums(!is.na(fj_all_data[,pick_5]))==3,] # select entries with triple ZIKV samples
  fj_data_by_denv <- fj_data_by_denv[!is.na(fj_data_by_denv$D3s17),] # remove line with missing entry
  
  fj_data_by_denv <- fj_data_by_denv %>%
    mutate_at(.funs=~pmax(.,0), .vars=20:34) %>%
    mutate(Zs_rise=Zs15-Zs13, Zs_drop=Zs17-Zs15) %>%
    drop_na("Zs13", "Zs15", "Zs17")

  # Examine pattern of mean DENV
  dXc <- 2
  fj_data_by_denv <- fj_data_by_denv %>% mutate(total_denv = (D1s17>=dXc) + (D2s17>=dXc) + (D3s17>=dXc) + (D4s17>=dXc))
  fj_data_by_denv <- fj_data_by_denv %>% mutate(Z13_pos = Zs13>=2,Z15_pos = Zs15>=2,Z17_pos = Zs17>=2)
  fj_data_by_denv <- fj_data_by_denv %>% mutate(mean_denv = ((D1s13) + (D2s13) + (D3s13) + (D4s13))/4)
  #fj_data_by_denv <- fj_data_by_denv %>% mutate(mean_denv = ((D1s15) + (D2s15) + (D3s15) + (D4s15))/4)
  #fj_data_by_denv <- fj_data_by_denv %>% mutate(mean_denv = ((D1s17) + (D2s17) + (D3s17) + (D4s17))/4)

  col1 <- rgb(0,0,0); col1f <- rgb(0,0,0,0.4)
  
  # Plot by scenario
  par(mfrow=c(1,3),mar=c(3,3,1,1),las=0,mgp=c(2,0.7,0))
  kk_s <- 5 # Have also tested with 2-6 as sensitivity
            # gam.check() used to test for k-index<1
  
  xmax <- 4
  shift1 <- 0.006 # shift points for plotting

  for(ii in 1:6){
    # Define colours
    if(ii%% 2==0){
      col1 <- rgb(0.9,0.5,0); col1f <- rgb(1,0.5,0,0.2)
      title(main=LETTERS[ii/2],adj=0)
    }else{
      col1 <- rgb(0.5,0.5,0.5); col1f <- rgb(0,0,0,0.2)
    }
    
    # Define GAMs
    if(ii==1){
      plot(fj_data_by_denv$mean_denv,fj_data_by_denv$Z13_pos,pch=19,col="white",xlim=c(0,xmax),ylim=c(-0.03,1.03),ylab="ZIKV seropositive in 2013",yaxs="i",xlab="Mean DENV log titre in 2013")
      modelB.P <- gam(ZIKV13 ~ s(mean_denv,k=kk_s) , data = fj_data_by_denv,family = "binomial")
      points(fj_data_by_denv$mean_denv,fj_data_by_denv$ZIKV13+shift1,col=col1f,pch=19,lwd=0)
      } # Pick MIA
    if(ii==2){
      modelB.P <- gam(Z13_pos ~ s(mean_denv,k=kk_s) , data = fj_data_by_denv,family = "binomial")
      points(fj_data_by_denv$mean_denv,fj_data_by_denv$Z13_pos-shift1,col=col1f,pch=19,lwd=0)
      } # Pick NT
    if(ii==3){
      plot(fj_data_by_denv$mean_denv,fj_data_by_denv$Z13_pos,pch=19,col="white",xlim=c(0,xmax),ylim=c(-0.03,1.03),ylab="ZIKV seropositive in 2015",yaxs="i",xlab="Mean DENV log titre in 2013")
      modelB.P <- gam(ZIKV15 ~ s(mean_denv,k=kk_s) , data = fj_data_by_denv,family = "binomial")
      points(fj_data_by_denv$mean_denv,fj_data_by_denv$ZIKV15+shift1,col=col1f,pch=19,lwd=0)
      } # Pick MIA
    if(ii==4){
      modelB.P <- gam(Z15_pos ~ s(mean_denv,k=kk_s) , data = fj_data_by_denv,family = "binomial")
      points(fj_data_by_denv$mean_denv,fj_data_by_denv$Z15_pos-shift1,col=col1f,pch=19,lwd=0)
      } # Pick NT
    if(ii==5){
      plot(fj_data_by_denv$mean_denv,fj_data_by_denv$Z13_pos,pch=19,col="white",xlim=c(0,xmax),ylim=c(-0.03,1.03),ylab="ZIKV seropositive in 2017",yaxs="i",xlab="Mean DENV log titre in 2013")
      modelB.P <- gam(ZIKV17 ~ s(mean_denv,k=kk_s) , data = fj_data_by_denv,family = "binomial")
      points(fj_data_by_denv$mean_denv,fj_data_by_denv$ZIKV17+shift1,col=col1f,pch=19,lwd=0)
      } # Pick MIA
    if(ii==6){
      modelB.P <- gam(Z17_pos ~ s(mean_denv,k=kk_s) , data = fj_data_by_denv,family = "binomial")
      points(fj_data_by_denv$mean_denv,fj_data_by_denv$Z17_pos-shift1,col=col1f,pch=19,lwd=0)
      } # Pick NT
    
    # Predict prevalence
    xx_pred <- seq(0,xmax,0.1); preds <- predict(modelB.P, newdata = list(mean_denv=xx_pred), type = "link", se.fit = TRUE)
    critval <- 1.96; upperCI <- preds$fit + (critval * preds$se.fit); lowerCI <- preds$fit - (critval * preds$se.fit)
    fit <- preds$fit
    fitPlotF <- modelB.P$family$linkinv(fit); CI1plotF <- modelB.P$family$linkinv(upperCI);  CI2plotF <- modelB.P$family$linkinv(lowerCI)
    
    # Plot GAM fits
    polygon(c(xx_pred,rev(xx_pred)),c(CI1plotF,rev(CI2plotF)),col=col1f,lty=0)
    lines(xx_pred, fitPlotF ,col=col1,lwd=2)

    # Add labels
    if(ii==2){
      text(x=1,y=0.9,labels="NT",adj=0,col=rgb(0.9,0.5,0))
      text(x=1,y=0.84,labels="MIA",adj=0,col=rgb(0.5,0.5,0.5))
    }
    
  }

  dev.copy(pdf, paste0("outputs/Compare_NT_MIA.pdf"), width=7, height=2.5)
  dev.off()

  
}



# Figure 3 - supplement 1 - Plot individual-level neutralisation titres ---------------------------------------------------


plot.individual.2017 <- function(){
  
  pick.name = c("DENV1","DENV2","DENV3","DENV4","ZIKV")
  
  unique.part = fj_neut_data$id
  
  collection_dates <- c(2013.9,2015.9,2017.4)
  
  # - - - 
  # Individual level titres
  store.data = NULL
  
  # Plot results
  par(mfrow = c(9,5),mar=c(2,2,1,1),mgp=c(1.5,0.5,0))

  for(ii in 1:length(unique.part)){
    
    data.pick = fj_neut_data[ii,]

    if( sum(!is.na(data.pick[,pick_5]))==3  ){
      
      plot(c(1,3),ylim=c(0,8),col="white",xlab="",ylab="",xlim=c(2013,2018),main=paste0("Participant ",data.pick$id),xaxt="n")

      axis(1, at=c(2013,2015,2017),labels = c(2013,2015,2017),col = "black") 
      
      titre.pick = data.pick[,pick_sera]
      titre.pick[titre.pick==-Inf]=0
      
      for(jj in 1:5){
        log.titre = titre.pick[,(3*(jj-1)+1):(3*jj)]
        points(collection_dates,log.titre,col=col_def[[jj]],pch=19)
        lines(collection_dates,log.titre,col=col_def[[jj]])
      }
      
      for(ii in 1:length(pick.name)){ text(labels = lab_names[ii], x = 2013,y= 8.6-0.6*ii ,col=col_def[[ii]],adj=0,cex=1) }
      

    }
    
  }

  dev.copy(pdf, "outputs/supp_fig_individual_titres.pdf", width=10, height=17)
  dev.off()
  
}


# Figure 3 - supplement 2 - neutralization titre correlation between DENV1-4 and ZIKV -----
plot_supp_titres <- function(denv_prior=F){
  
  ## Only use seronegative to all 5 arboviruses in 2013
  if(denv_prior==F){
    fj_neut_corr <- fj_neut_data_extra_13_15 %>% mutate_at(.funs=funs(pmax(.,0)),.vars=vars(ends_with("13"))) %>%
      filter_at(vars(ends_with("13")), .vars_predicate = all_vars(. < 2))
  }
  
  ## Only use seropositive to at least one dengue virus in 2013
  if(denv_prior==T){
    fj_neut_corr <- fj_neut_data_extra_13_15 %>% mutate_at(.funs=funs(pmax(.,0)),.vars=vars(ends_with("13"))) %>%
      filter_at(vars(ends_with("13")), .vars_predicate = any_vars(. >= 2))
  }
  
  par(mfrow = c(5,5),mar=c(2.5,2.5,1,1),mgp=c(1.5,0.5,0))
  # list colnames needed
  pick.part13 = c("D1s13","D2s13","D3s13","D4s13","Zs13")
  pick.part = c("D1s15","D2s15","D3s15","D4s15","Zs15")
  col.list =c("green","orange","purple","blue","red")
  pick.name = c("DENV1","DENV2","DENV3","DENV4","ZIKV")
  
  for(ii in 1:5){
    for(jj in 1:5){
      log.titre = as.numeric(fj_neut_corr[,pick.part[ii]]) - as.numeric(fj_neut_corr[,pick.part13[ii]])
      log.titre[log.titre==-Inf]=0
      log.titreJ = as.numeric(fj_neut_corr[,pick.part[jj]]) - as.numeric(fj_neut_corr[,pick.part13[jj]])
      log.titreJ[log.titreJ==-Inf]=0
      if(ii<jj){
        c.test = cor.test(log.titreJ,log.titre)
        plot(log.titreJ,log.titre,main="",xlab=paste(pick.name[jj]),ylab=paste(pick.name[ii]),xlim=c(0,9),ylim=c(0,9),pch=19)
        lines(c(-1,9),c(-1,9),col="grey")
        text(x=0,y=8.5,adj=0,cex=0.7,labels=paste("r=",signif(c.test$estimate,2)," p=",signif(c.test$p.value,2),sep=""))
      }
      if(ii>jj){
        plot(log.titreJ,log.titreJ,pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
      }
      if(ii==jj){
        hist(log.titreJ,xlab=pick.name[ii],main=NULL,breaks=seq(-2.5,8.5,1)) 
      }
    }
  }
  
}


# Table 3 - seroconversion and seroreversion to ZIKV (Fiji) ----------
supp_tab1_create <- function(){
  fj_mcnemars <- fj_data %>%
    mutate(PP=as.numeric(ZIKV15+ZIKV17==2),
           PN=as.numeric(ZIKV15==1 & ZIKV17==0),
           NP=as.numeric(ZIKV17==1 & ZIKV15==0),
           NN=as.numeric(ZIKV15+ZIKV17==0)) %>%
    group_by(age) %>%
    summarise(PP=sum(PP), PN=sum(PN), NP=sum(NP), NN=sum(NN)) 
  fj_mcnemars_total <- summarise(fj_mcnemars, PP=sum(PP), PN=sum(PN), NP=sum(NP), NN=sum(NN)) %>%
    mutate(age=NA)
  fj_mcnemars <- rbind(fj_mcnemars, fj_mcnemars_total)
  fj_mcnemars_output <- matrix(NA, nrow = 6, ncol=6)
  fj_mcnemars_output[1,] <- c(fj_mcnemars[1,2],fj_mcnemars[1,3], rep(NA,4)) %>% as.numeric()
  fj_mcnemars_output[2,] <- c(fj_mcnemars[1,4],fj_mcnemars[1,5], rep(NA,4)) %>% as.numeric()
  fj_mcnemars_output[3,] <- c(NA,NA,fj_mcnemars[2,2],fj_mcnemars[2,3],NA,NA) %>% as.numeric()
  fj_mcnemars_output[4,] <- c(NA,NA,fj_mcnemars[2,4],fj_mcnemars[2,5],NA,NA) %>% as.numeric()
  fj_mcnemars_output[5,] <- c(rep(NA,4), fj_mcnemars[3,2],fj_mcnemars[3,3]) %>% as.numeric()
  fj_mcnemars_output[6,] <- c(rep(NA,4), fj_mcnemars[3,4],fj_mcnemars[3,5]) %>% as.numeric()
  fj_mcnemars_output
}


# Table 4 and 5 - Age and DENV exposure profile (French Poly.) & bootstrap estimates -------------
supp_tab2_tab3_create <- function(){
  
  fp_ind_data_14 <- fp_ind_data %>%
    filter(year=="14") 
  fp_ind_data_15 <- fp_ind_data %>% 
    filter(year=="15") 
  freq_14 <- signif(table(fp_ind_data_14$DENVpos)/length(fp_ind_data_14$DENVpos),2)
  freq_15 <- signif(table(fp_ind_data_15$DENVpos)/length(fp_ind_data_15$DENVpos),2)
  sample_bias <- t(rbind(paste0(table(fp_ind_data_14$DENVpos), " [",
                                freq_14, "]"),
                         paste0(table(fp_ind_data_15$DENVpos), " [",
                                freq_15, "]")))
  colnames(sample_bias) <- c("2014 [pc]","2015 [pc]")
  ages <- fp_ind_data_14$Age
  
  age_dist <- function(x){
    bp1=signif(c(median(x),quantile(x,0.25),quantile(x,0.75)),2)
    paste(bp1[1]," [",bp1[2],"-",bp1[3],"]",sep="")
  }
  
  sample_bias <- rbind(cbind(age_dist(fp_ind_data_14$Age),age_dist(fp_ind_data_15$Age)),
                        sample_bias)
  
  
  # Bootstrap sample from dengue exposure histories
  viruses_5 <- c("DENV1","DENV2","DENV3","DENV4","ZIKV")

  ln_14 <- length(fp_ind_data_14$Age)
  btstrp_positive_14 <- NULL
  ln.btstrp <- 1000
  set.seed(012346789)
  
  for(ii in 1:ln.btstrp){
    
    store_positive <- NULL
    
    for(jj in 1:5){
      
      # Pick viruses other than chosen one
      virus_pick <- viruses_5[viruses_5!=viruses_5[jj]]
      if(jj<5){virus_pick <- virus_pick[1:3]}
      
      sampl_freq <- table(rowSums(fp_ind_data_15[,virus_pick]))
      btstrp_dist <- sample(x = 0:(length(sampl_freq)-1), size = ln_14, prob = as.numeric(sampl_freq), replace=T) # Pick bootstrap samples
      
      data_14 <- rowSums(fp_ind_data_14[,virus_pick])
      
      # Obtain bootstrap sample
      btstrp_s <- t(sapply(btstrp_dist, function(x){id <- sample(which(data_14==x), size=1); as.numeric(fp_ind_data_14[id,])}))
      colnames(btstrp_s) <- names(fp_ind_data_14)
      
      store_positive <- c(store_positive, sum(btstrp_s[,viruses_5[jj]]))
      store_positive <- as.data.frame(store_positive)
    }
    
    names(store_positive) <- viruses_5
    
    btstrp_positive_14 <- rbind(btstrp_positive_14,store_positive)

  }

  save(btstrp_positive_14,file="store_tables/bootstrapFP14.RData",precheck=T)
  
  btstrp_sum <- function(v, sigF=3){
    btsrp_pos <- btstrp_positive_14[,v]/ln_14
    bp1=signif(c(mean(btsrp_pos),quantile(btsrp_pos,0.025),quantile(btsrp_pos,0.975)),sigF)
    paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
  }
  
  param1=cbind(
    btstrp_sum('DENV1',2),
    btstrp_sum('DENV2',2),
    btstrp_sum('DENV3',2),
    btstrp_sum('DENV4',2),
    btstrp_sum('ZIKV',2)
  )
  rownames(param1)="2014 bootstrap estimates"
  colnames(param1)=c('DENV1',
                     'DENV2',
                     'DENV3',
                     'DENV4',
                     'ZIKV')
  
  ## summarise 2014 and 2015 real data
  summarise_FP <- function(VIRUS){
    fp_adults_14 <- fp_data %>%
      filter(year=="14",pop=="general") %>%
      mutate(sample_size=get(paste0(VIRUS,"P"))+get(paste0(VIRUS,"N"))) %>%
      summarise(n=sum(sample_size), pos=sum(get(paste0(VIRUS,"P")),na.rm=T))
    CI.calc <- binom.test(x=fp_adults_14$pos, n=fp_adults_14$n)          
    fp_adults_14 <- cbind(fp_adults_14, m=CI.calc$estimate, lci=CI.calc$conf.int[1], uci=CI.calc$conf.int[2], year=2014.3)
    #2015
    fp_adults_15 <- fp_data %>%
      filter(year=="15",pop=="general") %>%
      mutate(sample_size=get(paste0(VIRUS,"P"))+get(paste0(VIRUS,"N"))) %>%
      summarise(n=sum(sample_size), pos=sum(get(paste0(VIRUS,"P")),na.rm=T))
    CI.calc <- binom.test(x=fp_adults_15$pos, n=fp_adults_15$n)          
    fp_adults_15 <- cbind(fp_adults_15, m=CI.calc$estimate, lci=CI.calc$conf.int[1], uci=CI.calc$conf.int[2], year=2015.9)
    fp_adults <- rbind(fp_adults_14,fp_adults_15)
    rownames(fp_adults) <- NULL
    fp_adults <- signif(fp_adults,2)
    fp_adults_supp_table <- t(rbind(paste0(fp_adults[1,3]," (",fp_adults[1,4],"-",fp_adults[1,5],")"),
                                    paste0(fp_adults[2,3]," (",fp_adults[2,4],"-",fp_adults[2,5],")")))
    rownames(fp_adults_supp_table) <- VIRUS
    colnames(fp_adults_supp_table) <- c("2014","2015")
    fp_adults_supp_table
  }
  
    param_fp <- rbind.data.frame(summarise_FP("DENV1"),
                               summarise_FP("DENV2"),
                               summarise_FP("DENV3"),
                               summarise_FP("DENV4"),
                               summarise_FP("ZIKV")
  )
  
  supp_table_bootstrap <- cbind.data.frame("2014"=param_fp$`2014`,
                                           t(param1),
                                           "2015"=param_fp$`2015`)
  
    
  s <- btstrp_positive_14 %>% as.data.frame() %>% 
    mutate_at(.funs=~./ln_14, .vars=1:5) %>%
    summarise_at(.funs=~mean(.), .vars=1:5) %>%
    mutate_at(.funs=~round(.*ln_14), .vars=1:5) %>%
    as.numeric()
  s2 <- round(c(0.8,0.18,0.55,0.42,0.22)*672) ## 672 adults in the general population sampled in 2015
  retest <- sapply(1:5,function(x){matrix(c(s[x],49-s[x],s2[x],672-s2[x]), byrow=T, nrow=2) %>% chisq.test()}) 
  p.vals <- signif(as.numeric(c(retest['p.value',1],retest['p.value',2],retest['p.value',3],retest['p.value',4],retest['p.value',5])),2) 
  supp_table_bootstrap <- cbind(supp_table_bootstrap, p.vals)
  supp_table_bootstrap
  
  return(list(supp_tab2=sample_bias,supp_tab3=supp_table_bootstrap))
}

# Table 6 - Change in neutralization titre to DENV3 and ZIKV  --------
supp_tab4_create <- function(virus1,virus2){
  # virus1 <- "Zs"; virus2 <- "D3s"
  # virus1 <- "Zs"; virus2 <- "D1s"
  fj_neut_virus1 <- fj_neut_data %>%
    mutate_at(.funs=~pmax(.,0), .vars=3:17) %>%
    filter(!is.na(get(paste0(virus1,"13"))) & get(paste0(virus1,"13"))<2 & (get(paste0(virus1,"15"))-get(paste0(virus1,"13")))>=2) %>%
    mutate(rise_1=get(paste0(virus1,"15"))-get(paste0(virus1,"13")),
           drop_1=get(paste0(virus1,"17"))-get(paste0(virus1,"15")))
  fj_neut_virus2 <- fj_neut_data %>%
    mutate_at(.funs=~pmax(.,0), .vars=3:17) %>%
    filter(!is.na(get(paste0(virus2,"13"))) & get(paste0(virus2,"13"))<2 & (get(paste0(virus2,"15"))-get(paste0(virus2,"13")))>=2) %>%
    mutate(rise_2=get(paste0(virus2,"15"))-get(paste0(virus2,"13")),
           drop_2=get(paste0(virus2,"17"))-get(paste0(virus2,"15")))
  
  ## Tests
  test.vals3a <- fj_neut_virus2$rise_2
  test.vals3b <- fj_neut_virus2$drop_2
  test.vals5a <- fj_neut_virus1$rise_1
  test.vals5b <- fj_neut_virus1$drop_1
  t_test2 = c(t.test(test.vals3a)$p.value,t.test(test.vals5a)$p.value,
              t.test(test.vals3a)$estimate,t.test(test.vals5a)$estimate)
  t_test2 = matrix(t_test2,ncol=2)
  rownames(t_test2) = c("d3_rise","z_rise") 
  colnames(t_test2) = c("p","mean")
  t_test2
  
  # T_test for waning of DENV-3 and ZIKV
  t_test1 = c(t.test(test.vals3b)$p.value,t.test(test.vals5b)$p.value)
  t_test1 = matrix(t_test1,ncol=1)
  rownames(t_test1) = c("d3_wane","z_wane") #"d2_d3_rise","d2_z_rise","d3_z_rise",
  
  d15test <- t.test(test.vals3a)
  z15test <- t.test(test.vals5a)
  (d17test <- t.test(test.vals3b))
  (z17test <- t.test(test.vals5b))
  
  # averages 
  get_avg <- function(x){
    paste0(signif(mean(x),2), " [",
           signif(quantile(x,0.025),2),",",
           signif(quantile(x,0.975),2),"]")
  }
  get_testCIs <- function(test){
    paste0(signif(test$estimate,2), " [",
           signif(test$conf.int[1],2),",",
           signif(test$conf.int[2],2),"]"
    )
  }
  
  
  ## table for manuscript
  supptable1 <- data.frame("Virus"=c(virus1,virus2),
                           "change 2013-2015"=rep(NA,2),
                           "p1"=rep(NA,2),
                           "change 2015-2017"=rep(NA,2),
                           "p2"=rep(NA,2))
  supptable1$change.2013.2015 <- 
    c(get_testCIs(z15test),
      get_testCIs(d15test))
  supptable1$change.2015.2017 <- 
    c(get_testCIs(z17test),
      get_testCIs(d17test))
  supptable1$p1 <- c(signif(t.test(test.vals3a,test.vals5a)$p.value,3),
                     "")
  supptable1$p2 <- c(signif(t.test(test.vals3b,test.vals5b)$p.value,3),
                     "")
  supptable1
  
}



