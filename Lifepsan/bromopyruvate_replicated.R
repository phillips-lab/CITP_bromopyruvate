# Set Up ------------------------------------------------------------------

# load packages
library(readxl); library(coxme); library(multcomp); 
library(lme4); library(lubridate)

# set active directory
activedir <- file.choose()
setwd(dirname(activedir))

# load data
comps <- read_excel("bromopyruvate_replicated.xlsx")
date <- format(today(), "%m%d%y")

comps$Lab <- as.factor(comps$Lab)
comps$StartDate <- as.factor(comps$StartDate)
comps$Tech <- as.factor(comps$Tech)
comps$plate_id <- as.factor(comps$plate_id)
comps$Species <- as.factor(comps$Species)
comps$Strain <- as.factor(comps$Strain)
comps$Concentration <- as.factor(comps$Concentration)
comps$Rep <- as.factor(comps$Rep)
comps$Comp_Conc <- droplevels(with(comps, interaction(Compound,Concentration, sep="x")))

levels(comps$Strain)
levels(comps$Comp_Conc)

strain <- c("AF16", "ED3092", "HK104","JU775", "MY16", "N2_PD1073")
header <- c("af16", "ed3092", "hk104", "ju775", "my16", "n2pd1073")

comps_cc_pairs <- c("bromopyruvatex500 - CTRL_H2Ox0 = 0")

ivector <- c(1:6)

# subset data
for (i in ivector) {
  assign(paste(header[i],sep=""),droplevels(comps[ which(comps$Strain==strain[i]), ]))
}


# Cox ---------------------------------------------------------------------

for (i in ivector) {
  assign(paste(header[i],"cox",sep="_"),coxme(Surv(DeathAge,Dead) ~ Comp_Conc + (1|Lab/StartDate/Tech/plate_id), data=get(paste(header[i],sep="_"))))
  assign(paste(header[i],"cox","pairs",sep="_"),summary(glht(get(paste(header[i],"cox",sep="_")), linfct = mcp(Comp_Conc = get(paste("comps","cc","pairs",sep="_"))))))
}


# GLM ---------------------------------------------------------------------

for (i in ivector) {
  assign(paste(header[i],"lm",sep="_"),lmer(DeathAge ~ Comp_Conc + (1|Lab/StartDate/Tech/plate_id), data=get(paste(header[i],sep="_"))))
  assign(paste(header[i],"lm","prof",sep="_"),profile(get(paste(header[i],"lm",sep="_"))))
  assign(paste(header[i],"lm","prof","CI",sep="_"),confint(get(paste(header[i],"lm","prof",sep="_"))))
  assign(paste(header[i],"lm","pairs",sep="_"),summary(glht(get(paste(header[i],"lm",sep="_")), linfct = mcp(Comp_Conc = get(paste("comps","cc","pairs",sep="_"))))))
}


# text output -------------------------------------------------------------

ivector <- 1:6

for (i in ivector) {
  sink(paste(header[i],date,".txt",sep="_"), append = FALSE)
  writeLines(strain[i])  
  writeLines("\nCox\n")
  print(get(paste(header[i],"cox",sep="_")))
  print(get(paste(header[i],"cox","pairs",sep="_")))
  writeLines("\nLinear Model\n")
  print(summary(get(paste(header[i],"lm",sep="_"))))
  print(get(paste(header[i],"lm","prof","CI",sep="_")))
  print(get(paste(header[i],"lm","pairs",sep="_")))
  sink()
} 


# Quantiles and summary --------------------------------------------------------------

summary_list <- list()

for (i in ivector) {
  survfit_obj <- survfit(Surv(DeathAge, Dead) ~ Comp_Conc, data = get(paste(header[i], sep = "_")))
  assign(paste(header[i], "survfit", sep = "_"), survfit_obj)
  
  # calculate quantiles
  quantiles <- summary(survfit_obj)$table
  
  print(paste("Strain:", header[i]))
  print(quantiles)
  
  # calculate the 90% quantile
  quantile_90_list <- sapply(levels(get(paste(header[i], sep = "_"))$Comp_Conc), function(comp_conc) {
    survfit_obj_comp <- survfit(Surv(DeathAge, Dead) ~ 1, data = get(paste(header[i], sep = "_"))[get(paste(header[i], sep = "_"))$Comp_Conc == comp_conc, ])
    quantile(survfit_obj_comp, probs = 0.9)$quantile[1]
  })
  
  
  summary_df <- data.frame(
    Strain = header[i],
    Comp_Conc = rownames(quantiles),
    n_dead = quantiles[,"events"],
    n_censored = quantiles[,"records"] - quantiles[,"events"],
    n_total = quantiles[,"records"],
    median = quantiles[,"median"],
    LCL_95 = quantiles[,"0.95LCL"],
    UCL_95 = quantiles[,"0.95UCL"],
    rmean = quantiles[,"rmean"],
    se_rmean = quantiles[,"se(rmean)"],
    quantile_90 = quantile_90_list
  )
  
  summary_list[[i]] <- summary_df
}

summary_results <- do.call(rbind, summary_list)

# export to a CSV file
write.csv(summary_results, "summary_results.csv", row.names = FALSE)


# Linear Model ------------------------------------------------------------

# For overall sources of variance

lm_all <-lmer(DeathAge ~ Compound + (1|Lab/StartDate/Tech/plate_id) + (1|Species/Strain) + (1|Lab:Strain) + (1|Lab:Species) + (1|Compound:Species) + (1|Compound:Strain) + (1|Compound:Lab), data=comps)
summary(lm_all)

sink(paste("bromop_variance_",date,".txt",sep=""), append = FALSE)
writeLines("\nSources of Variance - Overall\n")
summary(lm_all)
#lm_all_prof_CI
sink()


# Functions for Output ----------------------------------------------------

# multiple comparisons, data table
pairsoutput <- function(modelmultcomp,model,source){
  est <- as.data.frame(modelmultcomp[["test"]][["coefficients"]])
  stderr <- as.data.frame(modelmultcomp[["test"]][["sigma"]])
  zval <- as.data.frame(modelmultcomp[["test"]][["tstat"]])
  pval <- as.data.frame(modelmultcomp[["test"]][["pvalues"]])
  output <- cbind(est,stderr,zval,pval,model,source)
  output$row <- row.names(output)
  output <- setNames(output,c("est","stderr","zval","pval","model","strain","conc_ctrl"))
  return(output)
}


# random effects variance 
lmoutput <- function(origlmmodel,model,source){
  output <- cbind(as.data.frame(summary(origlmmodel)[["varcor"]]),model,source)
  output <- setNames(output,c("group","intercept","na","variance","stddev","model","strain"))
  return(output)
}

# random effects standard deviation confidence intervals
lmcioutput <- function(origlmci,model,source){
  output <- cbind(as.data.frame(origlmci),model,source)
  output <- setNames(output,c("2.5%","97.5%","model","strain"))
  output$row <- row.names(output)
  return(output)
}


# random effects variance
coxoutput <- function(origcoxmodel,model,source){
  output <- cbind(as.data.frame(origcoxmodel[["vcoef"]]),model,source)
  output <- setNames(output,c("Lab. Date.Tech.Plate","Lab.Date.Tech","Lab.Date","Lab", "model","strain")) 
  return(output)
}


# Output by Function ------------------------------------------------------
header <- c("af16", "ed3092", "hk104", "ju775", "my16", "n2pd1073")

## Multiple Comparisons
af16pairs <- rbind(pairsoutput(af16_cox_pairs,"cox","af16"),pairsoutput(af16_lm_pairs,"lm","af16"))
ed3092pairs <- rbind(pairsoutput(ed3092_cox_pairs,"cox","ed3092"),pairsoutput(ed3092_lm_pairs,"lm","ed3092"))
hk104pairs <- rbind(pairsoutput(hk104_cox_pairs,"cox","hk104"),pairsoutput(hk104_lm_pairs,"lm","hk104"))
ju775pairs <- rbind(pairsoutput(ju775_cox_pairs,"cox","ju775"),pairsoutput(ju775_lm_pairs,"lm","ju775"))
my16pairs <- rbind(pairsoutput(my16_cox_pairs,"cox","my16"),pairsoutput(my16_lm_pairs,"lm","my16"))
n2pd1073pairs <- rbind(pairsoutput(n2pd1073_cox_pairs,"cox","n2pd1073"),pairsoutput(n2pd1073_lm_pairs,"lm","n2pd1073"))


allpairs <- rbind(af16pairs,ed3092pairs,hk104pairs,ju775pairs,my16pairs,n2pd1073pairs)
allpairs
write.csv(allpairs,"MultipleComparisonsOutput.csv")

## Variance of the random effects, cox model only
coxset1 <- rbind(coxoutput(af16_cox,"cox","af16"),coxoutput(ed3092_cox,"cox","ed3092"),coxoutput(hk104_cox,"cox","hk104"))
coxset2 <- rbind(coxoutput(my16_cox,"cox","my16"),coxoutput(n2pd1073_cox,"cox","n2pd1073"),coxoutput(ju775_cox,"cox","ju775"))
allcoxvar <- t(rbind(coxset1,coxset2))
allcoxvar
write.csv(allcoxvar,"VarianceRandomEffectsCoxPH.csv")

## Variance of the random effects, linear model only (SINGLE FILE)
lmset1 <- rbind(lmoutput(af16_lm,"lm","af16"),lmoutput(ed3092_lm,"lm","ed3092"),lmoutput(hk104_lm,"lm","hk104"))
lmset2 <- rbind(lmoutput(my16_lm,"lm","my16"),lmoutput(n2pd1073_lm,"lm","n2pd1073"),lmoutput(ju775_lm,"lm","ju775"))
alllmvar <- rbind(lmset1,lmset2) 
alllmvar
write.csv(alllmvar,"VarianceRandomEffectsLM.csv")

## Confidence Intervals of the variance of the random effects, linear model only
# only care about rows sig01 through sigma
lmciset1 <- rbind(lmcioutput(af16_lm_prof_CI,"lm_CI","af16"),lmcioutput(ed3092_lm_prof_CI,"lm_CI","ed3092"),lmcioutput(hk104_lm_prof_CI,"lm_CI","hk104"))
lmciset2 <- rbind(lmcioutput(my16_lm_prof_CI,"lm_CI","my16"),lmcioutput(n2pd1073_lm_prof_CI,"lm_CI","n2pd1073"),lmcioutput(ju775_lm_prof_CI,"lm_CI","ju775"))
alllmvarci <- rbind(lmciset1,lmciset2) 
alllmvarci
write.csv(alllmvarci,"VarianceRandomEffectsLMCI.csv")
