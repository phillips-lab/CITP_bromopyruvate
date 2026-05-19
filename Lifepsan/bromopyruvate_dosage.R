# Set up ------------------------------------------------------------------

library(readxl); library(coxme); library(multcomp); 
library(lme4); library(lubridate); 

activedir <- file.choose()
setwd(dirname(activedir))

comps <- read_excel("bromopyruvate_replicated.xlsx")
date <- format(today(), "%m%d%y")

origcomps <- comps

comps$Lab <- as.factor(comps$Lab)
comps$StartDate <- as.factor(comps$StartDate)
comps$Tech <- as.factor(comps$Tech)
comps$Compound <- as.factor(comps$Compound)
comps$Concentration <- as.factor(comps$Concentration)
comps$Rep <- as.factor(comps$Rep)
comps$death_id <- as.factor(comps$death_id)
comps$observation_id <- as.factor(comps$observation_id)
comps$plate_id <- as.factor(comps$plate_id)
comps$experiment_id <- as.factor(comps$experiment_id)

comps$Comp_Conc <- droplevels(with(comps, interaction(Compound,Concentration, sep="x")))

strain <- c("ED3092","JU1373","N2_PD1073")
header <- c("ed3092","ju1373","n2_pd1073")

levels(comps$Comp_Conc)

# set up multiple comparisons
comps_cc <- levels(comps$Comp_Conc)
assign(paste("comps","cc","pairs",sep="_"), c())
assign(paste("comps","cc","len",sep="_"), length(get(paste("comps","cc",sep="_"))))
for (k in 1:(comps_cc_len-1)) {
  comps_cc_pairs[k] <- c(paste(comps_cc[comps_cc_len+1-k],"-",comps_cc[1], " = 0"))
}



ivector <- c(1:3)

# subset data
for (i in ivector) {
  assign(paste(header[i],sep=""),droplevels(comps[ which(comps$Strain==strain[i]), ]))
}

# Cox ---------------------------------------------------------------------

for (i in ivector) {
  assign(paste(header[i],"cox",sep="_"),coxme(Surv(DeathAge,Dead) ~ Comp_Conc + (1|Tech/plate_id), data=get(paste(header[i],sep="_"))))
  assign(paste(header[i],"cox","pairs",sep="_"),summary(glht(get(paste(header[i],"cox",sep="_")), linfct = mcp(Comp_Conc = get(paste("comps","cc","pairs",sep="_"))))))
}


# GLM ---------------------------------------------------------------------

for (i in ivector) {
  assign(paste(header[i],"lm",sep="_"),lmer(DeathAge ~ Comp_Conc + (1|Tech/plate_id), data=get(paste(header[i],sep="_"))))
  assign(paste(header[i],"lm","prof",sep="_"),profile(get(paste(header[i],"lm",sep="_"))))
  assign(paste(header[i],"lm","prof","CI",sep="_"),confint(get(paste(header[i],"lm","prof",sep="_"))))
  assign(paste(header[i],"lm","pairs",sep="_"),summary(glht(get(paste(header[i],"lm",sep="_")), linfct = mcp(Comp_Conc = get(paste("comps","cc","pairs",sep="_"))))))
}

# text output -------------------------------------------------------------

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
  # Create the survfit object and assign it to a variable
  survfit_obj <- survfit(Surv(DeathAge, Dead) ~ Comp_Conc, data = get(paste(header[i], sep = "_")))
  assign(paste(header[i], "survfit", sep = "_"), survfit_obj)
  
  # Calculate quantiles on the survfit object
  quantiles <- summary(survfit_obj)$table
  
  print(paste("Strain:", header[i]))
  print(quantiles)
  
  # Calculate the 90% quantile
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

# Export to a CSV file
write.csv(summary_results, "summary_results.csv", row.names = FALSE)
