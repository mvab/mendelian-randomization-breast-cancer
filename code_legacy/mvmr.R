library(tidyr)
library(tibble)
library(dplyr)
#install.packages("devtools")
library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
library(MendelianRandomization)
#use Wes's MVMR package
library(devtools)
#install_github("WSpiller/MVMR")
library(MVMR)


# MVMR with exposure1, exposure2 and outcome
# Using TwoSampleMR package 

# Exposures
id_exposures <- c("ukb-b-3768", "ukb-b-12405")

#Read in / extract summary data for top SNPs in exposure 1 
exposure1 <- extract_instruments(c("ukb-b-3768"))
#Read in / extract summary data for top SNPs in exposure 2 
exposure2 <- extract_instruments(c("ukb-b-12405"))
exposures <- rbind(exposure1, exposure2)

# First obtain the instruments for each exposure
exposure_dat <- mv_extract_exposures(id_exposures)

#Next, also extract those SNPs from the outcome.
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 'ieu-a-1126')

#Once the data has been obtained, harmonise so that all are on the same reference allele.
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

#Finally, perform the multivariable MR analysis
res <- mv_multiple(mvdat)

mv_res<- res$result %>%
  split_outcome() %>% 
  split_exposure() %>% 
  separate(outcome, "outcome", sep="[(]") %>% 
  generate_odds_ratios()
mv_res


# To be able to use other packages for MVMR, need to convert data to a different format
# below is the shortcut to do it by using `exposure_dat` from 2SMR 
# but all manual steps (the long way) to do it are also below

# * beginning of shortcut * 
## NB: shortcut from `exposure_dat` in 2SMR to `exposures_joined` (here it is `exposures_joined_auto`)

# restructure exposure data to wide format
exposures_joined_auto <- exposure_dat %>%
  select(SNP, beta.exposure, se.exposure, exposure) %>% 
  # convert to wider format
  pivot_wider(names_from = exposure, values_from = c(beta.exposure, se.exposure)) %>% 
  # rearrange cols to be SNP, betaX1, seX1, betaX2, seX2
  select(1,2,4,3,5)
colnames(exposures_joined_auto)<-c('SNP', 'betaX1', 'seX1', 'betaX2', 'seX2')

# extract SNPs for both exposures from outcome dataset 
outcome_dat <- extract_outcome_data(snps = exposures_joined_auto$SNP,
                                    outcomes ="ieu-a-1126")
# harmonize datasets 
outcome_harmonised <- harmonise_data(exposures, outcome_dat)


# Create variables for MV
# remove factors structure in SNPs, and sort by SNP (YGs and XGs must have SNPs in the same order)
YG <- outcome_harmonised %>% 
  select("SNP", "beta.outcome", "se.outcome") %>%
  distinct() %>%
  mutate(SNP = as.character(SNP)) %>%
  arrange((SNP)) 
XGs <- exposures_joined_auto %>% 
  filter(SNP %in% outcome_harmonised$SNP) %>%
  mutate(SNP = as.character(SNP)) %>%
  arrange((SNP)) 

# some checks
stopifnot(dim(XGs)[1]==dim(YG)[1])
unique(YG$SNP %in% XGs$SNP)
unique(XGs$SNP %in% YG$SNP)
all.equal(YG$SNP, XGs$SNP)

# *end of shortcut *



# MVMR
mvmr_out <- format_mvmr(BXGs = XGs[,c("betaX1", "betaX2")],
                        BYG = YG$beta.outcome,
                        seBXGs = XGs[,c("seX1", "seX2")],
                        seBYG = YG$se.outcome,
                        RSID = XGs$SNP)
mvmr_res <- mvmr(mvmr_out, 0, 1)

res_df <- mvmr_res$coef %>% as.data.frame() %>% 
  rownames_to_column("exposure") %>% 
  rename(b=Estimate) %>% 
  rename(se="Std. Error") %>% 
  rename(pval="Pr(>|t|)") %>% 
  generate_odds_ratios()

res_df$exposure<-unique(exposure_dat$exposure)
res_df


# MVIVW
mr_mvivw_res <- mr_mvivw(mr_mvinput(bx = as.matrix(XGs[,c("betaX1", "betaX2")]), 
                                    bxse = as.matrix(XGs[,c("seX1", "seX2")]),
                                    by = YG$beta.outcome, 
                                    byse =YG$se.outcome, 
                                    exposure = unique(exposure_dat$exposure)))

get_odds_ratios_mvivw<-function(mr_result){
  # create empty df
  df<-setNames(data.frame(matrix(ncol = 9, nrow = 2)),
               c("exposure", 'b', 'se',  'pval' , 'lo_ci', 'up_ci' , 'or', 'or_lci95', 'or_uci95'))
  # extract results from MR object
  df$exposure <- mr_result$Exposure
  df$b <- mr_result$Estimate
  df$se <- mr_result$StdError
  df$pval <- mr_result$Pvalue
  # calculated OR and CIs
  df$lo_ci <- df$b  - 1.96 * df$se
  df$up_ci <- df$b  + 1.96 * df$se
  df$or <- round(exp(df$b),4)
  df$or_lci95<-round(exp(df$lo_ci),4)
  df$or_uci95<-round(exp(df$up_ci),4)
  return(df)
}
print(get_odds_ratios_mvivw(mr_mvivw_res))

#MVMR-Egger
mr_mvegger_res <- mr_mvegger(mr_mvinput(bx = as.matrix(XGs[,c("betaX1", "betaX2")]), 
                                        bxse = as.matrix(XGs[,c("seX1", "seX2")]),
                                        by = YG$beta.outcome, 
                                        byse =YG$se.outcome,
                                        exposure = unique(exposure_dat$exposure)),
                             orientate = 1)
get_odds_ratios_mvegger<-function(mr_result){
  # create empty df
  df<-setNames(data.frame(matrix(ncol = 9, nrow = 2)),
               c("exposure", 'b', 'se',  'pval' , 'lo_ci', 'up_ci' , 'or', 'or_lci95', 'or_uci95'))
  # extract results from MR object
  df$exposure <- mr_result$Exposure
  df$b <- mr_result$Estimate
  df$se <- mr_result$StdError.Est
  df$pval <- mr_result$Pvalue.Est
  # calculated OR and CIs
  df$lo_ci <- df$b  - 1.96 * df$se
  df$up_ci <- df$b  + 1.96 * df$se
  df$or <- round(exp(df$b),4)
  df$or_lci95<-round(exp(df$lo_ci),4)
  df$or_uci95<-round(exp(df$up_ci),4)
  return(df)
}
print(get_odds_ratios_mvegger(mr_mvegger_res))








## This section shows how to create exposure data without using 2SMR package
# might be useful for undestanding what happens internally 




# Exposures

#Read in / extract summary data for top SNPs in exposure 1 
exposure1 <- extract_instruments(c("ukb-b-3768"))
#Read in / extract summary data for top SNPs in exposure 2 
exposure2 <- extract_instruments(c("ukb-b-12405"))
exposures <- rbind(exposure1, exposure2)


# manually extracting all data to use with MVMR and MendelianRandomization packages. 
extract_exp_out_pair_SNPs <- function(exposure1.data, exposure2.id){
  ### this function gets GWAS1 from GWAS 2

  # extract exposure 1 SNPs from exposure 2 
  outcome <- extract_outcome_data(snps = exposure1.data$SNP,
                                  outcomes = exposure2.id) 
  # identify SNPs that occur >1 and drop them (i.e exclude multiallelic SNPs)
  n_occur <- data.frame(table(outcome$SNP))
  occur_multiple <- n_occur[n_occur$Freq > 1,]$Var1
  outcome = outcome[!outcome$SNP %in% occur_multiple,]
  
  # harmonize datasets 
  dat <- harmonise_data(exposure1.data, outcome)
  dat <- dat[!(dat$remove==TRUE),]
  return(dat)
}

# get exposure1 SNPs in exposure2  
dat1 <- extract_exp_out_pair_SNPs(exposure1.data = exposure1, exposure2.id = "ukb-b-12405")
# get exposure2 SNPs in exposure1 
dat2 <- extract_exp_out_pair_SNPs(exposure1.data = exposure2, exposure2.id =  "ukb-b-3768")


# get exposure1 SNPs and exposure2 SNPs in GWAS1
# exp1 SNPs in GWAS1
exp1_gwas1 <-dat1 %>%
  select(SNP, beta.exposure, se.exposure) %>%
  rename(betaX1 = beta.exposure) %>% 
  rename(seX1 = se.exposure)
# exp2 SNPs in GWAS1
exp2_gwas1<-dat2 %>%
  select(SNP, beta.outcome, se.outcome) %>%
  rename(betaX1 = beta.outcome) %>%
  rename(seX1 = se.outcome)
dat_gwas1 <- rbind(exp1_gwas1, exp2_gwas1) %>% distinct()

# get exposure1 SNPs and exposure2 SNPs in GWAS2
# exp2 SNPs in GWAS2
exp2_gwas2<-dat2 %>%
  select(SNP, beta.exposure, se.exposure) %>%
  rename(betaX2 = beta.exposure) %>%
  rename(seX2 = se.exposure)
# exp1 SNPs in GWAS2
exp1_gwas2<-dat1 %>%
  select(SNP, beta.outcome, se.outcome) %>%
  rename(betaX2 = beta.outcome) %>%
  rename(seX2 = se.outcome)
dat_gwas2<-rbind(exp2_gwas2, exp1_gwas2) %>% distinct()

# join data
exposures_joined <- full_join(dat_gwas1, dat_gwas2, by="SNP") %>% distinct()
dim(exposures_joined)


# Outcomes

# extract SNPs for both exposures from outcome dataset 
outcome_dat <- extract_outcome_data(snps = exposures_joined$SNP,
                                    outcomes ="ieu-a-1126")
# harmonize datasets 
outcome_harmonised <- harmonise_data(exposures, outcome_dat)



# Create variables for MV
YG <- outcome_harmonised %>% 
  select("SNP", "beta.outcome", "se.outcome") %>% 
  distinct() %>%
  mutate(SNP = as.character(SNP)) %>% 
  arrange((SNP)) 
XGs <- exposures_joined %>% 
  filter(SNP %in% outcome_harmonised$SNP) %>% 
  mutate(SNP = as.character(SNP)) %>% 
  arrange((SNP))


stopifnot(dim(XGs)[1]==dim(YG)[1])
unique(YG$SNP %in% XGs$SNP)
unique(XGs$SNP %in% YG$SNP)

# MVMR
mvmr_out <- format_mvmr(BXGs = XGs[,c("betaX1", "betaX2")],
                        BYG = YG$beta.outcome,
                        seBXGs = XGs[,c("seX1", "seX2")],
                        seBYG = YG$se.outcome,
                        RSID = XGs$SNP)
mvmr_res <- mvmr(mvmr_out, 0, 1)

res_df <- mvmr_res$coef %>% as.data.frame() %>% 
  rownames_to_column("exposure") %>% 
  rename(b=Estimate) %>% 
  rename(se="Std. Error") %>% 
  rename(pval="Pr(>|t|)") %>% 
  generate_odds_ratios()
res_df


#MVMR-IVW
mr_mvivw_res <- mr_mvivw(mr_mvinput(bx = as.matrix(XGs[,c("betaX1", "betaX2")]), 
                                bxse = as.matrix(XGs[,c("seX1", "seX2")]),
                                by = YG$beta.outcome, 
                                byse =YG$se.outcome, 
                                exposure = c("ukb-b-3768", "ukb-b-12405")))

get_odds_ratios_mvivw<-function(mr_result){
  # create empty df
  df<-setNames(data.frame(matrix(ncol = 9, nrow = 2)),
               c("exposure", 'b', 'se',  'pval' , 'lo_ci', 'up_ci' , 'or', 'or_lci95', 'or_uci95'))
  # extract results from MR object
  df$exposure <- mr_result$Exposure
  df$b <- mr_result$Estimate
  df$se <- mr_result$StdError
  df$pval <- mr_result$Pvalue
  # calculated OR and CIs
  df$lo_ci <- df$b  - 1.96 * df$se
  df$up_ci <- df$b  + 1.96 * df$se
  df$or <- round(exp(df$b),4)
  df$or_lci95<-round(exp(df$lo_ci),4)
  df$or_uci95<-round(exp(df$up_ci),4)
  return(df)
}
print(get_odds_ratios_mvivw(mr_mvivw_res))

#MVMR-Egger
mr_mvegger_res <- mr_mvegger(mr_mvinput(bx = as.matrix(XGs[,c("betaX1", "betaX2")]), 
                                    bxse = as.matrix(XGs[,c("seX1", "seX2")]),
                                    by = YG$beta.outcome, 
                                    byse =YG$se.outcome,
                                    exposure = c("ukb-b-3768", "ukb-b-12405")),
                                    orientate = 1)
get_odds_ratios_mvegger<-function(mr_result){
  # create empty df
  df<-setNames(data.frame(matrix(ncol = 9, nrow = 2)),
               c("exposure", 'b', 'se',  'pval' , 'lo_ci', 'up_ci' , 'or', 'or_lci95', 'or_uci95'))
  # extract results from MR object
  df$exposure <- mr_result$Exposure
  df$b <- mr_result$Estimate
  df$se <- mr_result$StdError.Est
  df$pval <- mr_result$Pvalue.Est
  # calculated OR and CIs
  df$lo_ci <- df$b  - 1.96 * df$se
  df$up_ci <- df$b  + 1.96 * df$se
  df$or <- round(exp(df$b),4)
  df$or_lci95<-round(exp(df$lo_ci),4)
  df$or_uci95<-round(exp(df$up_ci),4)
  return(df)
}
print(get_odds_ratios_mvegger(mr_mvegger_res))



