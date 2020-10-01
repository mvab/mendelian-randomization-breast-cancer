library(data.table)
library(dplyr)
library(TwoSampleMR)

#data_path<-"/projects/MRC-IEU/research/projects/ieu2/p1/055/working/data/"
#data_path <- "/Volumes/055/working/data/"
#file_gwas <- paste0(data_path, "GWAS_results_new/BMI_10_adj_month_f_imputed.txt.gz")
data_path<-"/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/Data/"
file_gwas <- paste0(data_path, "GWAS_results/BMI_10_adj_month_f_imputed.txt.gz")

gwas_data_raw <-vroom::vroom(file_gwas, col_select  = c("SNP","CHR",'BP',"BETA","SE","ALLELE1","ALLELE0","A1FREQ","P_BOLT_LMM_INF"))
dim(gwas_data_raw) #12,320,838

# filter to 10e-8
gwas_data <-  gwas_data_raw %>%  dplyr::filter(P_BOLT_LMM_INF < 5e-8) 
dim(gwas_data) 

# format data to 'outcome' and clump with defaults
gwas_data_tidy_out<- format_data(gwas_data, 
                                 type="outcome",
                                 snp_col = "SNP",
                                 chr_col = "CHR",
                                 pos_col = "BP",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "ALLELE1",
                                 other_allele_col = "ALLELE0",
                                 eaf_col = "A1FREQ",
                                 pval_col = "P_BOLT_LMM_INF")

instruments_OUT <- gwas_data_tidy_out %>%
                    clump_data() 
dim(instruments_OUT) # 150

instruments_OUT <- convert_outcome_to_exposure(instruments_OUT)


# format data to 'exposure' and clump with defaults
gwas_data_tidy_exp<- format_data(gwas_data, 
                                 type="exposure",
                                 snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "ALLELE1",
                                 other_allele_col = "ALLELE0",
                                 eaf_col = "A1FREQ",
                                 pval_col = "P_BOLT_LMM_INF")


instruments_EXP <- gwas_data_tidy_exp %>%  
                    clump_data() 
dim(instruments_EXP) # 132

# clump 'exposure former with old r2=0.01
instruments_EXP2 <- gwas_data_tidy_exp %>%
                    clump_data(., clump_r2 = 0.01) 
dim(instruments_EXP2) #176


# read in data from Supplemetary table S10 (early BMI adj)
#real_early_file <-paste0(data_path,"GWAS_tophits/earlyBMI_adj_UKB_tophits.csv") 
real_early_file <-paste0(data_path,"GWAS_tophits/all_good/adultBMI_UKB_tophits.csv") 

early_bmi <- read_exposure_data(
  filename = real_early_file,
  sep = ",",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "Standard Error",
  effect_allele_col = "Effect allele",
  other_allele_col = "Other allele",
  eaf_col = "Effect allele frequency",
  pval_col = "P"
)
dim(early_bmi) # 138


## compare SNP overlaps
length(intersect(instruments_EXP$SNP, instruments_OUT$SNP)) # 22 

length(intersect(early_bmi$SNP, instruments_OUT$SNP)) # 11
length(intersect(early_bmi$SNP, instruments_EXP$SNP)) # 108
length(intersect(early_bmi$SNP, instruments_EXP2$SNP)) # 121 <-- most similar


# check if early_bmi can be further clumped
early_bmi_clumped<-clump_data(early_bmi, clump_r2 = 0.01)
dim(early_bmi_clumped) # 126

length(intersect(early_bmi_clumped$SNP, instruments_EXP2$SNP)) # 121 


# run MR on all

quick_MR_w_BC<- function(instruments){
  instr_in_breast_cancer <- extract_outcome_data(
            snps = instruments$SNP,
            outcome = 'ieu-a-1126') 
  dat <- harmonise_data(exposure_dat = instruments, 
                          outcome_dat = instr_in_breast_cancer)
  res <- mr(dat, method_list=c('mr_ivw')) %>% 
        generate_odds_ratios()
  
  return(res)
}

result_OUT <- quick_MR_w_BC(instruments_OUT)
result_EXP <- quick_MR_w_BC(instruments_EXP)
result_EXP2 <- quick_MR_w_BC(instruments_EXP2)
result_early_bmi <- quick_MR_w_BC(early_bmi)
result_early_bmi_cl <- quick_MR_w_BC(early_bmi_clumped)

# join into one df and add desctiptions
dat<- bind_rows(result_OUT,
               result_EXP,
               result_EXP2,
               result_early_bmi,
               result_early_bmi_cl)
dat$source <-
    c("instr. from Outcome format",
      "instr. from Exposure format",
      "instr. from Exposure format, clumped w/ r2=0.01 ",
      "instr. from Supl Table (as is)",
      "instr. from Supl Table, clumped w/ r2=0.01")
 
View(dat)     



