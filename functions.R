
add_trait_type_exp <- function(dat, category){
  
  if (category %in% c("hormones", "hormones_splitsample")){
    dat <- mutate(dat, 
                  trait_type =ifelse(grepl("IGF", exposure), "IGF",
                                     ifelse(grepl("Oestradiol", exposure), "Oestradiol",
                                            ifelse(grepl("SHBG", exposure), "SHBG",
                                                   ifelse(grepl("Testosterone", exposure), "Testosterone", NA
                                                   )))))
    
  } else if (category == "testosterone"){
    
    dat <- mutate(dat, 
                  measure =  ifelse(grepl("free", exposure), "Free", 
                                    ifelse(grepl("bioavailable", exposure), "Bioavailable",
                                           ifelse(grepl("total", exposure), "Total", "Total (old)")))) %>% 
      mutate(menopause = ifelse(grepl("pre", exposure), "pre-menopause", 
                                ifelse(grepl("post", exposure), "post-menopause", "all"))) %>% 
      #mutate(measure = factor(measure, levels=c("Total", "Total (old)", "Free", "Bioavailable"))) %>% 
      mutate(menopause = factor(menopause, levels=c( "pre-menopause", "post-menopause", "all"))) %>% 
      arrange(menopause)
    
    # creatting dummy data rows fro total,old, all
    tmp<-dat %>% filter(measure=="Total" & menopause == "all") %>% 
      mutate(or=NA, or_lci95=NA , or_uci95=NA) %>% 
      mutate(measure = "Total (old)" )
    dat <- rbind(dat, tmp)
    
  } else if (category == "reproductive_traits"){
    dat <- mutate(dat, 
                  trait_type = ifelse(grepl("menarche", exposure), "menarche",
                                      ifelse(grepl("menopause", exposure), "menopause",
                                             ifelse(grepl("first live birth", exposure), "first birth",
                                                    ifelse(grepl("Number of", exposure), "births no.", NA
                                                    )))))
  } else if (category == "glycemic_traits"){
    
    dat <- mutate(dat, 
                  trait_type = ifelse(grepl("insulin", exposure), "insulin",
                                      ifelse(grepl("glucose", exposure), "glucose",
                                             ifelse(grepl("HOMA-B", exposure), "HOMA-B",
                                                    ifelse(grepl("HOMA-IR", exposure), "HOMA-IR",
                                                           ifelse(grepl("HBa1c", exposure, ignore.case = T), "HBa1c", NA
                                                           ))))))
  } else if (category == "physical_traits"){
    dat <- mutate(dat, 
                  trait_type = ifelse(grepl("Childhood height", exposure), "Childhood height",
                                  ifelse(grepl("Breast", exposure), "breast",
                                      ifelse(grepl("Percent", exposure), "percent",
                                             ifelse(grepl("non-dense", exposure), "non-dense",
                                                    ifelse(grepl("dense", exposure), "dense", NA
                                                    ))))))
  } else{
    print("unknown categogy; returning as is")
  }
  return(dat)
}


add_trait_type_out <- function(dat, category){
  
  if (category %in% c("hormones", "hormones_splitsample")){
    dat <- mutate(dat, 
                  trait_type =ifelse(grepl("IGF", outcome), "IGF",
                                     ifelse(grepl("Oestradiol", outcome), "Oestradiol",
                                            ifelse(grepl("SHBG", outcome), "SHBG",
                                                   ifelse(grepl("Testosterone", outcome), "Testosterone", NA
                                                   )))))
    
  } else if (category == "testosterone"){
    
    dat <- mutate(dat, 
                  measure =  ifelse(grepl("free", outcome), "Free", 
                                    ifelse(grepl("bioavailable", outcome), "Bioavailable",
                                           ifelse(grepl("total", outcome), "Total", "Total (old)")))) %>% 
      mutate(menopause = ifelse(grepl("pre", outcome), "pre-menopause", 
                                ifelse(grepl("post", outcome), "post-menopause", "all"))) %>% 
      #mutate(measure = factor(measure, levels=c("Total", "Total (old)", "Free", "Bioavailable"))) %>% 
      mutate(menopause = factor(menopause, levels=c( "pre-menopause", "post-menopause", "all"))) %>% 
      arrange(menopause)
    
    # creatting dummy data rows fro total,old, all
    tmp<-dat %>% filter(measure=="Total" & menopause == "all") %>% 
      mutate(b=NA, lo_ci=NA , up_ci=NA) %>% 
      mutate(measure = "Total (old)" )
    dat <- rbind(dat, tmp)
    
  } else if (category == "reproductive_traits"){
    dat <- mutate(dat, 
                  trait_type = ifelse(grepl("menarche", outcome), "menarche",
                                      ifelse(grepl("menopause", outcome), "menopause",
                                             ifelse(grepl("first live birth", outcome), "first birth",
                                                    ifelse(grepl("Number of", outcome), "births no.", NA
                                                    )))))
  } else if (category == "glycemic_traits"){
    
    dat <- mutate(dat, 
                  trait_type = ifelse(grepl("insulin", outcome), "insulin",
                                      ifelse(grepl("glucose", outcome), "glucose",
                                             ifelse(grepl("HOMA-B", outcome), "HOMA-B",
                                                    ifelse(grepl("HOMA-IR", outcome), "HOMA-IR",
                                                           ifelse(grepl("HBa1c", outcome, ignore.case = T), "HBa1c", NA
                                                           ))))))
  } else if (category == "physical_traits"){
    
    dat <- mutate(dat, 
                  trait_type = ifelse(grepl("Childhood height", outcome), "Childhood height",
                                      ))
  }
  return(dat)
}

# supl function to display the results as a table to make it easy to copy values

kable_it<-function(df){
  library(kableExtra)
  df %>% 
    tidy_pvals %>% 
    kable(.) %>%
    kable_styling()
}
#dat %>% kable_it()


tidy_pvals<-function(df){
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=2) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}


clump_data_local <- function(dat, path){
#https://github.com/MRCIEU/TwoSampleMR/issues/173  
  dat %>% 
  rename(rsid = SNP, 
         pval = pval.exposure,
         id = id.exposure) %>% 
  ieugwasr::ld_clump(
    dat = .,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = paste0(path, "Data/reference/1kg.v3/EUR")) %>% 
  rename(SNP = rsid, 
         pval.exposure = pval,
         id.exposure = id) 
  
}






