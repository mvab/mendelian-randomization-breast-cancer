require("purrr")
require("tidyr")
require("tibble")
require("TwoSampleMR")

get_mv_exposures <- function(exposure_list, full_gwas_list, clump_exposures=FALSE) {
  
  # Here we are using a re-written source code of `mv_extract_exposures` function from 2SMR package.
  # It was neccesary to do it, as it only works with MRBase database input at the moment (but we have external exposures)
  
  
  # Get effects of each instrument from each exposure
  # all tophit SNPs in both exposures, (clumped optionally)
  exposures <- exposure_list %>% purrr::reduce(bind_rows) %>% flip_same_snp()

  if (clump_exposures) {
    # ***optional*** : clump exposures 
    temp <- exposures
    temp$id.exposure <- 1
    temp <- clump_data(temp)
    #temp <- clump_data_local(temp, local_path)
    exposures <- filter(exposures, SNP %in% temp$SNP)
  }
  
  
  # merge exposures (in 'outcomes' ~ full gwas format)
  # extract all instruments from full GWAS of exposures
  for (i in 1:length(full_gwas_list)){
    full_gwas_list[[i]] <- full_gwas_list[[i]] %>% filter(SNP %in% exposures$SNP)
  }
  d1 <- full_gwas_list %>%
    purrr::reduce(bind_rows) %>% 
    distinct()

  # get ids
  id_exposure <- unique(d1$id.outcome) 
  
  # convert first trait to exposure format  -- exp1 is exposure
  tmp_exposure <- d1 %>% filter(id.outcome == id_exposure[1]) %>% convert_outcome_to_exposure()
  # keep other traits as outcome -- exp2+ are outcomes
  tmp_outcome <- d1 %>% filter(id.outcome != id_exposure[1])
  
  # Harmonise against the first trait
  d <- harmonise_data(exposure_dat = tmp_exposure, 
                      outcome_dat = tmp_outcome, action=2)
  
  # Only keep SNPs that are present in all
  snps_not_in_all <- d %>% dplyr::count(SNP)  %>% 
                    filter(n < length(exposure_list)-1) %>%
                    pull(SNP)
  d <- filter(d, !SNP %in% snps_not_in_all)

  
  # Subset and concat data
  
  # for exp1 get exposure cols
  dh1x <- d %>% filter(id.outcome == id.outcome[1]) %>% 
    select(SNP, contains("exposure"))
  # for exp2 get outcome cols
  dh2x <-d %>%  select(SNP, contains("outcome"))
  # rename outcome to exposure in these
  names(dh2x) <- gsub("outcome", "exposure", names(dh2x) )
  # join together (drop not needed cols)
  exposure_dat <- bind_rows(dh1x, dh2x) %>%  
    select(-c("samplesize.exposure" ,"mr_keep.exposure", "pval_origin.exposure")) %>% 
    distinct()
  
  return(exposure_dat)
}



#  function to convert 2SMR format into MVMR format
make_mvmr_input <- function(exposure_dat, outcome.id.mrbase="", outcome.data=""){
  # provide exposure_dat created in the same way as for TwoSampleMR 
  # also specify the outcome argument [only ONE!] (MR-base ID or full gwas data in .outcome format)
  
  # extract SNPs for both exposures from outcome dataset
  # (for the selected option mr.base or local outcome data)
  if (outcome.id.mrbase != "") {
    # if mrbase.id is provided
    outcome_dat <- extract_outcome_data(snps = unique(exposure_dat$SNP),
                                        outcomes = outcome.id.mrbase)
  } else if (outcome.data != ""){
    # if outcome df is provided
    outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  }
  
  # harmonize datasets
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  # Create variables for the analysis 
  
  ### works for many exposures
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  # add beta/se names
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}


tidy_mvmr_output <- function(mvmr_res) {
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate) %>% 
    rename(se="Std. Error") %>% 
    rename(pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}






flip_same_snp <- function(dat){
  # test if there are any unharmonised dups
  x <- dat %>% dplyr::count(SNP) %>% filter(n > 1) %>% arrange(-n)
  if (dim(x)[1] !=0 ){
    # all SNPs that occur > 1
    y <- dat %>% filter(SNP %in% x$SNP) 
    dat.upd<-tibble()
    
    for ( i in unique(y$SNP) ){
      # subset to current snp
      y.sub<- filter(y, SNP == i)
      snp_n <- dim(y.sub)[1]
      
      if (snp_n == 3){
        # likely 2/3 are in sync, need to drop 1 of the 'in sync' ones
        if (y.sub$effect_allele.exposure[1] == y.sub$effect_allele.exposure[2]) {
          dat.upd<-rbind(dat.upd, y.sub[1,])
          y.sub<-y.sub[2:3,]
        } else if (y.sub$effect_allele.exposure[1] == y.sub$effect_allele.exposure[3]){
          dat.upd<-rbind(dat.upd, y.sub[3,])
          y.sub<-y.sub[1:2,]
        } else if (y.sub$effect_allele.exposure[2] == y.sub$effect_allele.exposure[3]){
          dat.upd<-rbind(dat.upd, y.sub[3,])
          y.sub<-y.sub[1:2,]
        }
      }
      
      # if there are SNPs that seem to be exact opposite, flip them
      # flip and harmonise each pair of discordant SNPs
      if (y.sub$effect_allele.exposure[1] == y.sub$other_allele.exposure[2] & y.sub$effect_allele.exposure[2] == y.sub$other_allele.exposure[1] ){
        y.sub$other_allele.exposure[2] <-  y.sub$other_allele.exposure[1]
        y.sub$effect_allele.exposure[2] <- y.sub$effect_allele.exposure[1]
        y.sub$beta.exposure[2] <- y.sub$beta.exposure[2] * -1
        y.sub$eaf.exposure[2] <- 1 - y.sub$eaf.exposure[2]
      }
      # save new (or unchanged if if was FALSE) into tibbl dat.upd
      dat.upd <- rbind(dat.upd, y.sub)
    }

    # discard old, add new harmonised data
    dat <- dat %>% filter(!SNP %in% unique(y$SNP) )
    dat <- rbind(dat, dat.upd)
    return(dat)
  } else {
    # no dups
    return(dat)
  }  
}

qhet_mvmr_tmp<- function (r_input, pcor, CI, iterations) {
  library(boot)
  warning("qhet_mvmr() is currently undergoing development.")
  warning("This is Marina's local version with simple CIs.")
  if (missing(CI)) {
    CI <- F
    warning("95 percent confidence interval not calculated")
  }
  if (missing(iterations)) {
    iterations <- 1000
    warning("Iterations for bootstrap not specified. Default = 1000")
  }
  exp.number <- length(names(r_input)[-c(1, 2, 3)])/2
  
  Qtemp <- function(r_input, pcor) {
    exp.number <- length(names(r_input)[-c(1, 2, 3)])/2
    stderr <- as.matrix(r_input[, (exp.number + 4):length(r_input)])
    correlation <- pcor
    gammahat <- r_input$betaYG
    segamma <- r_input$sebetaYG
    pihat <- as.matrix(r_input[, c(4:(3 + exp.number))])
    PL_MVMR = function(a) {
      tau2 = a[1]
      PL2_MVMR = function(ab) {
        b <- ab
        cov = matrix(nrow = exp.number, ncol = exp.number)
        w = NULL
        for (l in 1:nrow(r_input)) {
          for (pp in 1:exp.number) {
            for (p2 in 1:exp.number) {
              cov[pp, p2] <- correlation[pp, p2] * stderr[l, 
                                                          pp] * stderr[l, p2]
            }
          }
          segamma <- r_input$sebetaYG
          w[l] <- segamma[l]^2 + t(b) %*% cov %*% b + 
            tau2
        }
        q = sum((1/w) * ((gammahat - pihat %*% b)^2))
        return(q)
      }
      st_PL2 = rep(0, exp.number)
      bc = optim(st_PL2, PL2_MVMR)
      bcresults <- bc$par
      cov = matrix(nrow = exp.number, ncol = exp.number)
      w = NULL
      for (l in 1:nrow(r_input)) {
        for (pp in 1:exp.number) {
          for (p2 in 1:exp.number) {
            cov[pp, p2] <- correlation[pp, p2] * stderr[l, 
                                                        pp] * stderr[l, p2]
          }
        }
        w[l] <- segamma[l]^2 + t(bcresults) %*% cov %*% 
          bcresults + tau2
      }
      q = (sum((1/w) * ((gammahat - pihat %*% bcresults)^2)) - 
             (nrow(r_input) - 2))^2
    }
    PL2_MVMR = function(ab) {
      b = ab
      w = NULL
      cov = matrix(nrow = exp.number, ncol = exp.number)
      for (l in 1:nrow(r_input)) {
        for (pp in 1:exp.number) {
          for (p2 in 1:exp.number) {
            cov[pp, p2] <- correlation[pp, p2] * stderr[l, 
                                                        pp] * stderr[l, p2]
          }
        }
        w[l] <- segamma[l]^2 + t(b) %*% cov %*% b + 
          tau_i
      }
      q = sum((1/w) * ((gammahat - pihat %*% b)^2))
    }
    limltauest = optimize(PL_MVMR, interval = c(-10, 10))
    tau_i = limltauest$objective
    tau = tau_i
    liml_het2 <- optim(rep(0, exp.number), PL2_MVMR)
    limlhets <- liml_het2$par
    Qexact_het <- liml_het2$value
    Effects <- limlhets
    Effects <- data.frame(Effects)
    names(Effects) <- "Effect Estimates"
    for (i in 1:exp.number) {
      rownames(Effects)[i] <- paste("Exposure", i, sep = " ")
    }
    return(Effects)
  }
  if (CI == F) {
    res <- Qtemp(r_input, pcor)
  }
  if (CI == T) {
    bootse <- function(data, indices) {
      bres <- Qtemp(data[indices, ], pcor)[, 1]
      return(bres)
    }
    
    b.results <- boot(data = r_input, statistic = bootse, 
                        R = iterations)
      
      
    out<-data.frame(b=b.results$t0,
                    se=apply(b.results$t, 2, sd),
                    mean = apply(b.results$t, 2, mean))
    out$bias = out$mean - out$b
    out$mean <-NULL
    
    rownames(out) <- paste("Exposure", c(1:exp.number), sep = " ")
      
    out<-TwoSampleMR::generate_odds_ratios(out)
  }
  return(out)

}





























