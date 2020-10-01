# set paths

set_paths<- function(currently_working_env){
  if (currently_working_env == "remote"){
    ## remote data paths (RDSF)
    data_path_tophits <<- "../../data/GWAS_tophits/"  # or clumped_instruments
    data_path_gwas <<- "../../data/GWAS_results_tidy/" # or simply GWAS_results 
    results_path <<-  "../../results/"
  
  } else if (currently_working_env == "local"){
    ## local data paths
    local_path <<-"XremovedX"
    data_path_tophits <<- paste0(local_path, "Data/GWAS_tophits/") 
    data_path_tophits_raw <<-   paste0(data_path_tophits, "unprocessed/")
    data_path_gwas_raw <<-  paste0(local_path,"Data/GWAS_results/") 
    data_path_gwas <<-  paste0(local_path,"Data/GWAS_results_tidy/") 
    results_path <<-   paste0(local_path, "Results/")
  }
}

# << makes it global variable