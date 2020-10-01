
library(data.table)

all.files.tsv <- list.files(path = data_path_tophits, pattern = "*tsv", full.names = T)

l <- lapply(all.files.tsv, fread)

get_names <- function(path_and_name){
  tmp<-path_and_name %>% 
    basename(.) %>% 
    strsplit(., ".", fixed = TRUE)
  return(tmp[[1]][1])
}


all.names <- unlist(lapply(all.files.tsv, get_names))


# update list names to file name
for (i in 1:length(l) ){
  #print(all.files[i])
  names(l)[i] <- all.names[i]
}

early_tophits <- l[[5]]$SNP


for (i in 1:length(l)){
  
  overlap <- length(intersect(early_tophits, l[[i]]$SNP))
  print(paste0("Overlap of Early BMI SNPs is ", overlap, " with ", names(l)[i]))
}
