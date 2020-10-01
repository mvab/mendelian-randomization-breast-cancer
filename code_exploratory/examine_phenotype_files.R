library(vroom)
library(dplyr)

path<-"/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project2/Data/phenotypes/"

okdata<-vroom(paste0(path, "testosterone2_female_invnormal_gwas.txt"), col_select=c("IID", "phenotype1")) %>% distinct()
dim(okdata)
sum(is.na(okdata$phenotype1))
okdata2<-vroom(paste0(path, "biotest2_female_invnormal_gwas.txt"), col_select=c("IID", "phenotype1")) %>% distinct() 
dim(okdata2)
sum(is.na(okdata2$phenotype1))

baddata<-vroom(paste0(path, "igf1_2_female_invnormal_gwas.txt"), col_select=c("IID", "phenotype1")) %>% distinct() 
dim(baddata)
sum(is.na(baddata$phenotype1))

bmidata<-vroom(paste0(path, "earlybmi1_female_gwas.txt"), col_select=c("IID", "phenotype1")) %>% distinct() 
dim(bmidata)
sum(is.na(bmidata$phenotype1))


length(intersect(okdata$IID, okdata2$IID))

okdata_compl<-okdata[complete.cases(okdata),]
dim(okdata_compl)#105374
okdata2_compl<-okdata2[complete.cases(okdata2),]
dim(okdata2_compl) #95267
baddata_compl<-baddata[complete.cases(baddata),]
dim(baddata_compl) #125431
bmidata_compl<-bmidata[complete.cases(bmidata),]
dim(bmidata_compl) #130013


length(intersect(okdata_compl$IID, okdata2_compl$IID)) #95267
length(intersect(baddata_compl$IID, okdata_compl$IID)) #104737
length(intersect(baddata_compl$IID, okdata2_compl$IID)) #94683

length(intersect(bmidata_compl$IID, okdata_compl$IID)) # 0
length(intersect(bmidata_compl$IID, okdata2_compl$IID)) # 0
length(intersect(bmidata_compl$IID, baddata_compl$IID)) # 0


test<-okdata[duplicated(okdata),]

m <- full_join( okdata, baddata, by="IID")
dim(m)
m2 <- full_join( okdata, baddata, by="IID")
dim(m2)
m3 <- full_join(  baddata,okdata, by="IID")
dim(m3)


l <- full_join( okdata, okdata2, by="IID")
lm<-full_join(l, baddata, by = "IID")
View(lm)

sum(is.na(okdata$IID))
sum(is.na(baddata$IID))
sum(is.na(m$IID))

sum(is.na(m$phenotype1.x))
sum(is.na(m$phenotype1.y))



### testing inclusion


okdata2<-vroom(paste0(path, "biotest2_female_invnormal_gwas.txt"), col_select=c("IID", "phenotype1")) %>% distinct() 
dim(okdata2)
sum(is.na(okdata2$phenotype1))

okdata2_compl<-okdata2[complete.cases(okdata2),]
dim(okdata2_compl) #95267

cov_file<-readr::read_delim(paste0(path,"data.covariates.bolt.wAge.txt"), delim=" ") %>% select(-FID)

test<-left_join(okdata2_compl, cov_file, by="IID")
dim(test)

test_compl<-test[complete.cases(test),]
dim(test_compl)

#### combing SHBG data into one

part1<-vroom(paste0(path, "shbg1_female_invnormal_gwas.txt")) %>% distinct()
dim(part1)
part2<-vroom(paste0(path, "shbg2_female_invnormal_gwas.txt")) %>% distinct()
dim(part2)

length(intersect(part1$IID, part2$IID))

part1_complete<-part1[complete.cases(part1),]
dim(part1_complete)
part2_complete<-part2[complete.cases(part2),]
dim(part2_complete)

length(intersect(part1_complete$IID, part2_complete$IID))

m <- full_join(part1, part2, by = c("FID", "IID"))

dim(m)

m2<- m %>% 
    mutate(phenotype1 = coalesce(phenotype1.x, phenotype1.y)) %>% 
    select(-phenotype1.x, -phenotype1.y)

m2_complete<-m2[complete.cases(m2),]
dim(m2_complete)

vroom_write(m2_complete, delim = " ", paste0(path, "shbg_all_female_invnormal_gwas.txt"))



