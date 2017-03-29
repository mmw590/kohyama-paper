#### R_Kohyama_3.DataPrep_SaplingRecruitment_BCI.R ####
#### Kohyama Hypothesis Project
#### Code by Maria Wang, Summer/Fall 2016

#### Load packages ####
library(dplyr)
library(tidyr)
library(date)
attach('C:/Users/wang/Dropbox/LFDP/data analysis for Seth/CTFSRPackage.rdata')
ls(2)

#### Set working directory ####
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/")
getwd()

#### Load tree census data ####
load('input data/bci/bci.stem.Rdata31Aug2012/bci.stem1.Rdata') 
load('input data/bci/bci.stem.Rdata31Aug2012/bci.stem2.Rdata') 
load('input data/bci/bci.stem.Rdata31Aug2012/bci.stem3.Rdata') 
load('input data/bci/bci.stem.Rdata31Aug2012/bci.stem4.Rdata') 
load('input data/bci/bci.stem.Rdata31Aug2012/bci.stem5.Rdata') 
load('input data/bci/bci.stem.Rdata31Aug2012/bci.stem6.Rdata') 
load('input data/bci/bci.stem.Rdata31Aug2012/bci.stem7.Rdata') 
str(bci.stem1)
unique(bci.stem1$sp)

# checking the R table
bci.anti <- anti_join(bci.stem1, bci.stem7, by= c("treeID", "stemID", "tag", "StemTag", "sp", "quadrat", "gx", "gy")) #zero obs is good


#### Calculate recruitment using CTFS R package function ####
brec1 <- recruitment.eachspp(bci.stem1, bci.stem2, mindbh=15) %>% as.data.frame 
brec2 <- recruitment.eachspp(bci.stem2, bci.stem3, mindbh=15) %>% as.data.frame 
brec3 <- recruitment.eachspp(bci.stem3, bci.stem4, mindbh=15) %>% as.data.frame 
brec4 <- recruitment.eachspp(bci.stem4, bci.stem5, mindbh=15) %>% as.data.frame 
brec5 <- recruitment.eachspp(bci.stem5, bci.stem6, mindbh=15) %>% as.data.frame 
brec6 <- recruitment.eachspp(bci.stem6, bci.stem7, mindbh=15) %>% as.data.frame 

## Rename column names
brec1$sp <- row.names(brec1)
brec2$sp <- row.names(brec2)
brec3$sp <- row.names(brec3)
brec4$sp <- row.names(brec4)
brec5$sp <- row.names(brec5)
brec6$sp <- row.names(brec6)

### Combine dataframes ####
brec.all <- bind_rows(brec1, brec2, brec3, brec4, brec5, brec6)
brec.all$interval <- rep(1:6, each=nrow(brec1)) #create column indicating interval no.
names(brec.all)

## Clean up 
brec.all2 <- brec.all[is.na(brec.all$rate) == FALSE, ] #omit 121 NAs
range(brec.all2$rate) #Inf will cause issues when log-transforming later
brec.all2 <- brec.all2[brec.all2$rate < Inf, ] #omit Inf

## Only incld spp where number of individuals alive and >15mm dbh in census 2 is at least 10 
brec.all2 <- brec.all2[brec.all2$N2 > 9, ] 
brec.all0 <- brec.all2[brec.all2$rate == 0, ] #20 zeros


write.csv(brec.all2, "output.sapling/bci_recruit_spp.csv", row.names=FALSE)

#### Load Height data ####
bciht <- read.csv("input data/bci/bci_heights_comb.csv")
str(bciht)
bciht2 <- bciht %>% group_by(SP6) %>% summarize(max.ht=mean(HT_M, na.rm=TRUE)) %>% rename(sp=SP6)
unique(bciht2$sp)
write.csv(bciht2, "input data/bci/bci_heights_sp.csv", row.names=FALSE)


#### Preparing data for model testing ####
brec.all2 <- read.csv("output.sapling/bci_recruit_spp.csv")
bciht2 <- read.csv("input data/bci/bci_heights_sp.csv") 
str(brec.all2)
str(bciht2)
unique(brec.all2$sp); unique(bciht2$sp)

#change lowercase letters to uppercase to match treeht sp
brec.all2 <- mutate_each(brec.all2, funs(toupper), sp)

#check what tree heights are missing
bneed_heights <- anti_join(brec.all2, bciht2) %>% distinct(sp)
write.csv(bneed_heights, "output.sapling/need_heights_bci.csv")

# join dataframes
brecruit <- inner_join(brec.all2, bciht2) #only keep sp with max heights
str(brecruit)
head(brecruit)

# write csv ####
write.csv(brecruit, "output.sapling/bci_recruit_lme.csv", row.names=FALSE)
#brecruit <- read.csv("output.sapling/bci_recruit_lme.csv")

# checking
brecruit[is.na(brecruit$rate) == TRUE, ] #zero is good
range(brecruit$rate) #no Inf