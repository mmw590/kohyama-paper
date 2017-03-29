#' ---
#' title: "R_Kohyama_4.ModelTesting_Seed_LFDP.R"
#' author: "Maria Wang"
#' date: "2017-03-28"
#' description: "This script tests for effects of seed mass and species tree height on the seed-to-seedling transition (STS)"
#' ---

#### Steps to install package glmmADMB (only need to do once) --
# install.packages("R2admb")
# install.packages("glmmADMB", 
#                   repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),
#                   type="source")

# ~ Loading required R packages ----
library(dplyr)
library(lme4)
library(glmmADMB)

# ~ Set working directory ----
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis")
getwd()

# ~ Create folders to store output ----
dir.create("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/figures.seed/")

# ~ Load required input files----
sshm2 <- read.csv("intermediate data/ldfp_sts_forglm.csv")
str(sshm2)
# sts = mean per-seed-success = n.new.sling/nseeds.x6


#### Correlation Between Seed Mass And Tree Height ####
cor(sshm2$seed.mass.z, sshm2$max.ht.z) #0.68

# ~ Partial correlation using ppcor----
library(ppcor)
# partial correlation between seed mass and STS given tree height
pcor.test(sshm2$seed.mass.z, sshm2$sts, sshm2$max.ht.z) 
# estimate      p.value statistic    n gp  Method
# 1 0.2018229 1.749778e-52  15.40659 5593  1 pearson
 
# partial correlation between tree height and STS given seed mass
pcor.test(sshm2$max.ht.z, sshm2$sts, sshm2$seed.mass.z) 
# estimate    p.value statistic    n gp  Method
# 1 -0.03147362 0.01859042 -2.354332 5593  1 pearson

# partial correlation between seed mass and tree height given STS
pcor.test(sshm2$seed.mass, sshm2$max.ht, sshm2$sts) 
# estimate p.value statistic    n gp  Method
# 1 0.6684136       0  67.18942 5593  1 pearson


#### Some exploratory plots ####
plot(sts ~ max.ht, data=sshm2)
plot(sts ~ seed.mass, data=sshm2)
plot(log1p(sts) ~ max.ht, data=sshm2)
plot(log1p(sts) ~ log10(seed.mass), data=sshm2)
plot(log10(sts) ~ log10(seed.mass), data=sshm2) #as in King et al 2006
plot(sts ~ seed.mass, data=sshm2, log="xy")
hist(sshm2$sts, breaks=30)


#### Hypothesis: Seed-to-Seedling (STS) ratio ####
# Alt: STS ratio is affected by max tree height and seed mass
# Null: STS ratio is not affected by max tree height and seed mass

# ~ Convert categories to factor ----
sshm2$sp <- as.factor(sshm2$sp)
sshm2$cohort <- as.factor(sshm2$cohort)
sshm2$plotno <- as.factor(sshm2$plotno)

# ~ Create 2-column matrix of success and failures for GLMM ----
seedsv <- cbind(sshm2$n.new.sling, sshm2$failures)

# ~ Read this on p-values in mixed models ----
help("pvalues",package="lme4")


#### (Initial Model-not in pub) Generalized Linear Mixed Model with Binomial Errors ####
# results not presented in paper because model is overdispersed.

#### ~ Testing interaction between seed mass and height ----
# use standardized data only for interaction model
glmm01 <- glmer(seedsv ~ seed.mass.z * max.ht.z + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial) #might take a minute
glmm02 <- glmer(seedsv ~ seed.mass.z + max.ht.z + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)

# use likelihood ratio test to compare models
anova(glmm01, glmm02) #interaction not significant, p=0.5793

# use Parametric Bootstrap to compare models (takes a long time) #not run
# library(afex)
# PBmodcomp(glmm01, glmm02, nsim=1000)

#### ~ Testing fixed effects ----
glmm12 <- glmer(seedsv ~ seed.mass + max.ht + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)
glmm13 <- glmer(seedsv ~ seed.mass + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)
glmm14 <- glmer(seedsv ~ max.ht + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)
anova(glmm12, glmm13) #max.ht not important after controlling for seedmass (p=0.6374)
anova(glmm12, glmm14) #seed mass important after controlling for hmax (p=0.00725)


#### ~ Overdispersion test ####
# if much larger than 1, the model may be overdispersed, esp if there are not many observations...
deviance(glmm14)/df.residual(glmm14) # = 2.12... likely overdispersed
sum(residuals(glmm14, type = "deviance")^2) # = 11853.16

#overdispersion function from http://glmm.wikidot.com/faq
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(glmm14) #p<0.00, this model is overdispersed.

###### (Final Model for Paper) Negative binomial glmm with offset ------- ######
# From Muscarella et al. 2012, pg.5:
# To evaluate specific factors associated with the STS transition (Question 2), we fit statistical models where the response variable was the number of seedlings recruited in individual seedling plots. 
# The log of the number of seeds observed in each associated nearby seed basket was included as an offset. 
# Initial model residuals exhibited over-dispersion so the results reported here are based on a generalized linear-mixed model with negative binomial errors.
?glmmadmb
str(sshm2)

# !!!! This takes a while (a few days?) to run so make sure you're prepared ####
glmm.nb0<- glmmadmb(n.new.sling ~ seed.mass.z * max.ht.z + offset(log(sshm2$nseeds)) + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, zeroInflation=FALSE, family="nbinom") 
glmm.nb2<- glmmadmb(n.new.sling ~ seed.mass.z + max.ht.z + offset(log(sshm2$nseeds)) + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, zeroInflation=FALSE, family="nbinom") 
glmm.nb3<- glmmadmb(n.new.sling ~ seed.mass.z + offset(log(sshm2$nseeds)) + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, zeroInflation=FALSE, family="nbinom") 
glmm.nb4<- glmmadmb(n.new.sling ~ max.ht.z + offset(log(sshm2$nseeds)) + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, zeroInflation=FALSE, family="nbinom") 
AIC(glmm.nb0, glmm.nb2, glmm.nb3, glmm.nb4)

# Results from 2016-12-16:
# df     AIC
# glmm.nb0  8 9803.46
# glmm.nb2  7 9801.92
# glmm.nb3  6 9799.98 #best model is where fixed effect is seed mass alone
# glmm.nb4  6 9806.78


anova(glmm.nb0, glmm.nb2)
# Results from 2016-12-16:
# Analysis of Deviance Table
# Model 1: n.new.sling ~ seed.mass.z + max.ht.z + offset(log(sshm2$nseeds))
# Model 2: n.new.sling ~ seed.mass.z * max.ht.z + offset(log(sshm2$nseeds))
# NoPar  LogLik Df Deviance Pr(>Chi)
# 1     7 -4894.0                     
# 2     8 -4893.7  1     0.46   0.4976

anova(glmm.nb2, glmm.nb3)
# Results from 2016-12-16:
# Analysis of Deviance Table
# Model 1: n.new.sling ~ seed.mass.z + offset(log(sshm2$nseeds))
# Model 2: n.new.sling ~ seed.mass.z + max.ht.z + offset(log(sshm2$nseeds))
# NoPar LogLik Df Deviance Pr(>Chi)
# 1     6  -4894                     
# 2     7  -4894  1     0.06   0.8065


anova(glmm.nb2, glmm.nb4)
# Results from 2016-12-16:
# Analysis of Deviance Table
# Model 1: n.new.sling ~ max.ht.z + offset(log(sshm2$nseeds))
# Model 2: n.new.sling ~ seed.mass.z + max.ht.z + offset(log(sshm2$nseeds))
# NoPar  LogLik Df Deviance Pr(>Chi)   
# 1     6 -4897.4                        
# 2     7 -4894.0  1     6.86 0.008815 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# RESULT: SEED MASS IS SIGNIFICANT FIXED EFFECT!

range(log(sshm2$nseeds))

##### Interpretation of Results ####
# Mean per seed success or STS ratio increases with seed mass, but not tree height
# Interaction between seed mass and tree height not significant


#### Visualizing Results  ####
library(ggplot2)
library(scales)

#### ~ Fig 1(a) Actual seed survival ~ seed mass + shaded by hmax FOR PUBLICATION ####
sshm2 <- read.csv("intermediate data/ldfp_sts_forglm.csv")
sshm2.sp <- sshm2 %>% group_by(sp, seed.mass, max.ht) %>% summarize(n.new.sling=sum(n.new.sling), nseeds = sum(nseeds), sts.mean=mean(sts))

range(log10(sshm2.sp$sts.mean))
table(sshm2.sp$sts.mean) #the value closest to zero is 0.00016175748513507

sshm2.sp[sshm2.sp$sts.mean > 1, ]$sts.mean <- 1 #max sts ratio is 1
sshm2.sp[sshm2.sp$sts.mean == 0, ]$sts.mean <- 1*10^-4 #so that log10(0) doesn't turn the value into -Inf

a3 <- ggplot(sshm2.sp, aes(x=seed.mass, y=sts.mean, fill=max.ht)) + geom_point(aes(x=seed.mass, y=sts.mean), data=sshm2.sp, shape=21, size=3) + ylab("STS ratio") + xlab("Seed mass (g)") 

a31 <- a3 + scale_x_log10(breaks = c(0, 0.001, 0.01,0.1,1), labels=c(NA, "0.001", "0.01", "0.1", "1")) + 
  scale_y_log10(breaks = c(0, 10^-4, 10^-2,1), labels=c(NA, "0","0.01", "1")) + #the labels are such that log(0.00001) (actual sts value is zero) can be plotted and displayed as "0",
  theme_classic(base_size = 9) + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), #line size
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"), #adjust margins
        legend.position="none", #no legends
        plot.title=element_text(hjust=0, size=9) ) + #plot.title for labelling subplots
  scale_fill_gradient(low=gray(1), high =gray(0)) + 
  annotation_logticks() +
  ggtitle('(a)')

png(file="figures.seed/lfdp_sts-seedmass_avg_shadedbyheight.png",width=82,height=82, units="mm", res=300)
a31
dev.off()


#
#
#

#### ***Can ignore the rest of code as not relevant to current analysis*** ####

##### Model evaluation and Diagnostic plots (Initial Model-not in pub) ####

# ~ Plot random effects against the predicted values from the fixed effect component of the model and check for no trend ----
sshm2$ft.fix   <- as.numeric(model.matrix(glmm14) %*% fixef(glmm14))

# random effects
rr.sp          <- ranef(glmm14, drop = TRUE)$sp 
rr.cohort      <- ranef(glmm14, drop = TRUE)$cohort 
rr.plotno      <- ranef(glmm14, drop = TRUE)$plotno #

# create df with random effect levels and values
rrr.sp    <- data.frame(sp = attr(rr.sp, "names"), ranef.sp = rr.sp)
rrr.cohort<- data.frame(cohort = attr(rr.cohort, "names"), ranef.cohort = rr.cohort) 
rrr.plotno<- data.frame(plotno = attr(rr.plotno, "names"), ranef.plotno = rr.plotno) 
sshm2$cohort <- as.factor(sshm2$cohort)
sshm2$plotno <- as.factor(sshm2$plotno)
sshm5 <- left_join(sshm2, rrr.sp) %>% left_join(rrr.cohort) %>% left_join(rrr.plotno) %>% arrange(sp, cohort, plotno)

# ~ The main diagnostic plot should be estimated random effects vs. fitted values from a GLMM fit.
png("figures.seed/modeleval-glmm14-ranef.vs.fitted.png", width = 25, height = 12, units = "cm", res=150)
par(mfcol=c(1,3))
with(sshm5, plot(ft.fix, ranef.sp, pch = 19, las = 1, cex = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
with(sshm5, plot(ft.fix, ranef.cohort, pch = 19, las = 1, cex = 1, main="ranef vs fitted - cohort"))
abline(0,0,lwd = 1.5)
with(sshm5, plot(ft.fix, ranef.plotno, pch = 19, las = 1, cex = 1, main="ranef vs fitted - plotno"))
abline(0,0,lwd = 1.5)
dev.off()

png("figures.seed/modeleval-glmm14-ranef.vs.fitted-sp.png", width = 18, height = 12, units = "cm", res=150)
with(sshm5, plot(ft.fix, ranef.sp, pch = 19, las = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
text(sshm5$ft.fix, sshm5$ranef.sp, labels=sshm5$sp, pos=4, cex = 0.5)
dev.off()

png("figures.seed/modeleval-glmm14-ranef.vs.fitted-cohort.png", width = 5, height = 5, units = "in", res=220)
with(sshm5, plot(ft.fix, ranef.cohort, pch = 19, las = 1, main="ranef vs fitted - cohort"))
abline(0,0,lwd = 1.5)
text(sshm5$ft.fix, sshm5$ranef.cohort, labels=unique(sshm5$cohort), pos=4, cex=0.5)
identify(sshm5$ft.fix, sshm5$ranef.cohort, labels=sshm5$cohort) #you can click on points on the plot to identify and label them
dev.off()

### Use identify tool to identify outlier points ###
x11()
with(sshm5, plot(ft.fix, ranef.plotno, pch = 19, las = 1, main="ranef vs fitted - plotno"))
abline(0,0,lwd = 1.5)
identify(sshm5$ft.fix, sshm5$ranef.plotno, labels=sshm5$sp)
dev.off()

x11()
with(sshm5, plot(ft.fix, ranef.plotno, pch = 19, las = 1, main="ranef vs fitted - plotno"))
abline(0,0,lwd = 1.5)
identify(sshm5$ft.fix, sshm5$ranef.plotno, labels=sshm5$sp)
dev.off()

#write.csv(sshm5, "sshm5.csv")


#### ~ Check for approximate normality of random effects: ####
png("figures.seed/modeleval-glmm14-qqplot-ranef.png", width = 24, height = 8, units = "cm", res=150)
par(mfcol=c(1,3))
qqnorm(ranef(glmm14, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(glmm)$sp")
qqline(ranef(glmm14, drop = TRUE)$sp, las = 1, cex = 1.4) #steps?!

qqnorm(ranef(glmm14, drop = TRUE)$cohort, las = 1, cex = 1.4, main="normality of ranef(glmm)$cohort")
qqline(ranef(glmm14, drop = TRUE)$cohort, las = 1, cex = 1.4) #seems fine except for the last year?

qqnorm(ranef(glmm14, drop = TRUE)$plotno, las = 1, cex = 1.4, main="normality of ranef(glmm)$plotno")
qqline(ranef(glmm14, drop = TRUE)$plotno, las = 1, cex = 1.4) #seems fine except for the last plot?
dev.off()

# use identify to label outliers
par(mfcol=c(1,1))
qq_sp <- qqnorm(ranef(glmm14, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(glmm)$sp")
identify(qq_sp)
qqline(ranef(glmm14, drop = TRUE)$sp, las = 1, cex = 1.4)
unique(sshm5$sp)[c(15,19) ] #outliers are GUAGUI, LAEPRO

qq_coh <- qqnorm(ranef(glmm14, drop = TRUE)$cohort, las = 1, cex = 1.4, main="normality of ranef(glmm)$cohort")
identify(qq_coh)
qqline(ranef(glmm14, drop = TRUE)$cohort, las = 1, cex = 1.4)
levels(sshm5$cohort)[c(4,7)] #2014 is very deviant year. drought?
ranef(glmm14)$cohort

qq_pno <- qqnorm(ranef(glmm14, drop = TRUE)$plotno, las = 1, cex = 1.4, main="normality of ranef(glmm)$plotno")
qqline(ranef(glmm14, drop = TRUE)$plotno, las = 1, cex = 1.4)
identify(qq_pno) #89 (far top)
unique(sshm5$plotno)[c(89)] #plot 106 is weird (highest ranef 1.95)

dev.off()

#qqplot for cohort ranef #
png("figures.seed/modeleval-glmm14-qqplot-ranef_cohort.png", width = 6, height = 5, units = "in", res=300)
qq_coh <- qqnorm(ranef(glmm14, drop = TRUE)$cohort, las = 1, cex = 1.4, main="normality of ranef(glmm)$cohort", ylim=c(-6,4))
qqline(ranef(glmm14, drop = TRUE)$cohort, las = 1, cex = 1.4)
text(qq_coh, labels=levels(sshm5$cohort), cex = 1, pos=3)
dev.off()


#### ~ Overdispersion test ####
# if much larger than one, it may be overdispersed, esp if there are not many observations...
deviance(glmm14)/df.residual(glmm14) # = ... 
sum(residuals(glmm14, type = "deviance")^2)

#overdispersion function from http://glmm.wikidot.com/faq
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(glmm14) #this model is overdispersed.


#### Possible issues ####
# overdispersion due to zero-inflation?
# non-normality in species random effects?
