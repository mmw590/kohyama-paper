#### R_Kohyama_4.ModelTesting_Seed_LFDP.R ####
#### Kohyama Hypothesis Project
#### Code by Maria Wang, Summer/Fall 2016

## Steps to install package glmmADMB properly (only need to do once) ####
# install.packages("R2admb")
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),
#                  type="source")

## Load packages ####
library(dplyr)
library(lme4)
library(glmmADMB)

#### Set working directory and load csv ####
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis")
getwd()

sshm2 <- read.csv("output.seed/seedsurvival_forglm.csv")
str(sshm2)

#### Convert categories to factor ####
sshm2$sp <- as.factor(sshm2$sp)
sshm2$cohort <- as.factor(sshm2$cohort)
sshm2$plotno <- as.factor(sshm2$plotno)

# Create 2-column matrix of success and failures for GLMM #
seedsv <- cbind(sshm2$n.new.sling, sshm2$failures)

# Correlation between seed mass and tree height ####
cor(sshm2$seed.mass, sshm2$max.ht) #0.68

# partial correlation for seed mass and seed survival given tree height  ####
library(ppcor)
pcor.test(sshm2$seed.mass, sshm2$seedsurvival, sshm2$max.ht) 
#estimate      p.value statistic    n gp  Method
#1 0.1332067 1.470499e-23  10.04893 5593  1 pearson

pcor.test(sshm2$max.ht, sshm2$seedsurvival, sshm2$seed.mass) 
#estimate   p.value statistic    n gp  Method
#1 -0.01701328 0.2033536 -1.272204 5593  1 pearson

pcor.test(sshm2$seed.mass, sshm2$max.ht, sshm2$seedsurvival) 
#estimate p.value statistic    n gp  Method
#1 0.6724488       0  67.92805 5593  1 pearson

#### Some exploratory plots ####
plot(seedsurvival ~ max.ht, data=sshm2)
plot(seedsurvival ~ seed.mass, data=sshm2)
plot(log1p(seedsurvival) ~ max.ht, data=sshm2)
plot(log1p(seedsurvival) ~ log10(seed.mass), data=sshm2)
plot(log10(seedsurvival) ~ log10(seed.mass), data=sshm2) #as in King et al 2006
plot(seedsurvival ~ seed.mass, data=sshm2, log="xy")
hist(sshm2$seedsurvival, breaks=30)

#### Hypothesis: Seed to Seedling ratio ####
# Alt: STS ratio is affected by max tree height and seed mass
# Null: STS ratio is not affected by max tree height and seed mass

# Read this on p-values in mixed models ####
help("pvalues",package="lme4")

#### Testing interaction between seed mass and height ####
# use standardized data only for interaction model
glmm01 <- glmer(seedsv ~ seed.mass.z * max.ht.z + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial) #might take a minute
glmm02 <- glmer(seedsv ~ seed.mass.z + max.ht.z + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)

# use likelihood ratio test to compare models
anova(glmm01, glmm02) #interaction not significant, p=0.579

# use Parametric Bootstrap to compare models (takes a long time) #not run
# library(afex)
# PBmodcomp(glmm01, glmm02, nsim=1000)

#### Testing fixed effects ####
glmm12 <- glmer(seedsv ~ seed.mass + max.ht + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)
glmm13 <- glmer(seedsv ~ seed.mass + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)
glmm14 <- glmer(seedsv ~ max.ht + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)
anova(glmm12, glmm13) #max.ht not important after controlling for seedmass (p=0.6374)
anova(glmm12, glmm14) #seed mass important after controlling for hmax (p=0.007)


##### Overdispersion test ####
# if much larger than 1, the model may be overdispersed, esp if there are not many observations...
deviance(glmm14)/df.residual(glmm14) # = 2.12... is it overdispersed?
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

overdisp_fun(glmm14) #p<0.05, this model is overdispersed.

#### Negative binomial glmm with offset ####
# From Muscarella et al. 2012, pg.5:
# To evaluate specific factors associated with the STS transition (Question 2), we fit statistical models where the response variable was the number of seedlings recruited in individual seedling plots. 
# The log of the number of seeds observed in each associated nearby seed basket was included as an offset. 
# Initial model residuals exhibited over-dispersion so the results reported here are based on a generalized linear-mixed model with negative binomial errors.
?glmmadmb
str(sshm2)
range(log(sshm2$totseeds))
ptm <- proc.time
glmm.nb0<- glmmadmb(n.new.sling ~ seed.mass.z * max.ht.z + offset(log(sshm2$totseeds)) + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, zeroInflation=FALSE, family="nbinom") 
glmm.nb2<- glmmadmb(n.new.sling ~ seed.mass.z + max.ht.z + offset(log(sshm2$totseeds)) + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, zeroInflation=FALSE, family="nbinom") 
glmm.nb3<- glmmadmb(n.new.sling ~ seed.mass.z + offset(log(sshm2$totseeds)) + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, zeroInflation=FALSE, family="nbinom") 
glmm.nb4<- glmmadmb(n.new.sling ~ max.ht.z + offset(log(sshm2$totseeds)) + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, zeroInflation=FALSE, family="nbinom") 
AIC(glmm.nb0, glmm.nb2, glmm.nb3, glmm.nb4)
proc.time() - ptm

# df     AIC
# glmm.nb0  8 9803.46
# glmm.nb2  7 9801.92
# glmm.nb3  6 9799.98 #best
# glmm.nb4  6 9806.78

anova(glmm.nb2, glmm.nb0)
# Analysis of Deviance Table
# Model 1: n.new.sling ~ seed.mass.z + max.ht.z + offset(log(sshm2$totseeds))
# Model 2: n.new.sling ~ seed.mass.z * max.ht.z + offset(log(sshm2$totseeds))
# NoPar  LogLik Df Deviance Pr(>Chi)
# 1     7 -4894.0                     
# 2     8 -4893.7  1     0.46   0.4976

anova(glmm.nb3, glmm.nb2)
# Analysis of Deviance Table
# Model 1: n.new.sling ~ seed.mass.z + offset(log(sshm2$totseeds))
# Model 2: n.new.sling ~ seed.mass.z + max.ht.z + offset(log(sshm2$totseeds))
# NoPar LogLik Df Deviance Pr(>Chi)
# 1     6  -4894                     
# 2     7  -4894  1     0.06   0.8065

anova(glmm.nb4, glmm.nb2)
# Analysis of Deviance Table
# Model 1: n.new.sling ~ max.ht.z + offset(log(sshm2$totseeds))
# Model 2: n.new.sling ~ seed.mass.z + max.ht.z + offset(log(sshm2$totseeds))
# NoPar  LogLik Df Deviance Pr(>Chi)   
# 1     6 -4897.4                        
# 2     7 -4894.0  1     6.86 0.008815 **


##### Model evaluation and Diagnostic plots ####

# Plot random effects against the predicted values from the fixed effect component of the model and check for no trend: 
sshm2$ft.fix   <- as.numeric(model.matrix(glmm.nb3) %*% fixef(glmm.nb3))

# random effects
rr.sp          <- ranef(glmm.nb3, drop = TRUE)$sp 
rr.cohort      <- ranef(glmm.nb3, drop = TRUE)$cohort 
rr.plotno      <- ranef(glmm.nb3, drop = TRUE)$plotno #

# create df with random effect levels and values
rrr.sp    <- data.frame(sp = row.names(rr.sp), ranef.sp = rr.sp[,1])
rrr.cohort<- data.frame(cohort = row.names(rr.cohort), ranef.cohort = rr.cohort[,1]) 
rrr.plotno<- data.frame(plotno = row.names(rr.plotno), ranef.plotno = rr.plotno[,1]) 
sshm2$cohort <- as.factor(sshm2$cohort)
sshm2$plotno <- as.factor(sshm2$plotno)
sshm5 <- left_join(sshm2, rrr.sp) %>% left_join(rrr.cohort) %>% left_join(rrr.plotno) %>% arrange(sp, cohort, plotno)

# The main diagnostic plot should be estimated random effects vs. fitted values from a GLMM fit
png("figures.seed/modeleval-glmm.nb3-ranef.vs.fitted.png", width = 25, height = 12, units = "cm", res=150)
par(mfcol=c(1,3))
with(sshm5, plot(ft.fix, ranef.sp, pch = 19, las = 1, cex = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
with(sshm5, plot(ft.fix, ranef.cohort, pch = 19, las = 1, cex = 1, main="ranef vs fitted - cohort"))
abline(0,0,lwd = 1.5)
with(sshm5, plot(ft.fix, ranef.plotno, pch = 19, las = 1, cex = 1, main="ranef vs fitted - plotno"))
abline(0,0,lwd = 1.5)
dev.off()

png("figures.seed/modeleval-glmm.nb3-ranef.vs.fitted-sp.png", width = 18, height = 12, units = "cm", res=150)
with(sshm5, plot(ft.fix, ranef.sp, pch = 19, las = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
text(sshm5$ft.fix, sshm5$ranef.sp, labels=sshm5$sp, pos=4, cex = 0.5)
dev.off()

png("figures.seed/modeleval-glmm.nb3-ranef.vs.fitted-cohort.png", width = 5, height = 5, units = "in", res=220)
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


## now check for approximate normality of random effects: ####
png("figures.seed/modeleval-glmm.nb3-qqplot-ranef.png", width = 24, height = 8, units = "cm", res=150)
par(mfcol=c(1,3))
qqnorm(ranef(glmm.nb3, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(glmm)$sp")
qqline(ranef(glmm.nb3, drop = TRUE)$sp, las = 1, cex = 1.4) #steps?!

qqnorm(ranef(glmm.nb3, drop = TRUE)$cohort, las = 1, cex = 1.4, main="normality of ranef(glmm)$cohort")
qqline(ranef(glmm.nb3, drop = TRUE)$cohort, las = 1, cex = 1.4) #seems fine except for the last year?

qqnorm(ranef(glmm.nb3, drop = TRUE)$plotno, las = 1, cex = 1.4, main="normality of ranef(glmm)$plotno")
qqline(ranef(glmm.nb3, drop = TRUE)$plotno, las = 1, cex = 1.4) #seems fine except for the last plot?
dev.off()

# use identify to label outliers
par(mfcol=c(1,1))
qq_sp <- qqnorm(ranef(glmm.nb3, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(glmm)$sp")
identify(qq_sp)
qqline(ranef(glmm.nb3, drop = TRUE)$sp, las = 1, cex = 1.4)
unique(sshm5$sp)[c(15,19) ] #outliers are GUAGUI, LAEPRO

qq_coh <- qqnorm(ranef(glmm.nb3, drop = TRUE)$cohort, las = 1, cex = 1.4, main="normality of ranef(glmm)$cohort")
identify(qq_coh)
qqline(ranef(glmm.nb3, drop = TRUE)$cohort, las = 1, cex = 1.4)
levels(sshm5$cohort)[c(4,7)] #2014 is very deviant year. drought?
ranef(glmm.nb3)$cohort

qq_pno <- qqnorm(ranef(glmm.nb3, drop = TRUE)$plotno, las = 1, cex = 1.4, main="normality of ranef(glmm)$plotno")
qqline(ranef(glmm.nb3, drop = TRUE)$plotno, las = 1, cex = 1.4)
identify(qq_pno) #89 (far top)
unique(sshm5$plotno)[c(89)] #plot 106 is weird (highest ranef 1.95)

dev.off()

#qqplot for cohort ranef #
png("figures.seed/modeleval-glmm.nb3-qqplot-ranef_cohort.png", width = 6, height = 5, units = "in", res=300)
qq_coh <- qqnorm(ranef(glmm.nb3, drop = TRUE)$cohort, las = 1, cex = 1.4, main="normality of ranef(glmm)$cohort", ylim=c(-6,4))
qqline(ranef(glmm.nb3, drop = TRUE)$cohort, las = 1, cex = 1.4)
text(qq_coh, labels=levels(sshm5$cohort), cex = 1, pos=3)
dev.off()

##### Interpretation ####
#seed survival increases with seed mass while controlling for max height
#no effect of max height
#seed mass is positively correlated with max.ht


#### Possible issues ####
# very high variance in cohort random effects
ranef(glmm.nb3) #2015 cohort is way different from other cohorts
summary(glmm.nb3)

#### Nicer Figures with ggplot2 ####
library(ggplot2)
library(scales)

sshm2$fit <- fitted(glmm.nb3)

## Get average seed survival value per species, and ranges ####
sshm6 <- sshm2 %>% group_by(sp, seed.mass, seed.mass.z, max.ht, max.ht.z) %>% summarize(ss.mean = mean(seedsurvival), ss.low = mean(seedsurvival)-sd(seedsurvival), ss.up = mean(seedsurvival)+sd(seedsurvival), fit.mean = mean(fit), fit.min = min(fit), fit.max = max(fit) )

write.csv(sshm6, "output.seed/lfdp_sshm6_forfig.csv", row.names=FALSE)

head(sshm6)

gl <- guide_legend(title=NULL)

## Scatterplot with predicted values ####
p1 <- ggplot(sshm6, aes(x=seed.mass, y=fit.mean)) + 
  geom_point(aes(x=seed.mass, y=fit.mean), data=sshm6, size = 3) + 
  ylab("Predicted STS ratio") + 
  xlab("Seed mass (g)")  +
  geom_linerange(aes(ymax=fit.max, ymin=fit.min)) # error bars are ranges of fitted values per species 

p11 <- p1 + scale_x_log10(breaks = c(0,0.001, 0.01,0.1,1)) + 
  scale_y_log10(breaks = c(0,10^-4, 10^-2,1,100), labels=c(NA, "0.0001", "0.01", "1", "100")) +   theme_classic() + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0,0.2,0.2),"cm"),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(linetype="solid") ) +
  guides(shape = gl, colour=gl) 

tiff(file="figures.seed/glmm.nb3_fit.vs.seedmass_avg_pred.tiff",width=5,height=3.5, units="in", res=300)
p11
dev.off()


#### Scatterplot of actual seed survival ####
#sshm6[sshm6$ss.mean == 0, ]$ss.mean <- 0.0001 #so that log10(0) doesn't turn the value into -Inf
range(sshm6$ss.low, na.rm=TRUE) # some negs
sshm6[is.na(sshm6$ss.low)==FALSE & sshm6$ss.low <= 0, ]$ss.low <- NA
sshm6[is.na(sshm6$ss.up)==FALSE & sshm6$ss.up > 1, ]$ss.up <- 1

p2 <- ggplot(sshm6, aes(x=seed.mass, y=ss.mean)) + geom_point(aes(x=seed.mass, y=ss.mean), data=sshm6, size = 5) + ylab("Seed survival (proportion)") + xlab("Seed mass (g)") 
#geom_linerange(aes(ymax=fit.max, ymin=fit.min))
#geom_linerange(aes(ymax=ss.up, ymin=ss.low))

p21 <- p2 + scale_x_log10(breaks = c(0, 0.001, 0.01,0.1,1), labels=c(NA, "0.001", "0.01", "0.1", "1")) + 
  scale_y_log10(breaks = c(0, 10^-4, 10^-2,1), labels=c(NA, "0", "0.01", "1")) + 
  theme_classic(base_size = 34) + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        axis.text=element_text(size=30),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(color="black", size=.35, linetype="solid"),
        legend.key.size=unit(1, "cm"),
        legend.margin=unit(0.01,"cm"),
        legend.text=element_text(size=26) ) +  
  guides(shape = gl, colour=gl)

tiff(file="figures.seed/glmm.nb3_seedsurvival.vs.seedmass_avg_actual.tiff",width=7.5,height=5.6, units="in", res=220)
p21
dev.off()

#### Scatterplot of actual seed survival ~ tree height w/ seed mass proportional to point size ######
p3 <- ggplot(sshm6, aes(x=max.ht, y=ss.mean, size=seed.mass.z)) + geom_point(aes(x=max.ht, y=ss.mean), data=sshm6) + ylab("STS ratio") +  xlab(expression(H[max]~(m))) #code not working currently...

p31 <- p3 + 
  scale_y_log10(breaks = c(0, 10^-4, 10^-2,1), labels=c(NA, "0", "0.01", "1")) + 
  theme_classic(base_size = 34) + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        axis.text=element_text(size=30),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(color="black", size=.35, linetype="solid"),
        legend.key.size=unit(1, "cm"),
        legend.margin=unit(0.01,"cm"),
        legend.text=element_text(size=26) ) +  
  guides(shape = gl, colour=gl)
# need to modify code to get better looking plot

tiff(file="figures.seed/glmm.nb3_seedsurvival.vs.treeht_avg_actual.tiff",width=7.5,height=5.6, units="in", res=220)
p31
dev.off()

#### Fig 1(a) Actual seed survival ~ seed mass + shaded by hmax FOR PUBLICATION ####
a3 <- ggplot(sshm6, aes(x=seed.mass, y=ss.mean, fill=max.ht)) + geom_point(aes(x=seed.mass, y=ss.mean), data=sshm6, shape=21, size=3) + ylab("STS ratio") + xlab("Seed mass (g)") 

a31 <- a3 + scale_x_log10(breaks = c(0, 0.001, 0.01,0.1,1), labels=c(NA, "0.001", "0.01", "0.1", "1")) + 
  scale_y_log10() + 
  theme_classic(base_size = 9) + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), #line size
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"), #adjust margins
        legend.position="none", #no legends
        plot.title=element_text(hjust=0, size=9) ) + #plot.title for labelling subplots
  scale_fill_gradient(low=gray(1), high =gray(0)) + 
  annotation_logticks() +
  ggtitle('(a)')

a31

a31
ggsave("lfdp_seedsurv-seedmass_avg_shadedbyheight.eps", path = "figures.seed", width = 82, height = 82, units = "mm", dpi=800) #save as eps