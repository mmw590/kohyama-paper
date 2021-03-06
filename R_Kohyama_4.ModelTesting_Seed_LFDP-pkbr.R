#### R_Kohyama_4.ModelTesting_Seed_LFDP.R ####
#### Kohyama Hypothesis Project
#### Code by Maria Wang, Summer/Fall 2016
install.packages("dplyr", "lme4", "afex")

# Load packages
library(dplyr)
library(lme4)
library(afex)
citation("glmmADMB")

#### Set working directory and load csv ####
setwd("C:/Users/user/Downloads/Kohyama Hypothesis")
getwd()

sshm2 <- read.csv("output.seed/seedsurvival_forglm.csv")
str(sshm2)

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

#### Hypothesis: Seed to Seedling ratio ####
# Alt: STS ratio is affected by max tree height and seed mass
# Null: STS ratio is not affected by max tree height and seed mass

#
#### Testing models with height as continuous var (not centered) ####
glmm01 <- glmer(seedsv ~ max.ht.z * seed.mass.z + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial) #might take a minute
glmm02 <- glmer(seedsv ~ max.ht.z + seed.mass.z + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)

# use likelihood ratio test to compare models
anova(glmm01, glmm02) #interaction not significant, p=0.53

summary(glmm01) #interaction not significant, p=0.59
summary(glmm02) #seed mass is sig (p=0.004), but not max.ht.z (0.64)

glmm03 <- glmer(seedsv ~ max.ht + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)
summary(glmm03) #max.ht not sig, p=0.089

library(pbkrtest)
PBmodcomp(glmm01, glmm02, nsim=1000) #40min with nsim=10, #3pm-10pm, not done after 17 hours
# Parametric bootstrap test; time: 119264.74 sec; samples: 1000 extremes: 606;
# Requested samples: 1000 Used samples: 980 Extremes: 606
# large : seedsv ~ max.ht.z * seed.mass.z + (1 | sp) + (1 | cohort) + (1 | 
#                                                                        plotno)
# small : seedsv ~ max.ht.z + seed.mass.z + (1 | sp) + (1 | cohort) + (1 | 
#                                                                        plotno)
# stat df p.value
# LRT    0.3074  1  0.5793
# PBtest 0.3074     0.6188

# Parametric bootstrap test; time: 1030.87 sec; samples: 10 extremes: 5;
# large : seedsv ~ max.ht.z * seed.mass.z + (1 | sp) + (1 | cohort) + (1 | 
#                                                                        plotno)
# small : seedsv ~ max.ht.z + seed.mass.z + (1 | sp) + (1 | cohort) + (1 | 
#                                                                        plotno)
# stat df p.value
# LRT    0.3074  1  0.5793
# PBtest 0.3074     0.5455

#### use package afex ####
mixed(seedsurvival ~ max.ht.z * seed.mass.z + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, family=binomial, method = "PB", weights=sshm2$totseeds.x3, args.test=list(nsim=1000)) #7.04pm-



#

# Model evaluation and Diagnostic plots ####

# Plot random effects against the predicted values from the fixed effect component of the model and check for no trend: 
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

# The main diagnostic plot should be estimated random effects vs. fitted values from a GLMM fit.
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


## now check for approximate normality of random effects: ####
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

##### Overdispersion test ####
# if much larger than one, it may be overdispersed, esp if there are not many observations...
deviance(glmm14)/df.residual(glmm14) # = 2.15... is it overdispersed?

##### Interpretation ####
#seed survival increases with seed mass
#no difference in seed survival between understory and canopy trees
#seed mass is positively correlated with max.ht

#Questions: should remove outliers observed in qqplots?

#### plot best model ####

## make a plot of predicted seed survival (fit) against seedmass
#glmm14 <- glmer(seedsv ~ seed.mass + (1|sp) + (1|cohort) + (1|plotno), data=sshm2, binomial)

sshm2$fit <- fitted(glmm14) 

# quick plot using base R package
png("figures.seed/glmm14-fit.vs.seedmass.png", width=18, height=18, units="cm", res=300)
plot(log10(fit) ~ log10(seed.mass), data = sshm2, type="n", ylab="log10(predicted seed survival)", xlab="log10(seed mass (g))")
points(log10(fit) ~ log10(seed.mass), data = sshm2[sshm2$htclass == "understory", ], col="brown", pch=1)
points(log10(fit) ~ log10(seed.mass), data = sshm2[sshm2$htclass == "canopy", ], col="skyblue", pch=2)
legend("bottomright", legend=c("understory", "canopy"), col=c("brown", "skyblue"), pch=c(1,2))
dev.off()

# use ggplot2 to make nice looking plot  ####
library(ggplot2)
library(scales)
#library(grid)
#library(gridExtra)

# Change order of levels so that understory comes first 
sshm2$htclass <- factor(sshm2$htclass,c("understory","canopy"))

# Get average seed survival value per species ####
sshm6 <- sshm2 %>% group_by(sp, seed.mass, seed.mass.z, max.ht, max.ht.z, htclass) %>% summarize(ss.mean = mean(seedsurvival), ss.low = mean(seedsurvival)-sd(seedsurvival), ss.up = mean(seedsurvival)+sd(seedsurvival), fit.mean = mean(fit), fit.min = min(fit), fit.max = max(fit) )

write.csv(sshm6, "output.seed/lfdp_sshm6_forfig.csv", row.names=FALSE)

head(sshm6)

gl <- guide_legend(title=NULL)

#### Scatterplot with predicted values ####
p1 <- ggplot(sshm6, aes(x=seed.mass, y=fit.mean, shape=htclass, colour=htclass)) + 
  geom_point(aes(x=seed.mass, y=fit.mean), data=sshm6, size = 3) + 
  ylab("Predicted seed survival") + 
  xlab("Seed mass (g)")  +
  geom_linerange(aes(ymax=fit.max, ymin=fit.min)) # error bars are ranges of fitted values per species 

# scale allows you to plot log-transformed values but label axes with untransformed values
p11 <- p1 + scale_x_log10(breaks = c(0,0.001, 0.01,0.1,1)) + 
  scale_y_log10(breaks = c(0,10^-6, 10^-4, 10^-2,1), labels=c(NA, "0.000001", "0.000100", "0.010000", "1.000000")) + 
  theme_classic() + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0,0.2,0.2),"cm"),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(linetype="solid") ) +
  guides(shape = gl, colour=gl) 

tiff(file="figures.seed/glmm14_fit.vs.seedmass_avg_pred.tiff",width=5,height=3.5, units="in", res=300)
p11
dev.off()


#### Scatterplot of actual seed survival ####
sshm6[sshm6$ss.mean > 1, ]$ss.mean <- 1
sshm6[sshm6$ss.mean == 0, ]$ss.mean <- 0.0001 #so that log10(0) doesn't turn the value into -Inf
range(sshm6$ss.low, na.rm=TRUE) # some negs
sshm6[is.na(sshm6$ss.low)==FALSE & sshm6$ss.low <= 0, ]$ss.low <- NA
sshm6[is.na(sshm6$ss.up)==FALSE & sshm6$ss.up > 1, ]$ss.up <- 1

p2 <- ggplot(sshm6, aes(x=seed.mass, y=ss.mean, shape=htclass, colour=htclass)) + geom_point(aes(x=seed.mass, y=ss.mean), data=sshm6, size = 5) + ylab("Seed survival (proportion)") + xlab("Seed mass (g)") 
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

tiff(file="figures.seed/glmm14_seedsurvival.vs.seedmass_avg_actual.tiff",width=7.5,height=5.6, units="in", res=220)
p21
dev.off()

#### Scatterplot of actual seed survival ~ tree height w/ seed mass proportional to point size ######
p3 <- ggplot(sshm6, aes(x=max.ht, y=ss.mean, size=seed.mass.z)) + geom_point(aes(x=max.ht, y=ss.mean), data=sshm6) + ylab("Seed survival (proportion)") + xlab("Max tree height (m)") #code not working currently...

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

tiff(file="figures.seed/glmm14_seedsurvival.vs.treeht_avg_actual.tiff",width=7.5,height=5.6, units="in", res=220)
p31
dev.off()

#### Actual seed survival ~ seed mass + shaded by hmax FOR PUBLICATION ####

#size of points = 3mm wide
a3 <- ggplot(sshm6, aes(x=seed.mass, y=ss.mean, fill=max.ht)) + geom_point(aes(x=seed.mass, y=ss.mean), data=sshm6, shape=21, size=3) + ylab("Seed survival (proportion)") + xlab("Seed mass (g)") 

a31 <- a3 + scale_x_log10(breaks = c(0, 0.001, 0.01,0.1,1), labels=c(NA, "0.001", "0.01", "0.1", "1")) + 
  scale_y_log10(breaks = c(0, 10^-4, 10^-2,1), labels=c(NA, "0", "0.01", "1")) + 
  theme_classic(base_size = 9) + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        legend.position="none") +
  scale_fill_gradient(low=gray(1), high =gray(0)) + annotation_logticks()

png(file="figures.seed/lfdp_seedsurv-seedmass_avg_shadedbyheight.png",width=82,height=82, units="mm", res=300)
a31
dev.off()

a31
ggsave("lfdp_seedsurv-seedmass_avg_shadedbyheight.eps", path = "figures.seed", width = 82, height = 82, units = "mm", dpi=800) #save as eps