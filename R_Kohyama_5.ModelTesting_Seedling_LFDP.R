#' ---
#' title: "R_Kohyama_5.ModelTesting_Seedling_LFDP.R"
#' author: "Maria Wang"
#' date: "2017-03-28"
#' description: "This script tests for effects of seed mass and species tree height on first year seedling survival"
#' ---


# ~ Loading required R packages ----
library(dplyr)
library(lme4)
#library(afex)

# ~ Set working directory ----
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis")
getwd()

# ~ Create folders to store output ----
dir.create("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/figures.seedling/")

# ~ Load required input files----
slinghm <- read.csv("intermediate data/lfdp_seedlinghm_forglm.csv")
str(slinghm)


#### Correlation Between Seed Mass And Tree Height ####
cor(slinghm$max.ht, slinghm$seed.mass) #0.787 high corr in this dataset


#### Hypothesis: Seedling survival ####
# Alt: Seedling survival is affected by max tree height and seed mass
# Null: Seedling survival is not affected by max tree height and seed mass

# Read this on p-values in mixed models ####
help("pvalues",package="lme4")

#### Generalized Linear Mixed Model with Binomial Errors ####
#### ~ Testing interaction between seed mass and height ####
glmm81 <- glmer(live.in.y1 ~ max.ht.z * seed.mass.z + (1|sp) + (1|cohort) + (1|plotno/subplot), data=slinghm, binomial) #might take a minute
glmm82 <- glmer(live.in.y1 ~ max.ht.z + seed.mass.z + (1|sp) + (1|cohort) + (1|plotno/subplot), data=slinghm, binomial)

# use likelihood ratio test to compare models
anova(glmm81, glmm82) #interaction not significant, p=0.2641

# use Parametric Bootstrap to compare models (takes a long time)
#PBmodcomp(glmm81, glmm82, nsim=1000)

# p-values from summary() are based on Wald's chi-square tests, which are worse than likelihood ratio test, so I'll use the lrt P-values for now. (or I can use parametric bootstrap if they finish running)

glmm92 <- glmer(live.in.y1 ~ max.ht + seed.mass + (1|sp) + (1|cohort) + (1|plotno/subplot), data=slinghm, binomial)
glmm93 <- glmer(live.in.y1 ~ seed.mass + (1|sp) + (1|cohort) + (1|plotno/subplot), data=slinghm, binomial)
glmm94 <- glmer(live.in.y1 ~ max.ht + (1|sp) + (1|cohort) + (1|plotno/subplot), data=slinghm, binomial)
anova(glmm92, glmm93) #max.ht not important after controlling for seedmass (p=0.3495)
anova(glmm92, glmm94) #seed mass not important after controlling for hmax (p=0.2598)

#Based on AIC, glmm93 is slightly better than others

#mixed(live.in.y1 ~ max.ht * seed.mass + (1|sp) + (1|cohort) + (1|plotno/subplot), data=slinghm, family=binomial, method="PB", args.test=list(nsim=1000)) #takes a long time

#### Interpretation of Results ####
# None of the models were significant.
# The "best" model (lowest AIC score) is one where the fixed effect is seed mass alone.
# Does this mean that the random effects were more important than the fixed effects?

#### ***Assuming one of the models was significant, you can run the following diagnostic tests*** ####

#### Model evaluation and Diagnostic plots ####

# ~ Plot random effects against the predicted values from the fixed effect component of the model and check for no trend ----
slinghm$ft.fix   <- as.numeric(model.matrix(glmm92) %*% fixef(glmm92))

# random effects
rr.sp          <- ranef(glmm92, drop = TRUE)$sp 
rr.cohort      <- ranef(glmm92, drop = TRUE)$cohort 
rr.plotno      <- ranef(glmm92, drop = TRUE)$plotno #

# create df with random effect levels and values
rrr.sp    <- data.frame(sp = attr(rr.sp, "names"), ranef.sp = rr.sp)
rrr.cohort<- data.frame(cohort = attr(rr.cohort, "names"), ranef.cohort = rr.cohort) 
rrr.plotno<- data.frame(plotno = attr(rr.plotno, "names"), ranef.plotno = rr.plotno) 
slinghm$cohort <- as.factor(slinghm$cohort)
slinghm$plotno <- as.factor(slinghm$plotno)

sshm15 <- left_join(slinghm, rrr.sp) %>% left_join(rrr.cohort) %>% left_join(rrr.plotno) %>% arrange(sp, cohort, plotno)

# The main diagnostic plot should be estimated random effects vs. fitted values from a GLMM fit
png("figures.seedling/modeleval-glmm92-ranef.vs.fitted.png", width = 27, height = 10, units = "cm", res=150)
par(mfcol=c(1,3))
with(sshm15, plot(ft.fix, ranef.sp, pch = 19, las = 1, cex = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
with(sshm15, plot(ft.fix, ranef.cohort, pch = 19, las = 1, cex = 1, main="ranef vs fitted - cohort"))
abline(0,0,lwd = 1.5)
with(sshm15, plot(ft.fix, ranef.plotno, pch = 19, las = 1, cex = 1, main="ranef vs fitted - plotno"))
abline(0,0,lwd = 1.5)
dev.off()
#this plot actually looks pretty good

#png("figures.seedling/modeleval-glmm83-ranef.vs.fitted-sp.png", width = 18, height = 12, units = "cm", res=150)
with(sshm15, plot(ft.fix, ranef.sp, pch = 19, las = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
text(sshm15$ft.fix, sshm15$ranef.sp, labels=sshm15$sp, pos=4, cex = 0.5)
dev.off()

#png("figures.seedling/modeleval-glmm83-ranef.vs.fitted-cohort.png", width = 18, height = 12, units = "cm", res=150)
with(sshm15, plot(ft.fix, ranef.cohort, pch = 19, las = 1, main="ranef vs fitted - cohort"))
abline(0,0,lwd = 1.5)
#text(sshm15$ft.fix, sshm15$ranef.cohort, labels=unique(sshm15$sp), pos=4, cex=0.5)
identify(sshm15$ft.fix, sshm15$ranef.cohort, labels=sshm15$cohort)
#cohort 2009 is a bit out there?
dev.off()

#png("figures.seedling/modeleval-glmm83-ranef.vs.fitted-plotno.png", width = 25, height = 12, units = "cm", res=150)
with(sshm15, plot(ft.fix, ranef.plotno, pch = 19, las = 1, main="ranef vs fitted - plotno"))
abline(0,0,lwd = 1.5)
identify(sshm15$ft.fix, sshm15$ranef.plotno, labels=sshm15$sp)
#text(sshm15$ft.fix, sshm15$ranef.plotno, labels=sshm15$plotno, cex=0.5, pos=4)
#dev.off()
#INGLAU is a bit high

x11()
with(sshm15, plot(ft.fix, ranef.plotno, pch = 19, las = 1, main="ranef vs fitted - plotno"))
abline(0,0,lwd = 1.5)
identify(sshm15$ft.fix, sshm15$ranef.plotno, labels=sshm15$sp)
dev.off()
#INGLAU is a bit high

write.csv(sshm15, "output.seedling/sshm15_sling_glmm82.csv", row.names=FALSE)
#seems to be some trend rather than no trend...


#### ~ Check for approximate normality of random effects: ####
png("figures.seedling/modeleval-glmm83-qqplot-ranef.png", width = 24, height = 8, units = "cm", res=150)
par(mfcol=c(1,3))
qqnorm(ranef(glmm83, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(glmm)$sp")
qqline(ranef(glmm83, drop = TRUE)$sp, las = 1, cex = 1.4) #steps?!

qqnorm(ranef(glmm83, drop = TRUE)$cohort, las = 1, cex = 1.4, main="normality of ranef(glmm)$cohort")
qqline(ranef(glmm83, drop = TRUE)$cohort, las = 1, cex = 1.4) #seems fine except for the last year?

qqnorm(ranef(glmm83, drop = TRUE)$plotno, las = 1, cex = 1.4, main="normality of ranef(glmm)$plotno")
qqline(ranef(glmm83, drop = TRUE)$plotno, las = 1, cex = 1.4) 
dev.off()

### # use identify to label outliers - only cohort has an outlier !
# par(mfcol=c(1,1))
# qq_sp <- qqnorm(ranef(glmm83, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(glmm)$sp")
# qqline(ranef(glmm83, drop = TRUE)$sp, las = 1, cex = 1.4)
# identify(qq_sp)
# unique(sshm15$sp)[c(2,14,15,8) ] 

qq_coh <- qqnorm(ranef(glmm83, drop = TRUE)$cohort, las = 1, cex = 1.4, main="normality of ranef(glmm)$cohort")
qqline(ranef(glmm83, drop = TRUE)$cohort, las = 1, cex = 1.4)
identify(qq_coh)
levels(sshm15$cohort)[6] #2013 is weird year
ranef(glmm83)$cohort

# qq_pno <- qqnorm(ranef(glmm83, drop = TRUE)$plotno, las = 1, cex = 1.4, main="normality of ranef(glmm)$plotno")
# qqline(ranef(glmm83, drop = TRUE)$plotno, las = 1, cex = 1.4)
# identify(qq_pno) #8, 32 (far bottom), 2 (bottom), 109 (top)
# [c(8,32,2,109)] #these plots are weird
# range(ranef(glmm83)$plotno)
# dev.off()

#### ~ Overdispersion test ####
# if much larger than one, it is maybe overdispersed, esp if there are not many observations...
deviance(glmm92)/df.residual(glmm92) #is < 1

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
overdisp_fun(glmm92) #this model is not overdispersed.


#### ~ Plot Predicted Values/Fit ####

## next, make a plot of predicted seedling survival (fit) against seedmass

slinghm$fit <- fitted(glmm92) 
slinghm$htclass <- factor(slinghm$htclass,c("understory","canopy"))

#### Possible issues ####
# very high variance in cohort random effects
ranef(glmm92) #2015 cohort is way different from other cohorts
summary(glmm92)

#### Visualizing Results ####

## ~ get average seedling survival value per species ####
str(slinghm)
slinghm.sp <- slinghm %>% group_by(sp, seed.mass, max.ht, cohort, plotno) %>% summarize(live.in.y1=sum(live.in.y1), nsling = n(), slingsurvival.y1=live.in.y1/nsling) %>% group_by(sp, seed.mass, max.ht) %>% summarize(sl.mean = mean(slingsurvival.y1))

range(slinghm.sp$sl.mean) #no zeroes

#### Fig1(b) Actual seedling survival ~ seed mass + shaded by hmax FOR PUBLICATION ####
names(slinghm.sp)
#size of points = 3mm wide

b3 <- ggplot(slinghm.sp, aes(x=seed.mass, y=sl.mean, fill=max.ht)) + geom_point(aes(x=seed.mass, y=sl.mean), data=slinghm.sp, shape=21, size=3) + ylab("Seedling survival (proportion)") + xlab("Seed mass (g)") 

b31 <- b3 + scale_x_log10(breaks = c(0, 0.001, 0.01,0.1,1), labels=c(NA, "0.001", "0.01", "0.1", "1")) + 
  scale_y_log10(breaks = c(0, 0.01, 0.1, 1), labels=c("0", "0.01", "0.1", "1")) + 
  theme_classic(base_size = 9) + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        legend.position="none",
        plot.title=element_text(hjust=0, size=9) ) +
  scale_fill_gradient(low=gray(1), high =gray(0)) + 
  annotation_logticks() +
  ggtitle('(b)')

b31

png(file="figures.seedling/lfdp_slingsurv-seedmass_avg_shadedbyheight.png",width=82,height=82, units="mm", res=300)
b31
dev.off()

# b31
# ggsave("lfdp_slingsurv-seedmass_avg_shadedbyheight.eps", path = "figures.seedling", width = 82, height = 82, units = "mm", dpi=800) #save as eps

#
#
#
#### ****the rest of this code has not been updated. exploratory code**** ####

# ~ quick plot using base R package ----

# png("figures.seedling/glmm82-fit.vs.seedmass.png", width=18, height=18, units="cm", res=300)
# plot(log10(fit) ~ log10(seed.mass), data = slinghm)
# points(log10(fit) ~ log10(seed.mass), data = slinghm[slinghm$htclass == "understory", ], col="brown")
# points(log10(fit) ~ log10(seed.mass), data = slinghm[slinghm$htclass == "canopy", ], col="skyblue")
# legend("bottomright", legend=c("understory", "canopy"), col=c("brown", "skyblue"), pch=1)
# #abline(fixef(glmm92), untf=FALSE)
# #abline(log10(fixef(glmm92)[1]), log10(fixef(glmm92)[2]) )
# dev.off()

png("figures.seedling/glmm92-fit.vs.max.ht.png", width=18, height=18, units="cm", res=300)
plot(log10(fit) ~ max.ht, data = slinghm, type="n", ylab="log10(predicted seedling survival)", xlab="max potential tree height (m)")
points(log10(fit) ~ max.ht, data = slinghm[slinghm$htclass == "understory", ], col="brown", pch=1)
points(log10(fit) ~ max.ht, data = slinghm[slinghm$htclass == "canopy", ], col="skyblue", pch=2)
legend("bottomright", legend=c("understory", "canopy"), col=c("brown", "skyblue"), pch=c(1,2))
dev.off()


png("figures.seedling/glmm92-fit.vs.max.ht_seedmass.png", width=18, height=18, units="cm", res=300)
plot(log10(fit) ~ log10(seed.mass), data = slinghm, type="n", ylab="log10(predicted seedling survival)", xlab="log10(seed mass (g))")
points(log10(fit) ~ log10(seed.mass), data = slinghm[slinghm$htclass == "understory", ], col="brown", pch=1)
points(log10(fit) ~ log10(seed.mass), data = slinghm[slinghm$htclass == "canopy", ], col="skyblue", pch=2)
legend("bottomright", legend=c("understory", "canopy"), col=c("brown", "skyblue"), pch=c(1,2))
dev.off()


#### ~ Scatterplots with ggplot2 ####
library(ggplot2)
library(scales)

names(slinghm)
slinghm$fit <- fitted(glmm92)

slinghm.sp <- slinghm %>% group_by(sp, seed.mass, max.ht) %>% summarize(live.in.y1=sum(live.in.y1), nsling = n(), slingsurvival.y1=live.in.y1/nsling) %>% mutate(sl.mean = mean(slingsurvival.y1), sl.min = min(slingsurvival.y1), sl.max = max(slingsurvival.y1), fit.mean = mean(fit), fit.min = min(fit), fit.max = max(fit) )

gl <- guide_legend(title=NULL)

## ~ Scatterplot with predicted values ####
p9 <- ggplot(slinghm.sp, aes(x=max.ht, y=fit.mean)) + geom_point(aes(x=max.ht, y=fit.mean), data=slinghm.sp, size = 3) + ylab("Predicted seedling survival") +  xlab(expression(H[max]~(m))) +
  geom_linerange(aes(ymax=fit.max, ymin=fit.min))

p19 <- p9 + theme_classic() + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0,0.2,0.2),"cm"),
        legend.position = c(1,1),
        legend.justification = c(1, 1),
        legend.background = element_rect(color="black", size=.25, linetype="solid") ) +
  guides(shape = gl, colour=gl) 

tiff(file="figures.seedling/glmm93_sling_fit.vs.maxht_avg_pred.tiff",width=5,height=3.5, units="in", res=300)
p19
dev.off()


#### ~ Scatterplot of actual sling survival ####
p8 <- ggplot(slinghm.sp, aes(x=max.ht, y=sl.mean, shape=htclass, colour=htclass)) + geom_point(aes(x=max.ht, y=sl.mean), data=slinghm.sp, size = 5) + ylab("Seedling survival (proportion)") +  xlab(expression(H[max]~(m)))

p18 <- p8 + theme_classic(base_size = 34) + 
  scale_x_continuous(breaks=c(0,10,20,30)) +
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        axis.text=element_text(size=30),
        legend.position = c(1,1),
        legend.justification = c(1, 1),
        legend.background = element_rect(color="black", size=.25, linetype="solid"),
        legend.key.size=unit(1, "cm"),
        legend.margin=unit(0.01,"cm"),
        legend.text=element_text(size=26) ) +
  guides(shape = gl, colour=gl) 

tiff(file="figures.seedling/glmm93_slingsurv.vs.maxht_avg_actual.tiff",width=7.5,height=5.6, units="in", res=220)
p18
dev.off()


