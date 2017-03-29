#### R_Kohyama_6.ModelTesting_SaplingRecruitment_LFDP.R ####
#### Kohyama Hypothesis Project
#### Code by Maria Wang, Summer/Fall 2016

#### Load packages ####
library(dplyr)
library(lme4)
library(afex)

# get citations for packages 
citation("lme4")
citation("afex")

#### Set working directory and load csv ####
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/")
getwd()

recruit <- read.csv("output.sapling/lfdp_recruit_lme.csv")
str(recruit)

#### Strategies for dealing with zeroes ####
#http://robjhyndman.com/hyndsight/transformations/
#http://www.mail-archive.com/r-sig-ecology@r-project.org/msg00655.html
#My folk-law guidelines on the c in log(x+c) are:
# 1. c should roughly be 1/2 of the smallest, non-zero value:  
# 2. c should be square of the first quantile devided by the third quantile (Stahel,  2002)

addtozero <- signif(0.5*sort(unique(recruit$rate))[2], 2)
addtozero2 <- (quantile(recruit$rate)[2]^2)/quantile(recruit$rate)[4]

#### Exploratory plots ####
hist(recruit$rate, breaks=30)
hist(log10(recruit$rate)) 
hist(log1p(recruit$rate)) 
hist(log(recruit$rate+addtozero))
hist(log(recruit$rate+addtozero2))

plot(rate ~ max.ht, data=recruit)
plot(log10(rate) ~ log10(max.ht), data=recruit)
plot(log1p(rate) ~ log1p(max.ht), data=recruit)

# Plot species level averages 
recruit.sp <- recruit %>% group_by(sp, max.ht) %>% summarize(rate=mean(rate),rate.sd=sd(rate))
plot((rate) ~ (max.ht), data=recruit.sp)
plot(log10(rate) ~ log10(max.ht), data=recruit.sp)
plot(log1p(rate) ~ log1p(max.ht), data=recruit.sp)


#### LINEAR MIXED EFFECTS MODELS ####
## HYPOTHESIS: Sapling recruitment increases as hmax decreases ---

lmm11 <- lmer((rate) ~ (max.ht) + (1|sp) + (1|interval), data=recruit) #slope = -
lmm0 <- lmer((rate) ~ 1 + (1|sp) + (1|interval), data=recruit)
anova(lmm11, lmm0) #p<0.0001 significant

lmm1 <- lmer(log1p(rate) ~ log1p(max.ht) + (1|sp) + (1|interval), data=recruit) #slope = -
lmm0 <- lmer(log1p(rate) ~ 1 + (1|sp) + (1|interval), data=recruit)
anova(lmm1, lmm0) #p<0.0001 significant

lmm2 <- lmer(log(rate+addtozero) ~ log(max.ht) + (1|sp) + (1|interval), data=recruit) 
lmm0 <- lmer(log(rate+addtozero) ~ 1 + (1|sp) + (1|interval), data=recruit)
anova(lmm2, lmm0) #p<0.0001 significant

lmm3 <- lmer(log(rate+addtozero2) ~ log(max.ht) + (1|sp) + (1|interval), data=recruit) 
lmm0 <- lmer(log(rate+addtozero2) ~ 1 + (1|sp) + (1|interval), data=recruit)
anova(lmm3, lmm0) #p<0.0001 significant

# Not sure if AIC is the best way to compare these models, but here goes:
AIC(lmm11, lmm1, lmm2, lmm3) #lmm1 with log1p(rate) has lowest AIC, use this model

#### Obtaining P-values using likelihood ratio test (anova) ####
lmm1 <- lmer(log1p(rate) ~ log1p(max.ht) + (1|sp) + (1|interval), data=recruit) #slope = -
lmm0 <- lmer(log1p(rate) ~ 1 + (1|sp) + (1|interval), data=recruit)
anova(lmm1, lmm0) #p<0.0001 significant
summary(lmm1)

#### Obtaining P-values using afex::mixed ####
# This P-value uses Kenward-Roger approximation for degrees of freedom . See ?mixed
mixed(log1p(rate) ~ log1p(max.ht) + (1|sp) + (1|interval), data=recruit, method = "KR")
#           Effect       df F.scaling         F p.value
# 1 log1p(max.ht) 1, 90.02      1.00 40.56 ***  <.0001




###############################################
#### Model Evaluation and Diagnostic Plots ####

# Plot random effects against the predicted values from the fixed effect component of the model and check for no trend: 
recruit$ft.fix   <- as.numeric(model.matrix(lmm1) %*% fixef(lmm1))

# random effects
rr.sp          <- ranef(lmm1, drop = TRUE)$sp 
rr.interval      <- ranef(lmm1, drop = TRUE)$interval 

# create df with random effect levels and values
rrr.sp    <- data.frame(sp = attr(rr.sp, "names"), ranef.sp = rr.sp)
rrr.interval<- data.frame(interval = attr(rr.interval, "names"), ranef.interval = rr.interval) 
recruit$interval <- as.factor(recruit$interval)

recruit5 <- left_join(recruit, rrr.sp) %>% left_join(rrr.interval) %>% arrange(sp, interval)

# The main diagnostic plot should be estimated random effects vs. fitted values from a GLMM fit
png("figures.sapling/lfdp-modeleval-lmm1-ranef.vs.fitted.png", width = 24, height = 12, units = "cm", res=150)
par(mfcol=c(1,2))
with(recruit5, plot(ft.fix, ranef.sp, pch = 19, las = 1, cex = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
with(recruit5, plot(ft.fix, ranef.interval, pch = 19, las = 1, cex = 1, main="ranef vs fitted - interval"))
abline(0,0,lwd = 1.5)
dev.off()
#one point in the spp plot that seems like an outlier

## use identify to find out which is outlier ##
with(recruit5, plot(ft.fix, ranef.sp, pch = 19, las = 1, cex = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
identify(recruit5$ft.fix, recruit5$ranef.sp, labels=recruit5$sp)
recruit5[133, ]$sp #GONSPI
recruit5[recruit5$sp == "GONSPI", ]
#GONSPI has a rather high random effect

write.csv(recruit5, "output.sapling/lfdp_recruit5.csv")

## now check for approximate normality of random effects: ####
png("figures.sapling/lfdp-modeleval-lmm1-qqplot-ranef.png", width = 24, height = 12, units = "cm", res=150)
par(mfcol=c(1,2))
qqnorm(ranef(lmm1, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(lmm)$sp")
qqline(ranef(lmm1, drop = TRUE)$sp, las = 1, cex = 1.4) #

qqnorm(ranef(lmm1, drop = TRUE)$interval, las = 1, cex = 1.4, main="normality of ranef(lmm)$interval")
qqline(ranef(lmm1, drop = TRUE)$interval, las = 1, cex = 1.4) #seems fine 
dev.off()
#hmm only 1 spp that's much higher than everyone else

# use identify to find out which is the outlier
qq_sp <- qqnorm(ranef(lmm1, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(lmm)$sp")
qqline(ranef(lmm1, drop = TRUE)$sp, las = 1, cex = 1.4) #
identify(qq_sp)
ranef(lmm1)$sp[36, ]
unique(recruit$sp)[c(36) ] #GONSPI has a rather high random effect

ranef(lmm1)$interval 


##### VISUALIZATION: MAKING FIGURES ####
#use ggplot2 scale transformation for log x axis without transforming numbers
library(ggplot2)
library(scales)

# get model fitted/predicted values 
recruit$fit <- fitted(lmm1)

# get avg values per species
recruit6 <- recruit %>% group_by(sp, max.ht) %>% summarize(rr.mean = mean(rate), rr.min = min(rate), rr.max=max(rate), fit.mean = mean(fit), fit.min = min(fit), fit.max=max(fit))
write.csv(recruit6, "output.sapling/lfdp_recruit6_forfig.csv")

## scatterplot of fitted (predicted) values ####
e2 <- ggplot(recruit6, aes(x=max.ht, y=rr.mean)) + geom_point(aes(x=max.ht, y=rr.mean), data=recruit6, shape=1) + ylab(expression(Sapling~recruitment~(no.~y^-1))) + 
  xlab(expression(H[max]~(m))) +
  geom_linerange(aes(ymax=fit.max, ymin=fit.min))

e32 <- e2 + theme_classic(base_size = 9) + scale_y_log10(breaks=c(0,0.05,0.1,0.5,1)) + scale_x_continuous(breaks=c(0,10,20,30)) +
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))

png(file="figures.sapling/lfdp-lmm1_recruit.vs.maxht_avg_pred.png", width = 82, height = 82, units = "mm", res=1200)
e32
dev.off()

## scatterplot of actual values ####
e3 <- ggplot(recruit6, aes(x=max.ht, y=rr.mean)) + geom_point(aes(x=max.ht, y=rr.mean), data=recruit6, shape=1) + ylab(expression(Sapling~recruitment~(no.~y^-1))) + 
  xlab(expression(H[max]~(m))) +
  geom_smooth(aes(x=max.ht, y=fit.mean), color="black", size=0.75)

e33 <- e3 + theme_classic(base_size = 9) + scale_y_log10(breaks=c(0,0.05,0.1,0.5,1)) + scale_x_continuous(breaks=c(0,10,20,30)) +
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        plot.title=element_text(hjust=0, size=9) ) + #plot.title for labelling subplots
  ggtitle('(c)')

e33
png(file="figures.sapling/lfdp-lmm1_recruit.vs.maxht_avg_actual_trendline.png", width = 82, height = 82, units = "mm", res=1200)
e33
dev.off()

# e33
# ggsave("lfdp-lmm1_recruit.vs.maxht_avg_actual.eps", path = "figures.sapling", width = 82, height = 82, units = "mm", dpi=1200) #save as eps


