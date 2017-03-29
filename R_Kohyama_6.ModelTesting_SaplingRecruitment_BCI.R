#### R_Kohyama_6.ModelTesting_SaplingRecruitment_BCI.R ####
#### Kohyama Hypothesis Project
#### Code by Maria Wang, Summer/Fall 2016

#### Load packages ####
library(dplyr)
library(lme4)
library(afex)

#### Set working directory and load csv ####
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis/")
getwd()

brecruit <- read.csv("output.sapling/bci_recruit_lme.csv")

#### Strategies for dealing with zeroes ####
#http://robjhyndman.com/hyndsight/transformations/
#http://www.mail-archive.com/r-sig-ecology@r-project.org/msg00655.html
#My folk-law guidelines on the c in log(x+c) are:
# 1. c should roughly be 1/2 of the smallest, non-zero value:  
# 2. c should be square of the first quantile devided by the third quantile (Stahel,  2002)

addtozero <- signif(0.5*sort(unique(brecruit$rate))[2], 2)
addtozero2 <- (quantile(brecruit$rate)[2]^2)/quantile(brecruit$rate)[4]


#### Exploratory plots ####
par(mfcol=c(2,2))
hist(brecruit$rate)
hist(log10(brecruit$rate)) #looks the most normal
hist(log(brecruit$rate)) #looks quite normal
hist(log1p(brecruit$rate)) #not normal
hist(log(brecruit$rate+addtozero))
hist(log(brecruit$rate+addtozero2))
par(mfcol=c(1,1))

plot(rate ~ max.ht, data=brecruit)
plot(log10(rate) ~ log10(max.ht), data=brecruit)
plot(log1p(rate) ~ log1p(max.ht), data=brecruit)

# Plot species level averages 
brecruit.sp <- brecruit %>% group_by(sp, max.ht) %>% summarize(rate=mean(rate))
plot((rate) ~ (max.ht), data=brecruit.sp)
plot(log10(rate) ~ log10(max.ht), data=brecruit.sp)
plot(log(rate) ~ log(max.ht), data=brecruit.sp)
plot(log1p(rate) ~ log1p(max.ht), data=brecruit.sp)


#### LINEAR MIXED EFFECTS MODELS ####
library(lme4)
str(recruit2$rate)

blmm11 <- lmer((rate) ~ (max.ht) + (1|sp) + (1|interval), data=brecruit) #slope = -
blmm0 <- lmer((rate) ~ 1 + (1|sp) + (1|interval), data=brecruit)
anova(blmm11, blmm0) #p<0.0001 significant
summary(blmm11) 

blmm1 <- lmer(log1p(rate) ~ log1p(max.ht) + (1|sp) + (1|interval), data=brecruit) #slope = -
blmm0 <- lmer(log1p(rate) ~ 1 + (1|sp) + (1|interval), data=brecruit)
anova(blmm1, blmm0) #p<0.0001 significant
summary(blmm1) 

blmm2 <- lmer(log(rate+addtozero) ~ log(max.ht) + (1|sp) + (1|interval), data=brecruit) 
blmm0 <- lmer(log(rate+addtozero) ~ 1 + (1|sp) + (1|interval), data=brecruit)
anova(blmm2, blmm0) #p<0.0001 significant
summary(blmm2)

blmm3 <- lmer(log(rate+addtozero2) ~ log(max.ht) + (1|sp) + (1|interval), data=brecruit) 
blmm0 <- lmer(log(rate+addtozero2) ~ 1 + (1|sp) + (1|interval), data=brecruit)
anova(blmm3, blmm0) #p<0.0001 significant
summary(blmm3)

# Not sure if AIC is the best way to compare these models, but here goes:
AIC(blmm11, blmm1, blmm2, blmm3) #blmm1 with log1p(rate) has lowest AIC, use this model

#### Obtaining P-values using likelihood ratio test (anova) ####
lmm1 <- lmer(log1p(rate) ~ log1p(max.ht) + (1|sp) + (1|interval), data=recruit) #slope = -
lmm0 <- lmer(log1p(rate) ~ 1 + (1|sp) + (1|interval), data=recruit)
anova(lmm1, lmm0) #p<0.0001 significant
summary(lmm1)

#### Obtaining P-values using afex::mixed ####
# This P-value uses Kenward-Roger approximation for degrees of freedom . See ?mixed
mixed(log1p(rate) ~ log1p(max.ht) + (1|sp) + (1|interval), data=brecruit, method = "KR")
# Effect        df F.scaling         F p.value
# 1 log1p(max.ht) 1, 250.04      1.00 75.59 ***  <.0001


###############################################
#### Model Evaluation and Diagnostic Plots ####

# Plot random effects against the predicted values from the fixed effect component of the model and check for no trend:
brecruit$ft.fix   <- as.numeric(model.matrix(blmm1) %*% fixef(blmm1))

# random effects
rr.sp          <- ranef(blmm1, drop = TRUE)$sp 
rr.interval      <- ranef(blmm1, drop = TRUE)$interval 

# create df with random effect levels and values
rrr.sp    <- data.frame(sp = attr(rr.sp, "names"), ranef.sp = rr.sp)
rrr.interval<- data.frame(interval = attr(rr.interval, "names"), ranef.interval = rr.interval) 
brecruit$interval <- as.factor(brecruit$interval)
brecruit5 <- left_join(brecruit, rrr.sp) %>% left_join(rrr.interval) %>% arrange(sp, interval)

# The main diagnostic plot should be estimated random effects vs. fitted values from a Gblmm fit
png("figures.sapling/bci-modeleval-blmm1-ranef.vs.fitted.png", width = 24, height = 12, units = "cm", res=150)
par(mfcol=c(1,2))
with(brecruit5, plot(ft.fix, ranef.sp, pch = 19, las = 1, cex = 1, main="ranef vs fitted - species"))
abline(0,0,lwd = 1.5)
with(brecruit5, plot(ft.fix, ranef.interval, pch = 19, las = 1, cex = 1, main="ranef vs fitted - interval"))
abline(0,0,lwd = 1.5)
dev.off()
# looks like no trend, yay!

write.csv(brecruit5, "output.sapling/bci_recruit5.csv", row.names=FALSE)
brecruit5 <- read.csv("output.sapling/bci_recruit5.csv")

## now check for approximate normality of random effects: ####
png("figures.sapling/bci-modeleval-blmm1-qqplot-ranef.png", width = 24, height = 12, units = "cm", res=150)
par(mfcol=c(1,2))
qqnorm(ranef(blmm1, drop = TRUE)$sp, las = 1, cex = 1.4, main="normality of ranef(blmm)$sp")
qqline(ranef(blmm1, drop = TRUE)$sp, las = 1, cex = 1.4) #

qqnorm(ranef(blmm1, drop = TRUE)$interval, las = 1, cex = 1.4, main="normality of ranef(blmm)$interval")
qqline(ranef(blmm1, drop = TRUE)$interval, las = 1, cex = 1.4) #seems fine 
dev.off()
#hmm there's a banana shape for the qqplot of spp random effects

ranef(blmm1)$interval 


##### VISUALIZATION: MAKING FIGURES ####
#use ggplot2 scale transformation for log x axis without transforming numbers
library(ggplot2)
library(scales)

# get model fitted/predicted values 
brecruit$fit <- fitted(blmm1)


#get avg values per spp  
brecruit6 <- brecruit %>% group_by(sp, max.ht) %>% summarize(rr.mean = mean(rate), rr.min = min(rate), rr.max=max(rate), fit.mean = mean(fit), fit.min = min(fit), fit.max=max(fit))
write.csv(brecruit6, "output.sapling/bci_recruit6_forfig.csv")

## scatterplot of fitted (predicted) values ####
f2 <- ggplot(brecruit6, aes(x=max.ht, y=rr.mean)) + geom_point(aes(x=max.ht, y=rr.mean), data=brecruit6, shape=1) + ylab(expression(Sapling~recruitment~(no.~y^-1))) + 
  xlab(expression(H[max]~(m))) +
  geom_linerange(aes(ymax=fit.max, ymin=fit.min))

f32 <- f2 + theme_classic(base_size = 9) + scale_y_log10(breaks=c(0,0.05,0.1,0.5,1)) + scale_x_continuous(breaks=c(0,10,20,30)) +
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))

png(file="figures.sapling/bci-lmm1_recruit.vs.maxht_avg_pred.png", width = 82, height = 82, units = "mm", res=1200)
f32
dev.off()

## scatterplot of actual values ####
f3 <- ggplot(brecruit6, aes(x=max.ht, y=rr.mean)) + geom_point(aes(x=max.ht, y=rr.mean), data=brecruit6, shape=1) + ylab(expression(Sapling~recruitment~(no.~y^-1))) + 
  xlab(expression(H[max]~(m))) +
  geom_smooth(aes(x=max.ht, y=fit.mean), color="black", size=0.75)

f33 <- f3 + theme_classic(base_size = 9) + scale_y_log10(breaks=c(0,0.05,0.1,0.5,1)) + scale_x_continuous(breaks=c(0,10,20,30)) +
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        plot.title=element_text(hjust=0, size=9) ) + #plot.title for labelling subplots
  ggtitle('(c)')

f33
png(file="figures.sapling/bci-lmm1_recruit.vs.maxht_avg_actual_trendline.png", width = 82, height = 82, units = "mm", res=1200)
f33
dev.off()

# f33
# ggsave("lfdp-lmm1_recruit.vs.maxht_avg_actual.eps", path = "figures.sapling", width = 82, height = 82, units = "mm", dpi=1200) #save as eps

