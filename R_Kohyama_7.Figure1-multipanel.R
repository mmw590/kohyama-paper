#' ---
#' title: "R_Kohyama_7.Figure1-multipanel.R"
#' author: "Maria Wang"
#' date: "2017-03-28"
#' description: "This script makes a 2 x 6 multipanel figure (Figure 1 for publication)"
#' ---


# ~ Loading required R packages ----
library(dplyr)
library(ggplot2)
library(scales)
library(gridExtra)

# ~ Set working directory ----
setwd("C:/Users/wang/Dropbox/LFDP/data analysis for Seth/Kohyama Hypothesis")
getwd()

# ~ Load required input files----
sshm2 <- read.csv("intermediate data/ldfp_sts_forglm.csv")
str(sshm2)
sshm2.sp <- sshm2 %>% group_by(sp, seed.mass, max.ht) %>% summarize(n.new.sling=sum(n.new.sling), nseeds = sum(nseeds), sts.mean=mean(sts))
sshm2.sp[sshm2.sp$sts.mean > 1, ]$sts.mean <- 1 #max sts ratio is 1
sshm2.sp[sshm2.sp$sts.mean == 0, ]$sts.mean <- 1*10^-4 #add a tiny number to 0 so that log10(0) doesn't turn the value into -Inf

slinghm <- read.csv("intermediate data/lfdp_seedlinghm_forglm.csv")
slinghm.sp <- slinghm %>% group_by(sp, seed.mass, max.ht, cohort, plotno) %>% summarize(live.in.y1=sum(live.in.y1), nsling = n(), slingsurvival.y1=live.in.y1/nsling) %>% group_by(sp, seed.mass, max.ht) %>% summarize(sl.mean = mean(slingsurvival.y1))
unique(slinghm$sp)

recruit5 <- read.csv("output.sapling/lfdp_recruit5.csv")
recruit6 <- recruit5 %>% group_by(sp, max.ht) %>% summarize(rr.mean = mean(rate), rr.min = min(rate), rr.max=max(rate))

brecruit5 <- read.csv("output.sapling/bci_recruit5.csv")
brecruit6 <- brecruit5 %>% group_by(sp, max.ht) %>% summarize(rr.mean = mean(rate), rr.min = min(rate), rr.max=max(rate))


#### Fig 1(a) Actual seed survival ~ seed mass + shaded by hmax FOR PUBLICATION ####
a3 <- ggplot(sshm2.sp, aes(x=seed.mass, y=sts.mean, fill=max.ht)) + geom_point(aes(x=seed.mass, y=sts.mean), data=sshm2.sp, shape=21, size=3) + ylab("STS ratio") + xlab("Seed mass (g)") 

range(sshm2.sp$sts.mean)

a31 <- a3 + scale_x_log10(breaks = c(0, 0.001, 0.01,0.1,1), labels=c(NA, "0.001", "0.01", "0.1", "1")) + 
  scale_y_log10(breaks = c(0, 10^-4, 10^-2,1), labels=c(NA, "0", "0.01", "1")) + 
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

#### Fig1(b) Actual seedling survival ~ seed mass + shaded by hmax FOR PUBLICATION ####
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

# Fig1(c). LFDP sapling rec vs height for PUBLICATION ####
c1 <- ggplot(recruit6, aes(x=max.ht, y=rr.mean)) + geom_point(aes(x=max.ht, y=rr.mean), data=recruit6, shape=1) + ylab(expression(Sapling~recruitment~(no.~y^-1))) + 
  xlab(expression(H[max]~(m)))

c31 <- c1 + theme_classic(base_size = 9) + scale_y_log10(breaks=c(0,0.05,0.1,0.5,1)) + scale_x_continuous(breaks=c(0,10,20,30)) +
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        plot.title=element_text(hjust=0, size=9) ) + #plot.title for labelling subplots
  ggtitle('(c)')

c31

#### Fig 1(d) BCI - seed to seedling ratio ~ seed mass + shaded by hmax PLACEHOLDER ####
d3 <- ggplot(sshm2, aes(x=seed.mass, y=sts.mean, fill=max.ht)) + geom_point(aes(x=seed.mass, y=sts.mean), data=sshm2, shape=21, size=3) + ylab("STS ratio") + xlab("Seed mass (g)") 

d31 <- d3 + scale_x_log10(breaks = c(0, 0.001, 0.01,0.1,1), labels=c(NA, "0.001", "0.01", "0.1", "1")) + 
  scale_y_log10(breaks = c(0, 10^-4, 10^-2,1), labels=c(NA, "0", "0.01", "1")) + 
  theme_classic(base_size = 9) + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), #line size
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"), #adjust margins
        legend.position="none", #no legends
        plot.title=element_text(hjust=0, size=9) ) + #plot.title for labelling subplots
  scale_fill_gradient(low=gray(1), high =gray(0)) + 
  annotation_logticks() +
  ggtitle('(d) ')

d31 <- ggplot()

#### Fig1(e) BCI seedling survival ~ seed mass + shaded by hmax PLACEHOLDER  ####
e3 <- ggplot(slinghm.sp, aes(x=seed.mass, y=sl.mean, fill=max.ht)) + geom_point(aes(x=seed.mass, y=sl.mean), data=slinghm.sp, shape=21, size=3) + ylab("Seedling survival (proportion)") + xlab("Seed mass (g)") 

e31 <- e3 + scale_x_log10(breaks = c(0, 0.001, 0.01,0.1,1), labels=c(NA, "0.001", "0.01", "0.1", "1")) + 
  scale_y_log10(breaks = c(0, 0.01, 0.1, 1), labels=c("0", "0.01", "0.1", "1")) + 
  theme_classic(base_size = 9) + 
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        legend.position="none",
        plot.title=element_text(hjust=0, size=9) ) +
  scale_fill_gradient(low=gray(1), high =gray(0)) + 
  annotation_logticks() +
  ggtitle('(e)  ')

e31 <- ggplot()


# Fig1(f). BCI sapling rec vs height for PUBLICATION ####

f1 <- ggplot(brecruit6, aes(x=max.ht, y=rr.mean)) + geom_point(aes(x=max.ht, y=rr.mean), data=brecruit6, shape=1) + ylab(expression(Sapling~recruitment~(no.~y^-1))) + 
  xlab(expression(H[max]~(m)))

f31 <- f1 + theme_classic(base_size = 9) + scale_y_log10(breaks=c(0,0.05,0.1,0.5,1)) + scale_x_continuous(breaks=c(0,10,20,30)) +
  theme(panel.border = element_rect(fill=NA, color="black"), 
        line=element_line(size=0.25), 
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        plot.title=element_text(hjust=0, size=9) ) + #plot.title for labelling subplots
  ggtitle('(f)')

f31



#### Save as eps ####
library(RGraphics)

png(file="Fig1a-c_LFDP_seed.seedling.sapling.png",width=6.27,height=2.25, units="in", res=300)
grid.arrange(a31, b31, c31, top=textGrob("(a)-(c) Luquillo Forest Dynamics Plot (LFDP)", gp=gpar(fontsize=9)), ncol=3)
dev.off()

png(file="Fig1d-f_BCI_seed.seedling.sapling.png",width=6.27,height=2.25, units="in", res=300)
grid.arrange(d31, e31, f31, top=textGrob("(d)-(f) Barro Colorado Island (BCI)", gp=gpar(fontsize=9)), ncol=3)
dev.off()

grid.arrange(a31, b31, c31, d31, e31, f31, ncol=3, nrow=2)
ggsave("Fig1_2by6.eps", path = "figures.seed", width =6.27, height = 4.5, units = "in", dpi=800) #save as eps, doesn't work

### Fig1a-c eps #
grid.arrange(a31, b31, c31, top=textGrob("(a)-(c) Luquillo Forest Dynamics Plot (LFDP)", gp=gpar(fontsize=9)), ncol=3)
gtop <- arrangeGrob(a31, b31, c31, top=textGrob("(a)-(c) Luquillo Forest Dynamics Plot (LFDP)", gp=gpar(fontsize=9)), ncol=3)
ggsave("Fig1a-c_LFDP_seed.seedling.sapling.eps", gtop, path = "figures.seed", width=6.27, height = 2.25, units = "in", dpi=800) #save as eps

### Fig1d-f eps #
grid.arrange(d31, e31, f31, top=textGrob("(d)-(f) Barro Colorado Island (BCI)", gp=gpar(fontsize=9)), ncol=3)
gbot <- arrangeGrob(d31, e31, f31, top=textGrob("(d)-(f) Barro Colorado Island (BCI)", gp=gpar(fontsize=9)), ncol=3)
ggsave("Fig1d-f_BCI_seed.seedling.sapling.eps", gbot, path = "figures.seed", width =6.27, height = 2.25, units = "in", dpi=800) #save as eps

?ggsave

