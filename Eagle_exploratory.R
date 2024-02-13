homewd <- "/Users/emilyruhs/Desktop/"
setwd(homewd)

library(ggplot2)

data <- read.csv(file = paste0(homewd,"BaldEagle_EmilyGavin_Master_13Feb.csv"))

names(data)

# data2 <- subset(data, AnalyteName!="Perflouro-1-Octanesulfonate")
# 
# ggplot(data=data2)+geom_point(aes(x=SampleYear, y=ResultNDZero))+theme_bw()+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_grid(StudyArea~AnalyteName)

data3 <- subset(data, StudyArea=="GBLM") #only green bay area

ggplot(data=data3)+geom_point(aes(x=SampleYear, y=ResultNDZero))+theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(StudyArea~AnalyteName)


### Make the toxicant data the long format
#start with 2023
library(reshape2)
data2023 <- subset(data, year==2023)
subset.2023.dat = subset(data2023, select = c(1,11,27:65))
str(subset.2023.dat)
library(dplyr)
subset.2023.dat$field_id <- as.factor(subset.2023.dat$field_id)
subset.2023.dat$county <- as.factor(subset.2023.dat$county)
subset.2023.dat <- subset.2023.dat %>% mutate_if(is.character, as.numeric)
str(subset.2023.dat)

long.2023.dat <- melt(subset.2023.dat, by="field_id")
head(long.2023.dat)
names(long.2023.dat)[names(long.2023.dat)=="variable"] <- "toxicant"
names(long.2023.dat)[names(long.2023.dat)=="value"] <- "value"
head(long.2023.dat)

long.2023.dat$county2 <- factor(long.2023.dat$county, levels=c("WINNEBAGO", "OUTAGAMIE", "ASHLAND", "DOUGLAS", "BAYFIELD", "BROWN"))

p1 <- ggplot(long.2023.dat) + 
  geom_boxplot(aes(x=county2, y=value, fill=county2)) + theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),
        strip.background = element_rect(fill="white"), axis.text.y = element_blank(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p1)
#there are a lot of contaminants that are not detected, so let's remove those
head(long.2023.dat)
str(long.2023.dat)
unique(long.2023.dat$toxicant)
long.2023.dat2 = subset(long.2023.dat, toxicant =="FTSA10_2"|toxicant =="FTSA8_2"|toxicant =="FOSA"|toxicant =="N_ETFOSAA"|
                          toxicant =="N_MEFOSAA"|toxicant =="PFBA"|toxicant =="PFDA"|toxicant =="PFDOA"|toxicant =="PFDS"|
                          toxicant =="PFHPA"|toxicant =="PFHPS"|toxicant =="PFHXS"|toxicant =="PFNA"|toxicant =="PFNS"|
                          toxicant =="PFOA"|toxicant =="PFOS"|toxicant =="PFTEDA"|toxicant =="PFTRDA"|toxicant =="PFUNDA"|
                          toxicant =="TOTAL_PFAS")
unique(long.2023.dat2$toxicant)
p1 <- ggplot(long.2023.dat2) + 
  geom_boxplot(aes(x=county2, y=value, fill=county2)) + theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p1)

ggsave(file = "Fig4/Output-figures/mountainplot_estpro.pdf",
       plot=mountain.plot,
       units="mm",
       width=70,
       height=90,
       scale=3,
       dpi=300)
