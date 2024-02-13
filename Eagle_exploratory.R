homewd <- "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/eagle-contaminant/"
setwd(homewd)

library(ggplot2)

data <- read.csv(file = paste0(homewd,"BaldEagle_EmilyGavin_Master_13Feb.csv"))

names(data)

# data2 <- subset(data, AnalyteName!="Perflouro-1-Octanesulfonate")
# 
# ggplot(data=data2)+geom_point(aes(x=SampleYear, y=ResultNDZero))+theme_bw()+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_grid(StudyArea~AnalyteName)

# data3 <- subset(data, StudyArea=="GBLM") #only green bay area
# 
# ggplot(data=data3)+geom_point(aes(x=SampleYear, y=ResultNDZero))+theme_bw()+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_grid(StudyArea~AnalyteName)


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

ggsave(file = "figures/2023-contaminants.pdf",
       plot=p1,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)

#do any of these contaminants correlate with immune function variables?
#need to add back in the immune data
lab.2023.dat = subset(data2023, select = c(1,66,68,70,72:74,84:90))
#ok remove a couple
lab.2023.dat = subset(lab.2023.dat, select = -c(2))
all.2023.dat <- merge(long.2023.dat2, lab.2023.dat, by="field_id")
str(all.2023.dat)

p2 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=T4_div_T3, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=T4_div_T3), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p2)

ggsave(file = "figures/2023-T4_div_T3.pdf",
       plot=p2,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)

mod2 <- lm(T4_div_T3 ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod2) #nothing

p3 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=Clysis, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=Clysis), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p3)

ggsave(file = "figures/2023-Clysis.pdf",
       plot=p3,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)

mod3 <- lm(Clysis ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod3) #nothing


p4 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=IgY_final, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=IgY_final), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p4)

ggsave(file = "figures/2023-IgY.pdf",
       plot=p4,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)

mod4 <- lm(IgY_final ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod4)

p5 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=H_L, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=H_L), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p5)

ggsave(file = "figures/2023-HL.pdf",
       plot=p5,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)

mod5 <- lm(H_L ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod5) #some marginal significance - close to 0.1


p6 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=Lperc, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=Lperc), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p6)

ggsave(file = "figures/2023-Lperc.pdf",
       plot=p6,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)

mod6 <- glm(Lperc ~ toxicant:value, family="poisson",data=all.2023.dat) #should probably put sex and age in the model
summary(mod6) #nothing


p7 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=Hperc, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=Hperc), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p7)

ggsave(file = "figures/2023-Hperc.pdf",
       plot=p7,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)

mod7 <- glm(Hperc ~ toxicant:value, family="poisson",data=all.2023.dat) #should probably put sex and age in the model
summary(mod7) #nothing



p8 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=Total_count, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=Total_count), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p8)

ggsave(file = "figures/2023-Total_count.pdf",
       plot=p8,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)

mod8 <- glm(Total_count ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod8) #almost all if poisson, some marginal if gaussian






##### OKAY NOW LETS TRY TO LOOK AT PATTERNS OVER TIME FOR EACH TOXICANT. 
### Make the toxicant data the long format
#start with 2023
library(reshape2)
subset.dat = subset(data, select = c(1,9,11,27:65))
str(subset.dat)
library(dplyr)
subset.dat$field_id <- as.factor(subset.dat$field_id)
subset.dat$county <- as.factor(subset.dat$county)
subset.dat$year <- as.factor(subset.dat$year)
subset.dat <- subset.dat %>% mutate_if(is.character, as.numeric)
str(subset.dat)

long.dat <- melt(subset.dat, by="field_id")
head(long.dat)
names(long.dat)[names(long.dat)=="variable"] <- "toxicant"
names(long.dat)[names(long.dat)=="value"] <- "value"
head(long.dat)

long.dat2 = subset(long.dat, toxicant =="FTSA10_2"|toxicant =="FTSA8_2"|toxicant =="FOSA"|toxicant =="N_ETFOSAA"|
                          toxicant =="N_MEFOSAA"|toxicant =="PFBA"|toxicant =="PFDA"|toxicant =="PFDOA"|toxicant =="PFDS"|
                          toxicant =="PFHPA"|toxicant =="PFHPS"|toxicant =="PFHXS"|toxicant =="PFNA"|toxicant =="PFNS"|
                          toxicant =="PFOA"|toxicant =="PFOS"|toxicant =="PFTEDA"|toxicant =="PFTRDA"|toxicant =="PFUNDA"|
                          toxicant =="TOTAL_PFAS")

#this dataset needs to be filled in more
p9 <- ggplot(long.dat) + 
  geom_point(aes(x=county, y=value, color=year)) + theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p9)

ggsave(file = "figures/2023-contaminants.pdf",
       plot=p1,
       units="mm",
       width=80,
       height=50,
       scale=3,
       dpi=300)
