
library(skimr) #data summary package
library(ggplot2)

homewd <- "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/eagle-contaminant/"
setwd(homewd)
data <- read.csv(file = paste0(homewd,"working-data/BaldEagle_EmilyGavin_Master_13Feb.csv"))

names(data)

#dev.new()
skim(data) #summarizes the data
#data2 <- subset(data, AnalyteName!="Perflouro-1-Octanesulfonate")
# 
# ggplot(data=data2)+geom_point(aes(x=SampleYear, y=ResultNDZero))+theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_grid(StudyArea~AnalyteName)

# data3 <- subset(data, StudyArea=="GBLM") #only green bay area
# 
# ggplot(data=data3)+geom_point(aes(x=SampleYear, y=ResultNDZero))+theme_bw()+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   facet_grid(StudyArea~AnalyteName)

#add in a body condition index
modBaseF <- lm(log10(mass_kg)~log10(foot_pad_mm), data=subset(data, DNA_sex=="F"))
summary(modBaseF)

modBaseM <- lm(log10(mass_kg)~log10(foot_pad_mm), data=subset(data, DNA_sex=="M"))
summary(modBaseM)

data$predicted_mass_kg <- NA
data$predicted_mass_kg[data$DNA_sex=="F" & !is.na(data$mass_kg) &!is.na(data$foot_pad_mm)] <- 10^predict(modBaseF) #you are predicting log10 mass, so you need to exponent it here
data$predicted_mass_kg[data$DNA_sex=="M"& !is.na(data$mass_kg) &!is.na(data$foot_pad_mm)] <- 10^predict(modBaseM)





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

# long.2023.dat <- melt(subset.2023.dat, by="field_id")
# head(long.2023.dat)
# names(long.2023.dat)[names(long.2023.dat)=="variable"] <- "toxicant"
# names(long.2023.dat)[names(long.2023.dat)=="value"] <- "value"
# head(long.2023.dat)

subset.2023.dat$carboxyl_sum <- (subset.2023.dat$PFBA) + (subset.2023.dat$PFHXA)+ 
  (subset.2023.dat$PFHXDA) + (subset.2023.dat$PFOA) + (subset.2023.dat$PFDA) +
  (subset.2023.dat$PFDOA) + (subset.2023.dat$PFODA) + (subset.2023.dat$PFPEA) +
  (subset.2023.dat$PFUNDA) + (subset.2023.dat$PFTRDA) + (subset.2023.dat$PFTEDA) +
  (subset.2023.dat$PFNA) + (subset.2023.dat$DONA) + (subset.2023.dat$PFHPA) +
  (subset.2023.dat$PFMBA) + (subset.2023.dat$PFMPA)

subset.2023.dat$sulfuric_sum <- (subset.2023.dat$FTSA10_2) +(subset.2023.dat$PF3OUDS11CL)+
  (subset.2023.dat$FTSA4_2) + (subset.2023.dat$FTSA6_2) + (subset.2023.dat$CL_9PF3ONS) +
  (subset.2023.dat$PFBS) + (subset.2023.dat$PFDS) + (subset.2023.dat$PFHPS) + 
  (subset.2023.dat$PFHXS) + (subset.2023.dat$PFNS) + (subset.2023.dat$PFOS) +
  (subset.2023.dat$PFPES) + (subset.2023.dat$FTSA8_2)

subset.2023.dat$phosphate_sum <- (subset.2023.dat$DIPAP6_2) + (subset.2023.dat$DIPAP8_2)

subset.2023.dat$amine_sum <- (subset.2023.dat$FOSA) + (subset.2023.dat$N_ETFOSA) + 
  (subset.2023.dat$N_ETFOSAA) + (subset.2023.dat$N_ETFOSE) + (subset.2023.dat$N_MEFOSA)+
  (subset.2023.dat$N_MEFOSAA) + (subset.2023.dat$N_MEFOSE) 

subset.2023.dat$shortchain_sum <- (subset.2023.dat$PFMPA) + (subset.2023.dat$FTSA4_2) + 
  (subset.2023.dat$PFBA) + (subset.2023.dat$PFBS) + (subset.2023.dat$PFDS)+
  (subset.2023.dat$PFMBA) + (subset.2023.dat$PFHXA) + (subset.2023.dat$PFPEA) + (subset.2023.dat$PFPES) +
  (subset.2023.dat$DIPAP6_2) + (subset.2023.dat$FTSA6_2) +(subset.2023.dat$PFHXS) + 
  (subset.2023.dat$PFHPA) + (subset.2023.dat$PFHPS)

subset.2023.dat$longchain_sum <- (subset.2023.dat$PF3OUDS11CL) +(subset.2023.dat$DIPAP8_2) + (subset.2023.dat$FTSA8_2)+
  (subset.2023.dat$CL_9PF3ONS) + (subset.2023.dat$FOSA) + (subset.2023.dat$N_ETFOSA) + (subset.2023.dat$N_ETFOSAA) +
  (subset.2023.dat$N_ETFOSE) + (subset.2023.dat$PFOA) + (subset.2023.dat$PFOS) + (subset.2023.dat$DONA) +
  (subset.2023.dat$N_MEFOSA) + (subset.2023.dat$N_MEFOSAA) + (subset.2023.dat$N_MEFOSE) + (subset.2023.dat$PFDA) +
  (subset.2023.dat$PFNA) + (subset.2023.dat$PFNS) + (subset.2023.dat$FTSA10_2) + (subset.2023.dat$PFDOA) +
  (subset.2023.dat$PFUNDA) + (subset.2023.dat$PFTRDA) + (subset.2023.dat$PFTEDA) + (subset.2023.dat$PFHXDA) +
  (subset.2023.dat$PFODA)


long.2023.dat <- melt(subset.2023.dat, by="field_id")
head(long.2023.dat)
names(long.2023.dat)[names(long.2023.dat)=="variable"] <- "toxicant"
names(long.2023.dat)[names(long.2023.dat)=="value"] <- "value"
head(long.2023.dat)

long.2023.dat$county2 <- factor(long.2023.dat$county, levels=c("WINNEBAGO", "OUTAGAMIE", "ASHLAND", "DOUGLAS", "BAYFIELD", "BROWN"))

chemical.dat <- subset(long.2023.dat, toxicant=="carboxyl_sum" | toxicant=="sulfuric_sum" |
  toxicant=="phosphate_sum" | toxicant=="amine_sum")
#there are not really any phosphates, so remove
chemical.dat <- subset(chemical.dat, toxicant!="phosphate_sum")

chain.dat <- subset(long.2023.dat, toxicant=="shortchain_sum" | toxicant=="longchain_sum")

reg.dat <- subset(long.2023.dat, toxicant!="longchain_sum" & toxicant!="shortchain_sum"&toxicant!="carboxyl_sum"&
                    toxicant!="sulfuric_sum"& toxicant!="amine_sum"&toxicant!="phosphate_sum")


total.aov <- aov(TOTAL_PFAS~ county2, data=data2023)
summary(total.aov)
TukeyHSD(total.aov)
# ASHLAND-WINNEBAGO    198.33992   49.478071 347.201762 0.0049059
# DOUGLAS-WINNEBAGO    224.89525   25.176127 424.614373 0.0213410
# ASHLAND-OUTAGAMIE    180.21867   -8.078329 368.515662 0.0659901
# DOUGLAS-OUTAGAMIE    206.77400  -23.841779 437.389779 0.0962950
# BAYFIELD-ASHLAND    -112.49694 -234.042132   9.048243 0.0805328
# BROWN-ASHLAND       -127.34117 -266.985956  12.303623 0.0876675

data2023$county2 <- factor(data2023$county, levels=c("WINNEBAGO", "OUTAGAMIE", "ASHLAND", "DOUGLAS", "BAYFIELD", "BROWN"))
ptotal <- ggplot(data2023) + 
  geom_boxplot(aes(x=county2, y=log10(TOTAL_PFAS), fill=county2)) + theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),
        axis.text.y = element_text(size=14), legend.position = "left") + ggtitle("Total PFAS contamination (ugg)")
print(ptotal)

ggsave(file = "figures/2023-totalpfas.pdf",
       plot=ptotal,
       units="mm",
       width=50,
       height=40,
       scale=3,
       dpi=300)

#
p1 <- ggplot(reg.dat) + 
  geom_boxplot(aes(x=county2, y=value, fill=county2)) + theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(size=12), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 13, face = "italic"), strip.text.y = element_text(size = 13),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p1)
#there are a lot of contaminants that are not detected, so let's remove those
# head(long.2023.dat)
# str(long.2023.dat)
# unique(long.2023.dat$toxicant)
long.2023.dat2 = subset(long.2023.dat, toxicant =="FTSA10_2"|toxicant =="FTSA8_2"|toxicant =="FOSA"|toxicant =="N_ETFOSAA"|
                          toxicant =="N_MEFOSAA"|toxicant =="PFBA"|toxicant =="PFDA"|toxicant =="PFDOA"|toxicant =="PFDS"|
                          toxicant =="PFHPA"|toxicant =="PFHPS"|toxicant =="PFHXS"|toxicant =="PFNA"|toxicant =="PFNS"|
                          toxicant =="PFOA"|toxicant =="PFOS"|toxicant =="PFTEDA"|toxicant =="PFTRDA"|toxicant =="PFUNDA"|
                          toxicant =="TOTAL_PFAS")
# unique(long.2023.dat2$toxicant)
# p1 <- ggplot(long.2023.dat2) + 
#   geom_boxplot(aes(x=county2, y=value, fill=county2)) + theme_bw()+
#   theme(axis.title = element_blank(),axis.text.x = element_blank(),
#         strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
#         panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
#         strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
#         legend.position = "left") + facet_wrap(~toxicant, scales = "free")
# print(p1)
# 
ggsave(file = "figures/2023-contaminants-subgroups2.pdf",
       plot=p1,
       units="mm",
       width=90,
       height=40,
       scale=3,
       dpi=300)

aov1 <- aov(value~toxicant*county2, data=chain.dat)
summary(aov1)
TukeyHSD(aov1)
#longchain_sum:ASHLAND-longchain_sum:WINNEBAGO      197.08183333   81.553534  312.610133 0.0000299
#longchain_sum:ASHLAND-longchain_sum:OUTAGAMIE      179.42383333   33.290809  325.556857 0.0057287
#longchain_sum:DOUGLAS-longchain_sum:OUTAGAMIE      206.02350000   27.047828  384.999172 0.0123284
#longchain_sum:BROWN-longchain_sum:ASHLAND         -127.85313333 -236.228284  -19.477982 0.0092984
#longchain_sum:BAYFIELD-longchain_sum:DOUGLAS      -137.93911111 -277.850830    1.972607 0.0564659

#do any of these contaminants correlate with immune function variables?
#need to add back in the immune data
lab.2023.dat = subset(data2023, select = c(1,66,68,70,72:74,84:91))
#ok remove a couple
lab.2023.dat = subset(lab.2023.dat, select = -c(2))
all.2023.dat <- merge(long.2023.dat2, lab.2023.dat, by="field_id")
all.2023.dat <- merge(chain.dat, lab.2023.dat, by="field_id")
all.2023.dat <- merge(chemical.dat, lab.2023.dat, by="field_id")
str(all.2023.dat)

pCon <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=predicted_mass_kg, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=predicted_mass_kg), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(pCon)

mod_Con <- lm(predicted_mass_kg ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod_Con) #nothing

p2 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=T4_div_T3, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=T4_div_T3), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p2)

ggsave(file = "figures/2023-T4_div_T3_subgroups2.pdf",
       plot=p2,
       units="mm",
       width=80,
       height=40,
       scale=3,
       dpi=300)

mod_t3 <- lm(T3_conc ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod_t3) #nothing

mod_t4 <- lm(T4_conc ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod_t4) #nothing

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

ggsave(file = "figures/2023-Clysis_subgroups2.pdf",
       plot=p3,
       units="mm",
       width=80,
       height=40,
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

ggsave(file = "figures/2023-IgY_subgroups2.pdf",
       plot=p4,
       units="mm",
       width=80,
       height=40,
       scale=3,
       dpi=300)

mod4 <- lm(IgY_final ~ toxicant:value, data=all.2023.dat) #should probably put sex and age in the model
summary(mod4)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   0.6530873  0.0708443   9.219 1.34e-12 ***
# toxicantshortchain_sum:value -0.0578988  0.0228969  -2.529   0.0145 *  
# toxicantlongchain_sum:value  -0.0008514  0.0002862  -2.975   0.0044 ** 

p5 <- ggplot(all.2023.dat) + 
  geom_point(aes(x=value, y=H_L, shape=county)) + theme_bw()+
  geom_smooth(aes(x=value, y=H_L), color="black", linewidth=0.5, method=lm, data=all.2023.dat)+
  theme(axis.title = element_blank(), axis.text.x = element_text(),
        strip.background = element_rect(fill="white"), axis.text.y = element_text(), 
        panel.spacing.x  =  unit(c(.2), "cm"), panel.spacing.y  =  unit(c(0), "cm"),
        strip.text.x = element_text(size = 12, face = "italic"), strip.text.y = element_text(size = 12),
        legend.position = "left") + facet_wrap(~toxicant, scales = "free")
print(p5)

ggsave(file = "figures/2023-HL_2.pdf",
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

mod6 <- glm(Lperc ~ toxicant:value, family="quasibinomial",data=all.2023.dat) #should probably put sex and age in the model
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

mod7 <- glm(Hperc ~ toxicant:value, family="quasibinomial",data=all.2023.dat) #should probably put sex and age in the model
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

ggsave(file = "figures/2023-Total_count_subgroups_2.pdf",
       plot=p8,
       units="mm",
       width=80,
       height=40,
       scale=3,
       dpi=300)

mod8 <- glm(Total_count ~ toxicant:value, data=all.2023.dat, family=poisson) #should probably put sex and age in the model
summary(mod8) #almost all if poisson





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
