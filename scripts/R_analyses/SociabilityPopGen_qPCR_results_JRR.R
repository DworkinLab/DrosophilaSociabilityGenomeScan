
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)

setwd("../data/")

#First I loaded in the data from each experiment and cleaned it, at the bottom I bind them together. 

###### Sec5

sec5_data <- read_excel('JR_GOI_Sec5_7142025_data.xls',sheet=1,range='A8:L95') #import data without metadata

sec5_data <- subset(sec5_data, Task == "UNKNOWN") #omit no Reverse Transcriptase controls and Master mix only controls

head(sec5_data)
dim(sec5_data)
str(sec5_data)

sec5_data$Cт <- as.numeric(sec5_data$Cт)

sec5_data <- sec5_data %>% 
  select(Well, 'Sample Name', 'Target Name', Cт)

#Adds a Reference CT column calculated 
sec5_data <- sec5_data %>% 
  group_by(`Sample Name`) %>%
  mutate(Cт_ref=mean(Cт[`Target Name` %in% c("elf1A", "Rap2")], na.rm = TRUE)) %>%
  ungroup()

sec5_data$delta_Cт <- sec5_data$Cт_ref - sec5_data$Cт #Calculate Delta CT

sec5_data <- subset(sec5_data, `Target Name` == "Sec5") #Remove housekeeping genes

sec5_data$Treatment = substr(sec5_data$`Sample Name`,1,nchar(sec5_data$`Sample Name`)-4) #Add new column with the same values for biological replicates 

ggplot(sec5_data, aes(x = Treatment, y = delta_Cт)) +
  geom_boxplot()

ggplot(sec5_data, aes(x = Treatment, y = delta_Cт, colour = `Sample Name`)) +
  geom_point()

sec5_data$Treatment <- as.factor(sec5_data$Treatment)
summary(sec5_data)

####### Cpsf100

cpsf100_data <- read_excel('JR_Cpsf100_7142025_data.xls',sheet=1,range='A8:L87') #import data without metadata

cpsf100_data <- subset(cpsf100_data, Task == "UNKNOWN") #omit no Reverse Transcriptase controls and Master mix only controls

head(cpsf100_data)
dim(cpsf100_data)
str(cpsf100_data)

cpsf100_data$Cт <- as.numeric(cpsf100_data$Cт)

cpsf100_data <- cpsf100_data %>% 
  select(Well, 'Sample Name', 'Target Name', Cт)

#Adds a Reference CT column calculated 
cpsf100_data <- cpsf100_data %>% 
  group_by(`Sample Name`) %>%
  mutate(Cт_ref=mean(Cт[`Target Name` %in% c("elf1A", "Rap2")], na.rm = TRUE)) %>%
  ungroup()

cpsf100_data$delta_Cт <- cpsf100_data$Cт_ref - cpsf100_data$Cт #Calculate Delta CT

cpsf100_data <- subset(cpsf100_data, `Target Name` == "Cpsf100") #Remove housekeeping genes

cpsf100_data$Treatment = substr(cpsf100_data$`Sample Name`,1,nchar(cpsf100_data$`Sample Name`)-4) #Add new column with the same values for biological replicates 

ggplot(cpsf100_data, aes(x = Treatment, y = delta_Cт)) +
  geom_boxplot()

ggplot(cpsf100_data, aes(x = Treatment, y = delta_Cт, colour = `Sample Name`)) +
  geom_point()

cpsf100_data$Treatment <- as.factor(cpsf100_data$Treatment)
summary(cpsf100_data)

####### Ptp99A

ptp99A_data <- read_excel('JR_Ptp99A_5122025_data.xls',sheet=1,range='A8:L94') #import data without metadata

ptp99A_data <- subset(ptp99A_data, Task == "UNKNOWN") #omit no Reverse Transcriptase controls and Master mix only controls

head(ptp99A_data)
dim(ptp99A_data)
str(ptp99A_data)

ptp99A_data$Cт <- as.numeric(ptp99A_data$Cт)

ptp99A_data <- ptp99A_data %>% 
  select(Well, 'Sample Name', 'Target Name', Cт)

#Adds a Reference CT column calculated 
ptp99A_data <- ptp99A_data %>% 
  group_by(`Sample Name`) %>%
  mutate(Cт_ref=mean(Cт[`Target Name` %in% c("elf1A", "Rap2")], na.rm = TRUE)) %>%
  ungroup()

ptp99A_data$delta_Cт <- ptp99A_data$Cт_ref - ptp99A_data$Cт #Calculate Delta CT

ptp99A_data <- subset(ptp99A_data, `Target Name` == "Ptp99A") #Remove housekeeping genes

ptp99A_data$Treatment = substr(ptp99A_data$`Sample Name`,1,nchar(ptp99A_data$`Sample Name`)-4) #Add new column with the same values for biological replicates 

ggplot(ptp99A_data, aes(x = Treatment, y = delta_Cт)) +
  geom_boxplot()

ggplot(ptp99A_data, aes(x = Treatment, y = delta_Cт, colour = `Sample Name`)) +
  geom_point()

ptp99A_data$Treatment <- as.factor(ptp99A_data$Treatment)
summary(ptp99A_data)

####### dpr8

dpr8_data <- read_excel('JR_5132025_GO_dpr8_data.xls',sheet=1,range='A8:L92') #import data without metadata

dpr8_data <- subset(dpr8_data, Task == "UNKNOWN") #omit no Reverse Transcriptase controls and Master mix only controls

head(dpr8_data)
dim(dpr8_data)
str(dpr8_data)

dpr8_data$Cт <- as.numeric(dpr8_data$Cт)

dpr8_data <- dpr8_data %>% 
  select(Well, 'Sample Name', 'Target Name', Cт)

#Adds a Reference CT column calculated 
dpr8_data <- dpr8_data %>% 
  group_by(`Sample Name`) %>%
  mutate(Cт_ref=mean(Cт[`Target Name` %in% c("elf1A", "Rap2")], na.rm = TRUE)) %>%
  ungroup()

dpr8_data$delta_Cт <- dpr8_data$Cт_ref - dpr8_data$Cт #Calculate Delta CT

dpr8_data <- subset(dpr8_data, `Target Name` == "dpr8") #Remove housekeeping genes

dpr8_data$Treatment = substr(dpr8_data$`Sample Name`,1,nchar(dpr8_data$`Sample Name`)-4) #Add new column with the same values for biological replicates 

ggplot(dpr8_data, aes(x = Treatment, y = delta_Cт)) +
  geom_boxplot()

ggplot(dpr8_data, aes(x = Treatment, y = delta_Cт, colour = `Sample Name`)) +
  geom_point()

dpr8_data$Treatment <- as.factor(dpr8_data$Treatment)
summary(dpr8_data)

####### kek6

kek6_data <- read_excel('JR_5142025_GOI_kek6_data.xls',sheet=1,range='A8:L93') #import data without metadata

kek6_data <- subset(kek6_data, Task == "UNKNOWN") #omit no Reverse Transcriptase controls and Master mix only controls

head(kek6_data)
dim(kek6_data)
str(kek6_data)

kek6_data$Cт <- as.numeric(kek6_data$Cт)

kek6_data <- kek6_data %>% 
  select(Well, 'Sample Name', 'Target Name', Cт)

#Adds a Reference CT column calculated 
kek6_data <- kek6_data %>% 
  group_by(`Sample Name`) %>%
  mutate(Cт_ref=mean(Cт[`Target Name` %in% c("elf1A", "Rap2")], na.rm = TRUE)) %>%
  ungroup()

kek6_data$delta_Cт <- kek6_data$Cт_ref - kek6_data$Cт #Calculate Delta CT

kek6_data <- subset(kek6_data, `Target Name` == "kek6") #Remove housekeeping genes

kek6_data$Treatment = substr(kek6_data$`Sample Name`,1,nchar(kek6_data$`Sample Name`)-4) #Add new column with the same values for biological replicates 

ggplot(kek6_data, aes(x = Treatment, y = delta_Cт)) +
  geom_boxplot()

ggplot(kek6_data, aes(x = Treatment, y = delta_Cт, colour = `Sample Name`)) +
  geom_point()

kek6_data$Treatment <- as.factor(kek6_data$Treatment)
summary(kek6_data)

####### rg

rg_data <- read_excel('rg_JRR_5122025_data.xls',sheet=1,range='A8:L92') #import data without metadata

rg_data <- subset(rg_data, Task == "UNKNOWN") #omit no Reverse Transcriptase controls and Master mix only controls

head(rg_data)
dim(rg_data)
str(rg_data)

rg_data$Cт <- as.numeric(rg_data$Cт)

rg_data <- rg_data %>% 
  select(Well, 'Sample Name', 'Target Name', Cт)

#Adds a Reference CT column calculated 
rg_data <- rg_data %>% 
  group_by(`Sample Name`) %>%
  mutate(Cт_ref=mean(Cт[`Target Name` %in% c("elf1A", "Rap2")], na.rm = TRUE)) %>%
  ungroup()

rg_data$delta_Cт <- rg_data$Cт_ref - rg_data$Cт #Calculate Delta CT

rg_data <- subset(rg_data, `Target Name` == "rg") #Remove housekeeping genes

rg_data$Treatment = substr(rg_data$`Sample Name`,1,nchar(rg_data$`Sample Name`)-4) #Add new column with the same values for biological replicates 

ggplot(rg_data, aes(x = Treatment, y = delta_Cт)) +
  geom_boxplot()

ggplot(rg_data, aes(x = Treatment, y = delta_Cт, colour = `Sample Name`)) +
  geom_point()

rg_data$Treatment <- as.factor(rg_data$Treatment)
summary(rg_data)

####### axo

axo_data <- read_excel('JR_05152025_GOI_axo_data.xls',sheet=1,range='A8:L95') #import data without metadata

axo_data <- subset(axo_data, Task == "UNKNOWN") #omit no Reverse Transcriptase controls and Master mix only controls

head(axo_data)
dim(axo_data)
str(axo_data)

axo_data$Cт <- as.numeric(axo_data$Cт)

axo_data <- axo_data %>% 
  select(Well, 'Sample Name', 'Target Name', Cт)

#Adds a Reference CT column calculated 
axo_data <- axo_data %>% 
  group_by(`Sample Name`) %>%
  mutate(Cт_ref=mean(Cт[`Target Name` %in% c("elf1A", "Rap2")], na.rm = TRUE)) %>%
  ungroup()

axo_data$delta_Cт <- axo_data$Cт_ref - axo_data$Cт #Calculate Delta CT

axo_data <- subset(axo_data, `Target Name` == "axo") #Remove housekeeping genes

axo_data$Treatment = substr(axo_data$`Sample Name`,1,nchar(axo_data$`Sample Name`)-4) #Add new column with the same values for biological replicates 

ggplot(axo_data, aes(x = Treatment, y = delta_Cт)) +
  geom_boxplot()

ggplot(axo_data, aes(x = Treatment, y = delta_Cт, colour = `Sample Name`)) +
  geom_point()

axo_data$Treatment <- as.factor(axo_data$Treatment)
summary(axo_data)

####### nlg1

nlg1_data <- read_excel('JR_5132025_GO_nlg1_data.xls',sheet=1,range='A8:L96') #import data without metadata

nlg1_data <- subset(nlg1_data, Task == "UNKNOWN") #omit no Reverse Transcriptase controls and Master mix only controls

head(nlg1_data)
dim(nlg1_data)
str(nlg1_data)

nlg1_data$Cт <- as.numeric(nlg1_data$Cт)

nlg1_data <- nlg1_data %>% 
  select(Well, 'Sample Name', 'Target Name', Cт)

#Adds a Reference CT column calculated 
nlg1_data <- nlg1_data %>% 
  group_by(`Sample Name`) %>%
  mutate(Cт_ref=mean(Cт[`Target Name` %in% c("elf1A", "Rap2")], na.rm = TRUE)) %>%
  ungroup()

nlg1_data$delta_Cт <- nlg1_data$Cт_ref - nlg1_data$Cт #Calculate Delta CT

nlg1_data <- subset(nlg1_data, `Target Name` == "nlg1") #Remove housekeeping genes

nlg1_data$Treatment = substr(nlg1_data$`Sample Name`,1,nchar(nlg1_data$`Sample Name`)-4) #Add new column with the same values for biological replicates 

ggplot(nlg1_data, aes(x = Treatment, y = delta_Cт)) +
  geom_boxplot()

ggplot(nlg1_data, aes(x = Treatment, y = delta_Cт, colour = `Sample Name`)) +
  geom_point()

nlg1_data$Treatment <- as.factor(nlg1_data$Treatment)
summary(nlg1_data)

#################################

all_data <- rbind(axo_data, cpsf100_data, dpr8_data,
                  kek6_data, nlg1_data, ptp99A_data,
                  rg_data, sec5_data)

all_data$Gene <- all_data$`Target Name` #Added Gene collumn
all_data$Sex <- ifelse(grepl("F", all_data$`Sample Name`), "F", "M") #Added Sex collumn 

head(all_data)

model1 <- lm(delta_Cт ~ Treatment * Gene * Sex, 
             data = all_data)

emm1 <- emmeans(model1, ~ Treatment | Gene + Sex)

contrasts1 <- contrast(emm1, 
                      method = "pairwise",
                      by = c("Gene", "Sex"))

summary(contrasts1)

contrast_df <- as.data.frame(contrasts1)

ggplot(contrast_df, aes(x = contrast, y = estimate, colour = Sex)) +
  geom_point() +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), width = 0.2) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.25) + 
  facet_grid(~ Gene, scales = "free_y") +
  labs(
    x = "",
    y = "Estimated ∆∆CT (control-RNAi)",
    title = "qPCR Results by Gene and Sex"
  ) + theme_light() + 
  theme(axis.text.x = element_blank())

