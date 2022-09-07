####Analysis of hybrid deterioration data####
#Updated 2022-09-05
#Analyses includes ANOVAs, general linear mixed effect models and general linear models
#R version 4.2.1 (2022-06-23)
#Packages used - car, DHARMa, emmeans and lme4

####Import Data ####
setwd("C:/TargetDirectory")
#Set to directory containing "Table S1.txt"
H2O2 <- read.delim("Table S1.txt")
#check data
str(H2O2)
#convert to factors
H2O2$Population <- as.factor(H2O2$Population)
#Which population (P1, P2, F1 or F2)
H2O2$Plate <- as.factor(H2O2$Plate)
#Replicate plates (R1-R4) for each population
H2O2$Well <- as.factor(H2O2$Well)
#Position on plate by population
H2O2$Position <- as.factor(H2O2$Position)
#Name for each population (P001-P960)

H2O2$Generation <- as.factor(H2O2$Generation)
#Which generation of adaptation/divergence - 100 or 1000 gens
str(H2O2)
H2O2$Survival <- as.factor(H2O2$Survival)
#Binary 1 = Alive, 0 = Dead
#Converted to survival if Growth >= 0.1
#See Growth_SurvivalByBlank.txt for survival based on blank
H2O2$Concentration <- as.factor(H2O2$Concentration)
#Change to factor to compare survival at each concentration

H2O2_G100 <- subset(H2O2, Generation == "100")
H2O2_G100 <- droplevels(H2O2_G100)

H2O2_G1K <- subset(H2O2, Generation == "1000")
H2O2_G1K <- droplevels(H2O2_G1K)
str(H2O2_G1K)

####Set Model####
#load required packages
library(lme4)
library(car)
library(DHARMa)
library(emmeans)

####Survival all replicates####
#####100 gens#####
#Generate models from complicated to simple
fm1 <- glmer(Survival ~ Population + Concentration + Population:Concentration + (1|Plate), data = H2O2_G100, family=binomial)
#make model - fixed effect = Gens, random effect = replicate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glmer(Survival ~ Population + Concentration + (1|Plate), data = H2O2_G100, family=binomial)

AIC(fm1)
#4306.8
AIC(fm2)
#4508.9
#fm1 has lowest AIC. Interaction kept

fm3 <- glmer(Survival ~ Population + Population:Concentration + (1|Plate), data = H2O2_G100, family=binomial)
fm4 <- glmer(Survival ~ Concentration + Population:Concentration + (1|Plate), data = H2O2_G100, family=binomial)
AIC(fm3)
#4306.8
AIC(fm4)
#4306.8
#fm3 & fm4 no change in AIC. both single effects dropped

fm5 <- glmer(Survival ~ Population:Concentration + (1|Plate), data = H2O2_G100, family=binomial)
AIC(fm5)
#4306.8
#No change in AIC vs. fm1
#fm5 simplest best model

summary(fm5)
Anova(fm5)
#chi-sq interaction = 1573.3, df = 35, p < 0.001

fm5a <- glm(Survival ~ Population:Concentration, data = H2O2_G100, family=binomial)
anova(fm5, fm5a)

#Get stats and check model. Currently using fm5, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm5)
plot(simulationOutput)
testDispersion(simulationOutput)

mmG100 <- emmeans(fm5, ~ Population:Concentration)
bonG100 <- contrast(mmG100, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG100, file = "Stats/Survival_stats_G100.txt", sep = "\t", row.names = F)
bonG100 <- read.delim("Stats/Survival_stats_G100.txt")
bonG100$Significance <- ifelse(bonG100$p.value <= 0.001, "***", ifelse(bonG100$p.value <= 0.01, "**", ifelse(bonG100$p.value <= 0.05, "*", "ns")))
write.table(bonG100, file = "Stats/Survival_stats_G100.txt", sep = "\t", row.names = F)

#####1000 gens#####

fm1 <- glmer(Survival ~ Population + Concentration + Population:Concentration + (1|Plate), data = H2O2_G1K, family=binomial)
#make model - fixed effect = Gens, random effect = replicate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glmer(Survival ~ Population + Concentration + (1|Plate), data = H2O2_G1K, family=binomial)

AIC(fm1)
#705.7
AIC(fm2)
#687.5
#fm2 has lowest AIC. Interaction dropped

fm3 <- glmer(Survival ~ Population + (1|Plate), data = H2O2_G1K, family=binomial)
fm4 <- glmer(Survival ~ Concentration + (1|Plate), data = H2O2_G1K, family=binomial)
AIC(fm3)
#3991.7
AIC(fm4)
#1399.7
#fm3 & fm4 both increase AIC. both single effects kept
#fm2 simplest best model

summary(fm2)
Anova(fm2)
#chi-sq population = 196.05, df = 3, p < 0.001
#chi-sq concentration = 303.72, df = 6, p < 0.001

fm2a <- glm(Survival ~ Population + Concentration, data = H2O2_G1K, family=binomial)
anova(fm2, fm2a)

#Get stats and check model. Currently using fm2, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm2)
plot(simulationOutput)
testDispersion(simulationOutput)

mmG1000 <- emmeans(fm2, ~ Population)
bonG1000 <- contrast(mmG1000, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG1000, file = "Stats/Survival_stats_G1000_Pop.txt", sep = "\t", row.names = F)
bonG1000 <- read.delim("Stats/Survival_stats_G1000_Pop.txt")
bonG1000$Significance <- ifelse(bonG1000$p.value <= 0.001, "***", ifelse(bonG1000$p.value <= 0.01, "**", ifelse(bonG1000$p.value <= 0.05, "*", "ns")))
write.table(bonG1000, file = "Stats/Survival_stats_G1000_Pop.txt", sep = "\t", row.names = F)

mmG1000 <- emmeans(fm2, ~ Concentration)
bonG1000con <- contrast(mmG1000, "pairwise", simple = "Concentration", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG1000con, file = "Stats/Survival_stats_G1000_Con.txt", sep = "\t", row.names = F)
bonG1000con <- read.delim("Stats/Survival_stats_G1000_Con.txt")
bonG1000con$Significance <- ifelse(bonG1000con$p.value <= 0.001, "***", ifelse(bonG1000con$p.value <= 0.01, "**", ifelse(bonG1000con$p.value <= 0.05, "*", "ns")))
write.table(bonG1000con, file = "Stats/Survival_stats_G1000_Con.txt", sep = "\t", row.names = F)


#####All gens#####
fm1 <- glmer(Survival ~ Population + Concentration + Generation + 
               Population:Concentration + Population:Generation + Concentration:Generation +
               Population:Concentration:Generation +
               (1|Plate), data = H2O2, family=binomial)
#make model - fixed effect = Gens, random effect = replicate plate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glmer(Survival ~ Population + Concentration + Generation + 
               Population:Concentration + Population:Generation + Concentration:Generation +
               (1|Plate), data = H2O2, family=binomial)

AIC(fm1)
#AIC = 5030.58
AIC(fm2)
#AIC = 5022.30
#fm2 has lowest AIC. Drop 3-way interaction

fm3 <- glmer(Survival ~ Population + Concentration + Generation + 
               Population:Concentration + Population:Generation +
               (1|Plate), data = H2O2, family=binomial)
fm4 <- glmer(Survival ~ Population + Concentration + Generation + 
               Population:Concentration + Concentration:Generation +
               (1|Plate), data = H2O2, family=binomial)
fm5 <- glmer(Survival ~ Population + Concentration + Generation + 
               Population:Generation + Concentration:Generation +
               (1|Plate), data = H2O2, family=binomial)

AIC(fm3)
#5211.267
AIC(fm4)
#5552.37
AIC(fm5)
#5213.62
#All higher than fm2. Keep all 2-way

fm6 <- glmer(Survival ~ Population + Concentration + 
               Population:Concentration + Population:Generation + Concentration:Generation +
               (1|Plate), data = H2O2, family=binomial)
fm7 <- glmer(Survival ~ Population + Generation + 
               Population:Concentration + Population:Generation + Concentration:Generation +
               (1|Plate), data = H2O2, family=binomial)
fm8 <- glmer(Survival ~ Concentration + Generation + 
               Population:Concentration + Population:Generation + Concentration:Generation +
               (1|Plate), data = H2O2, family=binomial)

AIC(fm6)
#5022.3
AIC(fm7)
#5022.3
AIC(fm8)
#5022.3

fm9 <- glmer(Survival ~ Population:Concentration + Population:Generation + Concentration:Generation +
               (1|Plate), data = H2O2, family=binomial)
AIC(fm9)
#5022.3

#AIC same as fm2. Drop single factors

fm10 <- glmer(Survival ~ Population:Generation + Concentration:Generation +
                (1|Plate), data = H2O2, family=binomial)

fm11 <- glmer(Survival ~ Population:Concentration + Concentration:Generation +
                (1|Plate), data = H2O2, family=binomial)

fm12 <- glmer(Survival ~ Population:Concentration + Population:Generation +
                (1|Plate), data = H2O2, family=binomial)

AIC(fm10)
#5213.6
AIC(fm11)
#5552.4
AIC(fm12)
#5211.3
#All higher AIC than fm9. fm9 = simplest, best model

summary(fm9)
Anova(fm9)
#chi-sq P:C = 1608.94, df = 35, p < 0.001
#chi-sq P:G = 189.30, df = 4, p < 0.001
#chi-sq C:G = 88.64, df = 6, p < 0.001

#Get stats and check model. Currently using fm10, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm9)
plot(simulationOutput)
testDispersion(simulationOutput)

mmPC <- emmeans(fm9, ~ Population:Concentration)
bonPC <- contrast(mmPC, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonAll, file = "Stats/Survival_stats_PC.txt", sep = "\t", row.names = F)
bonPC <- read.delim("Stats/Survival_stats_PC.txt")
bonPC$Significance <- ifelse(bonPC$p.value <= 0.001, "***", ifelse(bonPC$p.value <= 0.01, "**", ifelse(bonPC$p.value <= 0.05, "*", "ns")))
write.table(bonPC, file = "Stats/Survival_stats_PC.txt", sep = "\t", row.names = F)


mmPG <- emmeans(fm9, ~ Population:Generation)
bonPG <- contrast(mmPG, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonPG, file = "Stats/Survival_stats_PG.txt", sep = "\t", row.names = F)
bonPG <- read.delim("Stats/Survival_stats_PG.txt")
bonPG$Significance <- ifelse(bonPG$p.value <= 0.001, "***", ifelse(bonPG$p.value <= 0.01, "**", ifelse(bonPG$p.value <= 0.05, "*", "ns")))
write.table(bonPG, file = "Stats/Survival_stats_PG.txt", sep = "\t", row.names = F)

mmGC <- emmeans(fm9, ~ Generation:Concentration)
bonGC <- contrast(mmGC, "pairwise", simple = "Generation", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonGC, file = "Stats/Survival_stats_PG.txt", sep = "\t", row.names = F)
bonGC <- read.delim("Stats/Survival_stats_GC.txt")
bonGC$Significance <- ifelse(bonGC$p.value <= 0.001, "***", ifelse(bonGC$p.value <= 0.01, "**", ifelse(bonGC$p.value <= 0.05, "*", "ns")))
write.table(bonGC, file = "Stats/Survival_stats_GC.txt", sep = "\t", row.names = F)

####Survival by replicate####
#####100 gens#####
######R1######
#Generate models from complicated to simple
fm1 <- glm(Survival ~ Population + Concentration + Population:Concentration, data = subset(H2O2_G100, Plate == "R1"), family=binomial)
#make model - fixed effect = Gens, random effect = replicate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glm(Survival ~ Population + Concentration , data = subset(H2O2_G100, Plate == "R1"), family=binomial)

AIC(fm1)
#1047.8
AIC(fm2)
#1151.5
#fm1 has lowest AIC. Interaction kept

fm3 <- glm(Survival ~ Population + Population:Concentration, data = subset(H2O2_G100, Plate == "R1"), family=binomial)
fm4 <- glm(Survival ~ Concentration + Population:Concentration, data = subset(H2O2_G100, Plate == "R1"), family=binomial)
AIC(fm3)
#1047.8
AIC(fm4)
#1047.8
#fm3 & fm4 no change in AIC. both single effects dropped

fm5 <- glm(Survival ~ Population:Concentration, data = subset(H2O2_G100, Plate == "R1"), family=binomial)
AIC(fm5)
#1049.8
#Slight increase in AIC vs. fm1
#Small enough as dropping single factors does not increase AIC
#fm5 is simplest best model

summary(fm5)
Anova(fm5)
#chi-sq interaction = 2016.7, df = 36, p < 0.001

#Get stats and check model. Currently using fm5, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm5)
plot(simulationOutput)
testDispersion(simulationOutput)

mmG100 <- emmeans(fm5, ~ Population:Concentration)
bonG100 <- contrast(mmG100, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG100, file = "Stats/Survival_stats_G100_R1.txt", sep = "\t", row.names = F)
bonG100 <- read.delim("Stats/Survival_stats_G100_R1.txt")
bonG100$Significance <- ifelse(bonG100$p.value <= 0.001, "***", ifelse(bonG100$p.value <= 0.01, "**", ifelse(bonG100$p.value <= 0.05, "*", "ns")))
write.table(bonG100, file = "Stats/Survival_stats_G100_R1.txt", sep = "\t", row.names = F)

######R2######
#Generate models from complicated to simple
fm1 <- glm(Survival ~ Population + Concentration + Population:Concentration, data = subset(H2O2_G100, Plate == "R2"), family=binomial)
#make model - fixed effect = Gens, random effect = replicate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glm(Survival ~ Population + Concentration , data = subset(H2O2_G100, Plate == "R2"), family=binomial)

AIC(fm1)
#639.4
AIC(fm2)
#719.8
#fm1 has lowest AIC. Interaction kept

fm3 <- glm(Survival ~ Population + Population:Concentration, data = subset(H2O2_G100, Plate == "R2"), family=binomial)
fm4 <- glm(Survival ~ Concentration + Population:Concentration, data = subset(H2O2_G100, Plate == "R2"), family=binomial)
AIC(fm3)
#639.4
AIC(fm4)
#639.4
#fm3 & fm4 no change in AIC. both single effects dropped

fm5 <- glm(Survival ~ Population:Concentration, data = subset(H2O2_G100, Plate == "R2"), family=binomial)
AIC(fm5)
#641.4
#Slight increase in AIC vs. fm1
#Small enough as dropping single factors does not increase AIC
#fm5 is simplest best model

summary(fm5)
Anova(fm5)
#chi-sq interaction = 2193.3, df = 36, p < 0.001

#Get stats and check model. Currently using fm5, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm5)
plot(simulationOutput)
testDispersion(simulationOutput)

mmG100 <- emmeans(fm5, ~ Population:Concentration)
bonG100 <- contrast(mmG100, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG100, file = "Stats/Survival_stats_G100_R2.txt", sep = "\t", row.names = F)
bonG100 <- read.delim("Stats/Survival_stats_G100_R2.txt")
bonG100$Significance <- ifelse(bonG100$p.value <= 0.001, "***", ifelse(bonG100$p.value <= 0.01, "**", ifelse(bonG100$p.value <= 0.05, "*", "ns")))
write.table(bonG100, file = "Stats/Survival_stats_G100_R2.txt", sep = "\t", row.names = F)

######R3######
#Generate models from complicated to simple
fm1 <- glm(Survival ~ Population + Concentration + Population:Concentration, data = subset(H2O2_G100, Plate == "R3"), family=binomial)
#make model - fixed effect = Gens, random effect = replicate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glm(Survival ~ Population + Concentration , data = subset(H2O2_G100, Plate == "R3"), family=binomial)

AIC(fm1)
#736.9
AIC(fm2)
#756.6
#fm1 has lowest AIC. Interaction kept

fm3 <- glm(Survival ~ Population + Population:Concentration, data = subset(H2O2_G100, Plate == "R3"), family=binomial)
fm4 <- glm(Survival ~ Concentration + Population:Concentration, data = subset(H2O2_G100, Plate == "R3"), family=binomial)
AIC(fm3)
#736.9
AIC(fm4)
#736.9
#fm3 & fm4 no change in AIC. both single effects dropped

fm5 <- glm(Survival ~ Population:Concentration, data = subset(H2O2_G100, Plate == "R3"), family=binomial)
AIC(fm5)
#738.9
#Slight increase in AIC vs. fm1
#Small enough as dropping single factors does not increase AIC
#fm5 is simplest best model

summary(fm5)
Anova(fm5)
#chi-sq interaction = 2254, df = 36, p < 0.001

#Get stats and check model. Currently using fm5, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm5)
plot(simulationOutput)
testDispersion(simulationOutput)

mmG100 <- emmeans(fm5, ~ Population:Concentration)
bonG100 <- contrast(mmG100, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG100, file = "Stats/Survival_stats_G100_R3.txt", sep = "\t", row.names = F)
bonG100 <- read.delim("Stats/Survival_stats_G100_R3.txt")
bonG100$Significance <- ifelse(bonG100$p.value <= 0.001, "***", ifelse(bonG100$p.value <= 0.01, "**", ifelse(bonG100$p.value <= 0.05, "*", "ns")))
write.table(bonG100, file = "Stats/Survival_stats_G100_R3.txt", sep = "\t", row.names = F)

######R4######
#Generate models from complicated to simple
fm1 <- glm(Survival ~ Population + Concentration + Population:Concentration, data = subset(H2O2_G100, Plate == "R4"), family=binomial)
#make model - fixed effect = Gens, random effect = replicate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glm(Survival ~ Population + Concentration , data = subset(H2O2_G100, Plate == "R4"), family=binomial)

AIC(fm1)
#612.96
AIC(fm2)
#585.0
#fm2 has lowest AIC. Interaction dropped

fm3 <- glm(Survival ~ Population , data = subset(H2O2_G100, Plate == "R4"), family=binomial)
fm4 <- glm(Survival ~ Concentration, data = subset(H2O2_G100, Plate == "R4"), family=binomial)
AIC(fm3)
#2621.3
AIC(fm4)
#1192.1
#fm3 & fm4 increase AIC. both single effects kept
#fm2 best model

summary(fm2)
Anova(fm2)
#chi-sq pop = 613.1, df = 3, p < 0.001
#chi-sq conc = 2052.2, df = 8, p < 0.001

#Get stats and check model. Currently using fm5, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm2)
plot(simulationOutput)
testDispersion(simulationOutput)

mmG100 <- emmeans(fm2, ~ Population)
bonG100 <- contrast(mmG100, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG100, file = "Stats/Survival_stats_G100_R4.txt", sep = "\t", row.names = F)
bonG100 <- read.delim("Stats/Survival_stats_G100_R4.txt")
bonG100$Significance <- ifelse(bonG100$p.value <= 0.001, "***", ifelse(bonG100$p.value <= 0.01, "**", ifelse(bonG100$p.value <= 0.05, "*", "ns")))
write.table(bonG100, file = "Stats/Survival_stats_G100_R4.txt", sep = "\t", row.names = F)

#####1000 gens#####
######R2######
#Generate models from complicated to simple
fm1 <- glm(Survival ~ Population + Concentration + Population:Concentration, data = subset(H2O2_G1K, Plate == "R2"), family=binomial)
#make model - fixed effect = Gens, random effect = replicate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glm(Survival ~ Population + Concentration , data = subset(H2O2_G1K, Plate == "R2"), family=binomial)

AIC(fm1)
#314.3
AIC(fm2)
#280.1
#fm2 has lowest AIC. Interaction dropped

fm3 <- glm(Survival ~ Population, data = subset(H2O2_G1K, Plate == "R2"), family=binomial)
fm4 <- glm(Survival ~ Concentration, data = subset(H2O2_G1K, Plate == "R2"), family=binomial)
AIC(fm3)
#1925.2
AIC(fm4)
#591.4
#fm3 & fm4 increase AIC. both single effects kept
#fm2 best, simplest model

summary(fm2)
Anova(fm2)
#chi-sq pop = 317.2, df = 3, p < 0.001
#chi-sq conc = 1657.1, df = 6, p < 0.001

#Get stats and check model. Currently using fm5, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm2)
plot(simulationOutput)
testDispersion(simulationOutput)

mmG1K <- emmeans(fm2, ~ Population)
bonG1K <- contrast(mmG1K, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG1K, file = "Stats/Survival_stats_G1K_R2.txt", sep = "\t", row.names = F)
bonG1K_R2 <- read.delim("Stats/Survival_stats_G1K_R2.txt")
bonG1K$Significance <- ifelse(bonG1K$p.value <= 0.001, "***", ifelse(bonG1K$p.value <= 0.01, "**", ifelse(bonG1K$p.value <= 0.05, "*", "ns")))
write.table(bonG1K, file = "Stats/Survival_stats_G1K_R2.txt", sep = "\t", row.names = F)

######R3######
#Generate models from complicated to simple
fm1 <- glm(Survival ~ Population + Concentration + Population:Concentration, data = subset(H2O2_G1K, Plate == "R3"), family=binomial)
#make model - fixed effect = Gens, random effect = replicate, response variable = yield in novel environments
simulationOutput <- simulateResiduals(fittedModel = fm1)
plot(simulationOutput)
testDispersion(simulationOutput)
fm2 <- glm(Survival ~ Population + Concentration , data = subset(H2O2_G1K, Plate == "R3"), family=binomial)

AIC(fm1)
#386.4
AIC(fm2)
#371.3
#fm2 has lowest AIC. Interaction dropped

fm3 <- glm(Survival ~ Population, data = subset(H2O2_G1K, Plate == "R3"), family=binomial)
fm4 <- glm(Survival ~ Concentration, data = subset(H2O2_G1K, Plate == "R3"), family=binomial)
AIC(fm3)
#2059.0
AIC(fm4)
#793.0
#fm3 & fm4 increase AIC. both single effects kept
#fm2 is simplest best model

summary(fm2)
Anova(fm2)
#chi-sq pop = 427.7, df = 3, p < 0.001
#chi-sq conc = 1699.7, df = 6, p < 0.001

#Get stats and check model. Currently using fm5, but change to lowest AIC model
simulationOutput <- simulateResiduals(fittedModel = fm2)
plot(simulationOutput)
testDispersion(simulationOutput)

mmG1K <- emmeans(fm2, ~ Population)
bonG1K <- contrast(mmG1K, "pairwise", simple = "Population", combine = TRUE, adjust = "bonferroni")
#post-hoc test for comparisons of Survival between various generations
write.table(bonG1K, file = "Stats/Survival_stats_G1K_R3.txt", sep = "\t", row.names = F)
bonG1K_R3 <- read.delim("Stats/Survival_stats_G1K_R3.txt")
bonG1K$Significance <- ifelse(bonG1K$p.value <= 0.001, "***", ifelse(bonG1K$p.value <= 0.01, "**", ifelse(bonG1K$p.value <= 0.05, "*", "ns")))
write.table(bonG1K, file = "Stats/Survival_stats_G1K_R3.txt", sep = "\t", row.names = F)

####Hybrid deterioration Figures####
#Updated 2022-09-05
#Plots by ggplot2 with gridExtra
#R version 4.2.1 (2022-06-23)
#Packages used - dplyr, ggplot2 and gridExtra

#Get data
setwd("PhD/Experiments/Hybrid Novel Environment/Survival")
H2O2_G1000 <- read.delim("1000 gens/H2O2Growth.txt")
str(H2O2_G1000)
H2O2_G1000$Population <- factor(H2O2_G1000$Population, levels = c("P1", "P2", "F1", "F2"), ordered = TRUE)
H2O2_G1000$Plate <- as.factor(H2O2_G1000$Plate)
H2O2_G1000$Well <- as.factor(H2O2_G1000$Well)
H2O2_G1000 <- subset(H2O2_G1000, Plate != "R4" & Plate != "R1")
H2O2_G1000$Generation <- 1000
H2O2_G1000 <- droplevels(H2O2_G1000)
str(H2O2_G1000)

H2O2_G100 <- read.delim("100 gens/H2O2Growth.txt")
str(H2O2_G100)
H2O2_G100$Population <- factor(H2O2_G100$Population, levels = c("P1", "P2", "F1", "F2"), ordered = TRUE)
H2O2_G100$Plate <- as.factor(H2O2_G100$Plate)
H2O2_G100$Well <- as.factor(H2O2_G100$Well)
H2O2_G100$Generation <- 100
H2O2_G100 <- droplevels(H2O2_G100)
str(H2O2_G100)

H2O2 <- rbind(H2O2_G100, H2O2_G1000)
write.table(H2O2, file = "Growth.txt", sep = "\t", row.names = F)

H2O2Survivors_G1000 <- read.delim("1000 gens/H2O2Survivors.txt")
str(H2O2Survivors_G1000)
H2O2Survivors_G1000$Population <- factor(H2O2Survivors_G1000$Population, levels = c("P1", "P2", "F1", "F2"), ordered = TRUE)
H2O2Survivors_G1000$Concentration <- factor(H2O2Survivors_G1000$Concentration, levels = c("0.02M", "0.04M", "0.06M", "0.08M", "0.1M", "0.12M", "0.14M"), ordered = TRUE)
H2O2Survivors_G1000 <- subset(H2O2Survivors_G1000, Replicate != "R4" & Replicate != "R1")
H2O2Survivors_G1000$Generation <- 1000
H2O2Survivors_G1000 <- droplevels(H2O2Survivors_G1000)
str(H2O2Survivors_G1000)

H2O2Survivors_G100 <- read.delim("100 gens/H2O2Survivors.txt")
str(H2O2Survivors_G100)
H2O2Survivors_G100$Population <- factor(H2O2Survivors_G100$Population, levels = c("P1", "P2", "F1", "F2"), ordered = TRUE)
H2O2Survivors_G100$Concentration <- factor(H2O2Survivors_G100$Concentration, levels = c("0.02M", "0.04M", "0.06M", "0.08M", "0.1M", "0.12M", "0.14M", "0.16M", "0.18M"), ordered = TRUE)
H2O2Survivors_G100$Generation <- 100
H2O2Survivors_G100 <- droplevels(H2O2Survivors_G100)
str(H2O2Survivors_G100)

H2O2Survivors <- rbind(H2O2Survivors_G100, H2O2Survivors_G1000)
write.table(H2O2Survivors, file = "SurvivorsByReplicate.txt", sep = "\t", row.names = F)


H2O2PopSurvivors_G100 <- read.delim("100 gens/H2O2SurvivalPopulation.txt")
str(H2O2PopSurvivors_G100)
H2O2PopSurvivors_G100$Population <- factor(H2O2PopSurvivors_G100$Population, levels = c("P1", "P2", "F1", "F2"), ordered = TRUE)
H2O2PopSurvivors_G100$Generation <- 100
str(H2O2PopSurvivors_G100)

H2O2PopSurvivors_G1000 <- read.delim("1000 gens/H2O2SurvivalPopulation.txt")
str(H2O2PopSurvivors_G1000)
H2O2PopSurvivors_G1000$Population <- factor(H2O2PopSurvivors_G1000$Population, levels = c("P1", "P2", "F1", "F2"), ordered = TRUE)
H2O2PopSurvivors_G1000$Generation <- 1000
str(H2O2PopSurvivors_G1000)

H2O2PopSurvivors <- rbind(H2O2PopSurvivors_G100, H2O2PopSurvivors_G1000)
write.table(H2O2PopSurvivors, file = "SurvivorsByPopulation.txt", sep = "\t", row.names = F)

rm(list = ls())

####Survivors####

library(dplyr)
H2O2 <- read.delim("Growth.txt")
H2O2$Survival <- as.factor(H2O2$Survival)
H2O2$Plate <- as.factor(H2O2$Plate)
H2O2$Concentration <- as.factor(H2O2$Concentration)
H2O2$Generation <- as.factor(H2O2$Generation)
H2O2 <- droplevels(H2O2)
str(H2O2)
H2O2_Survivors <- H2O2 %>% group_by(Population, Plate, Concentration, Generation, Survival, .drop = FALSE) %>% count()
H2O2_Survivors <- subset(H2O2_Survivors, Survival == "1")
H2O2_Survivors$Proportion <- H2O2_Survivors$n/60

H2O2_Survivors_G1K <- subset(H2O2_Survivors, Generation == "1000")
H2O2_Survivors_G1K <- subset(H2O2_Survivors_G1K, Plate == "R2"|Plate == "R3")
H2O2_Survivors_G1K <- subset(H2O2_Survivors_G1K, Concentration != "0.16" & Concentration != "0.18")
H2O2_Survivors_G1 <- subset(H2O2_Survivors, Generation == "100")

H2O2_Survivors <- rbind(H2O2_Survivors_G1, H2O2_Survivors_G1K)
write.table(H2O2_Survivors, file = "SurvivorsByReplicate.txt", sep = "\t", row.names = F)

H2O2_Survivors_all <- H2O2 %>% group_by(Population, Concentration, Generation, Survival, .drop = FALSE) %>% count()
H2O2_Survivors_all <- subset(H2O2_Survivors_all, Survival == "1")
H2O2_Survivors_all$Proportion <- H2O2_Survivors_all$n/60

H2O2_Survivors_all_G1K <- subset(H2O2_Survivors_all, Generation == "1000")
H2O2_Survivors_all_G1K <- subset(H2O2_Survivors_all_G1K, Concentration != "0.16" & Concentration != "0.18")
H2O2_Survivors_all_G1K$Proportion <- H2O2_Survivors_all_G1K$n/120
H2O2_Survivors_all_G1 <- subset(H2O2_Survivors_all, Generation == "100")
H2O2_Survivors_all_G1$Proportion <- H2O2_Survivors_all_G1$n/240
H2O2_Survivors_all <- rbind(H2O2_Survivors_all_G1K, H2O2_Survivors_all_G1)
write.table(H2O2_Survivors_all, file = "SurvivorsByPopulation.txt", sep = "\t", row.names = F)

H2O2Survivors <- read.delim("SurvivorsByReplicate.txt")
H2O2Survivors$Population <- factor(H2O2Survivors$Population, levels = c("P1", "P2", "F1", "F2"), ordered = TRUE)
str(H2O2Survivors)

days <- c("2", "4", "6", "8", "10", "12", "14", "16", "18")

#Line by survival percentage
library(ggplot2)
library(gridExtra)

#AllRs
p1 <- ggplot(data = subset(H2O2Survivors, Generation == "100"), aes(x = Concentration, y = Proportion, group = Plate, colour = Population, fill = Population)) +
  geom_line(size=0.75, aes(group = Plate)) + labs(x = "H2O2 concentration (M)", y = "Survival %", title = "Survival under deterioraring conditions") +
  geom_point(shape = 21, colour = "black", size = 3, aes(group = Plate)) +
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  scale_colour_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) + scale_fill_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=10, face="bold"),axis.text=element_text(size=12), axis.ticks.length = unit(0.25, "cm"), 
                     axis.title=element_text(size=12,face="bold"))

png("Plots/Survival_H202_G100.png", res = 500, units = "in", height = 5, width = 10)
p1
dev.off()

p2 <- ggplot(data = subset(H2O2Survivors, Generation == "100"), aes(x = Concentration, y = Proportion, colour = Population, fill = Population, group = Population)) +
  geom_line(size=0.75, linetype = "dashed") + labs(x = "H2O2 concentration (M)", y = "Survival Proportion", title = "A) 100 Generations") +
  geom_point(shape = 21, colour = "black", size = 3) +
  facet_wrap(.~Plate, scales = "free", ncol = 1) +
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  scale_x_continuous(breaks=seq(0.02, 0.18, 0.02), limits = c(0.02, 0.18), sec.axis = dup_axis(name = "Days", labels = days)) +
  scale_colour_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) + scale_fill_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), title = element_text(size=22,face="bold"),  
                     axis.text=element_text(size=20), axis.ticks.length = unit(0.25, "cm"), axis.title=element_text(size=22,face="bold"), 
                     strip.text.x = element_text(size = 22, face ="bold"), strip.background = element_rect(colour="black", fill="white", size=1, linetype="blank"), strip.placement = "outside",
                     legend.position = c(1, 1), legend.justification = c(1, 1),legend.title=element_text(size=20), legend.text=element_text(size=20), legend.background = element_rect(fill = NA))


png("Plots/Survival_H202_G100_ByReplicate.png", res = 500, units = "in", height = 20, width = 7)
p2
dev.off()


p3 <- ggplot(data = subset(H2O2Survivors, Generation == "1000"), aes(x = Concentration, y = Proportion, group = Plate, colour = Population, fill = Population)) +
  geom_line(size=0.75) + labs(x = "H2O2 concentration (M)", y = "Survival %", title = "Survival under deterioraring conditions") +
  geom_point(shape = 21, colour = "black", size = 3) +
  scale_y_continuous(breaks=seq(0, 100, 20)) +
  scale_colour_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) + scale_fill_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=10, face="bold"),axis.text=element_text(size=12), axis.ticks.length = unit(0.25, "cm"), 
                     axis.title=element_text(size=12,face="bold"))

png("Plots/Survival_H202_G1000.png", res = 500, units = "in", height = 5, width = 10)
p3
dev.off()

p4 <- ggplot(data = subset(H2O2Survivors, Generation == "1000"), aes(x = Concentration, y = Proportion, colour = Population, fill = Population, group = Population)) +
  geom_line(size=0.75) + labs(x = "H2O2 concentration (M)", y = "Survival Proportion", title = "B) 1000 Generations") +
  geom_point(shape = 21, colour = "black", size = 3) +
  facet_wrap(.~Plate, scales = "free", ncol = 1) +
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  scale_x_continuous(breaks=seq(0.02, 0.18, 0.02), limits = c(0.02, 0.18), sec.axis = dup_axis(name = "Days", labels = days)) +
  scale_colour_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) + scale_fill_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), title = element_text(size=22,face="bold"),  
                     axis.text=element_text(size=20), axis.ticks.length = unit(0.25, "cm"), axis.title=element_text(size=22,face="bold"), 
                     strip.text.x = element_text(size = 22, face ="bold"), strip.background = element_rect(colour="black", fill="white", size=1, linetype="blank"), strip.placement = "outside",
                     legend.position = c(1, 1), legend.justification = c(1, 1),legend.title=element_text(size=20), legend.text=element_text(size=20), legend.background = element_rect(fill = NA))

png("Plots/Survival_H202_G1000_ByReplicate.png", res = 500, units = "in", height = 10, width = 7)
p4
dev.off()

png("Plots/Figure S1.png", res = 1000, units = "in", height = 20, width = 14)
grid.arrange(p2, p4, ncol = 2)
dev.off()

p5 <- ggplot(data = H2O2Survivors, aes(x = Concentration, y = Proportion, colour = Population, fill = Population, group = Population)) +
  geom_line(data = subset(H2O2Survivors, Generation == "100"), size=1, linetype = "dashed") + geom_line(data = subset(H2O2Survivors, Generation == "1000"), size=1) +
  labs(x = "H2O2 concentration (M)", y = "Survival proportion", title = "1000 Generations") +
  geom_point(data = subset(H2O2Survivors, Generation == "100"), shape = 21, colour = "black", size = 3.5) + geom_point(data = subset(H2O2Survivors, Generation == "1000"), shape = 21, colour = "black", size = 3.5) +
  facet_wrap(.~Plate, scales = "free", ncol = 2) +
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  scale_x_continuous(breaks=seq(0.02, 0.18, 0.02), limits = c(0.02, 0.18)) +
  scale_colour_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) + scale_fill_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), title = element_text(size=18, face="bold"),axis.text=element_text(size=14), axis.ticks.length = unit(0.25, "cm"), 
                     axis.title=element_text(size=18,face="bold"), strip.text.x = element_text(size = 16, face ="bold"), strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid"),
                     legend.title=element_text(size=16), legend.text=element_text(size=14))

png("Plots/Survival_H202_AllGens_ByReplicate.png", res = 500, units = "in", height = 10, width = 12)
p5
dev.off()


H2O2PopSurvivors <- read.delim("SurvivorsByPopulation.txt")
H2O2PopSurvivors$Population <- factor(H2O2PopSurvivors$Population, levels = c("P1", "P2", "F1", "F2"), ordered = TRUE)
str(H2O2PopSurvivors)
days <- c("2", "4", "6", "8", "10", "12", "14", "16", "18")


#Survival by pop
library(ggplot2)
library(gridExtra)

p6 <- ggplot(data = subset(H2O2PopSurvivors, Generation == "100"), aes(x = Concentration, y = Proportion, colour = Population, fill = Population, group = Population)) +
  geom_line(size=0.75, linetype = "dashed") + labs(x = "H2O2 concentration (M)", y = "Survival Proportion", title = "A) 100 Generations") +
  geom_point(shape = 21, colour = "black", size = 3) +
  geom_hline(yintercept=0, linetype="dashed", size = 0.5) +
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  scale_x_continuous(breaks=seq(0.02, 0.18, 0.02), limits = c(0.02, 0.18), sec.axis = dup_axis(name = "Days", labels = days)) +
  scale_colour_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) + scale_fill_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), title = element_text(size=14,face="bold"),  
                     axis.text=element_text(size=12), axis.ticks.length = unit(0.25, "cm"), axis.title=element_text(size=14,face="bold"), 
                     legend.position = c(1, 1), legend.justification = c(1, 1),legend.title=element_text(size=12), legend.text=element_text(size=12), legend.background = element_rect(fill = NA))


png("Plots/Survival_Pops_H202G100.png", res = 500, units = "in", height = 5, width = 10)
p6
dev.off()

p7 <- ggplot(data = subset(H2O2PopSurvivors, Generation == "1000"), aes(x = Concentration, y = Proportion, colour = Population, fill = Population, group = Population)) +
  geom_line(size=0.75) + labs(x = "H2O2 concentration (M)", y = "Survival Proportion", title = "B) 1000 Generations") +
  geom_point(shape = 21, colour = "black", size = 3) +
  geom_hline(yintercept=0, linetype="dashed", size = 0.5) +
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  scale_x_continuous(breaks=seq(0.02, 0.18, 0.02), limits = c(0.02, 0.18), sec.axis = dup_axis(name = "Days", labels = days)) +
  scale_colour_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) + scale_fill_manual(values=c("#D3D3D3", "#A9A9A9", "#FFA500", "#AC1016")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), title = element_text(size=18,face="bold"),  
                     axis.text=element_text(size=16), axis.ticks.length = unit(0.25, "cm"), axis.title=element_text(size=18,face="bold"), 
                     legend.position = "none", legend.justification = c(1, 1),legend.title=element_text(size=16), legend.text=element_text(size=16), legend.background = element_rect(fill = NA))


png("Plots/Survival_Pops_H202G1000.png", res = 500, units = "in", height = 5, width = 10)
p7
dev.off()

png("Plots/Figure 2.png", res = 1000, units = "in", height = 10, width = 7)
grid.arrange(p6, p7, nrow = 2)
dev.off()
