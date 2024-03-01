#Cox Survival Analysis

library(mlr3)
library(mlr3proba)
library(survival)
library(mlr3extralearners)
library(mlr3learners)
library(mlr3verse)
library(mlr3viz)
library(mlr3tuning)
library(ggplot2)
library(survival)
library(survminer)

#PM DATA
ClinicalFeatures <- Label_PM[,c(5:15)]
str(ClinicalFeatures)
ClinicalFeatures$Gender <- as.factor(ClinicalFeatures$Gender)
ClinicalFeatures$HT <- as.factor(ClinicalFeatures$HT)
ClinicalFeatures$Mulifocality <- as.factor(ClinicalFeatures$Mulifocality)
ClinicalFeatures$ETE <- as.factor(ClinicalFeatures$ETE)
ClinicalFeatures$SurgeryType <- as.factor(ClinicalFeatures$SurgeryType)
str(ClinicalFeatures)

#All Data
Data_ML <- as.data.frame(cbind(Label_PM$FUorRInterval, Label_PM$Recurrence, ClinicalFeatures))
colnames(Data_ML)[1] <- "time"
colnames(Data_ML)[2] <- "status"

#Univariate Cox
library(survival)
(cox1 <- coxph(Surv(time, status) ~ Gender, Data_ML))#p 0.347
summary(cox1)#GenderMale HR 0.4819 95%CI 0.1051-2.208
(cox2 <- coxph(Surv(time, status) ~ MaxNoduleSize, Data_ML))#p 0.997
summary(cox2)#MaxNoduleSize HR 0.9999 95%CI 0.9557-1.046
(cox3 <- coxph(Surv(time, status) ~ HT, Data_ML))#p 0.571
summary(cox3)#HT1 HR 1.394 95%CI 0.4413-4.404
(cox4 <- coxph(Surv(time, status) ~ Mulifocality, Data_ML))#p 0.126
summary(cox4)#Mulifocality1 HR 2.457 95%CI 0.7778-7.761
(cox5 <- coxph(Surv(time, status) ~ TLNR, Data_ML))#p 0.751
summary(cox5)#TLNR HR 1.313 95%CI 0.2443-7.051
(cox6 <- coxph(Surv(time, status) ~ ETE, Data_ML))#p 0.794
summary(cox6)#ETE1 HR 0.8164 95%CI 0.1785-3.735
(cox7 <- coxph(Surv(time, status) ~ LLNR, Data_ML))#p 0.128
summary(cox7)#LLNR HR 3.152 95%CI 0.7175-13.84
(cox8 <- coxph(Surv(time, status) ~ SurgeryType, Data_ML))#p 0.718
summary(cox8)#SurgeryTypeTotal_thyroidectomy HR 1.232 95%CI 0.397-3.824
(cox1 <- coxph(Surv(time, status) ~ Age, Data_ML))#p 0.0174
summary(cox1)#Age HR 0.7928 95%CI 0.6547-0.96
(cox2 <- coxph(Surv(time, status) ~ TLNN, Data_ML))#p 0.0225
summary(cox2)#TLNN HR 1.076 95%CI 1.01-1.146
(cox3 <- coxph(Surv(time, status) ~ LLNN, Data_ML))#p 0.0111
summary(cox3)#LLNN HR 1.101 95%CI 1.022-1.185

#K-M Curve
library(survival)
library(survminer)
summary(ClinicalFeatures$Age)
AgeGroup <- 1:85
AgeGroup[which(ClinicalFeatures$Age < 16)] <- 0
AgeGroup[which(ClinicalFeatures$Age >= 16)] <- 1
KMdata <- data.frame(Time = Label_PM$FUorRInterval, Status = (as.numeric(Label_PM$Recurrence)), AgeGroup = AgeGroup)
fit <- survfit(Surv(Time, Status) ~ AgeGroup, data = KMdata)
ggsurvplot(fit, data = KMdata, palette = c("#BC3C29FF", "#0072B5FF"), pval = TRUE, pval.method = T, linetype = 1, risk.table = T, size = 2, xlab = "Month", ylab = "Survival probability", legend = "right", legend.title = "Agegroup", pval.size = 7, censor.shape = 124, censor.size = 5, font.tickslab = c(12, "plain", "black"), font.y = c(15, "bold", "black"), font.2 = c(15, "bold", "black"), font.legend = c(12, "bold", "black"))
#p 0.02
rm(KMdata)
rm(fit)
rm(AgeGroup)

summary(ClinicalFeatures$TLNN)
library(survival)
library(survminer)
TLNNHappen <- 1:85
TLNNHappen[which(ClinicalFeatures$TLNN < 5)] <- 0
TLNNHappen[which(ClinicalFeatures$TLNN >= 5)] <- 1
TLNNHappen <- as.factor(TLNNHappen)
KMdata <- data.frame(Time = Label_PM$FUorRInterval, Status = (as.numeric(Label_PM$Recurrence)), TLNN = TLNNHappen)
fit <- survfit(Surv(Time, Status) ~ TLNN, data = KMdata)
ggsurvplot(fit, data = KMdata, palette = c("#BC3C29FF", "#0072B5FF"), pval = TRUE, pval.method = T, linetype = 1, size = 2, risk.table = T, xlab = "Month", ylab = "Survival probability", legend = "right", legend.title = "TLNN", pval.size = 7, censor.shape = 124, censor.size = 5, font.tickslab = c(12, "plain", "black"), font.y = c(15, "bold", "black"), font.2 = c(15, "bold", "black"), font.legend = c(12, "bold", "black"))
#p 0.004
rm(KMdata)
rm(fit)
rm(TLNNHappen)

summary(ClinicalFeatures$LLNN)
library(survival)
library(survminer)
LLNNHappen <- 1:85
LLNNHappen[which(ClinicalFeatures$LLNN < 1)] <- 0
LLNNHappen[which(ClinicalFeatures$LLNN >= 1)] <- 1
LLNNHappen <- as.factor(LLNNHappen)
KMdata <- data.frame(Time = Label_PM$FUorRInterval, Status = (as.numeric(Label_PM$Recurrence)), LLNN = LLNNHappen)
fit <- survfit(Surv(Time, Status) ~ LLNN, data = KMdata)
ggsurvplot(fit, data = KMdata, palette = c("#BC3C29FF", "#0072B5FF"), pval = TRUE, pval.method = T, linetype = 1, size = 2, risk.table = TRUE, xlab = "Month", ylab = "Survival probability", legend = "right", legend.title = "LLNN", pval.size = 7, censor.shape = 124, censor.size = 5, font.tickslab = c(12, "plain", "black"), font.y = c(15, "bold", "black"), font.2 = c(15, "bold", "black"), font.legend = c(12, "bold", "black"))
#p 0.023
rm(KMdata)
rm(fit)
rm(LLNNHappen)

#Multivariate Cox
res.cox <- coxph(Surv(time, status) ~ Gender + MaxNoduleSize + SurgeryType + Mulifocality + HT + LLNN + LLNR + TLNN + TLNR + ETE + Age, Data_ML)
res.cox
summary(res.cox)
ggforest(res.cox, data = Data_ML,
         main = "Hazard ratio",
         cpositions = c(0.05, 0.15, 0.35),
         fontsize = 1.5,
         noDigits = 3)
rm(Data_ML)
rm(res.cox)

AgeGroup16 <- 1:85
AgeGroup16[which(ClinicalFeatures$Age < 16)] <- 0
AgeGroup16[which(ClinicalFeatures$Age >= 16)] <- 1
AgeGroup16 <- as.factor(AgeGroup16)
(cox4 <- coxph(Surv(Label_PM$FUorRInterval, Label_PM$Recurrence) ~ AgeGroup16))
#p 0.0302
summary(cox4)#AgeGroup161 HR 0.2645 95%CI 0.07944-0.8804

#Multivariate Cox
ClinicalFeatures$AgeGroup <- AgeGroup16
ClinicalFeatures <- ClinicalFeatures[,-6]
str(ClinicalFeatures)
Data_ML <- as.data.frame(cbind(Label_PM$FUorRInterval, Label_PM$Recurrence, ClinicalFeatures))
colnames(Data_ML)[1] <- "time"
colnames(Data_ML)[2] <- "status"
res.cox <- coxph(Surv(time, status) ~ ., Data_ML)
res.cox
summary(res.cox)
ggforest(res.cox, data = Data_ML,
         main = "Hazard ratio",
         cpositions = c(0.05, 0.15, 0.35),
         fontsize = 1.5,
         noDigits = 3)
rm(Data_ML)
rm(res.cox)
rm(AgeGroup16)
