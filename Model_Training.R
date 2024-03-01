#Model Training

library(mlr3)
library(mlr3proba)
library(survival)
library(mlr3extralearners)
library(mlr3learners)
library(mlr3verse)
library(mlr3viz)
library(mlr3tuning)
library(ggplot2)

#Coxph Cli
Data_ML_TV <- as.data.frame(cbind(Label_PM$FUorRInterval[TrainingIndex], Label_PM$Recurrence[TrainingIndex], ClinicalFeatures[TrainingIndex,]))
colnames(Data_ML_TV)[1] <- "time"
colnames(Data_ML_TV)[2] <- "status"
task_TV <- as_task_surv(Data_ML_TV, id = "right_censored",
                       time = "time", event = "status", type = "right")
print(task_TV)
Data_ML <- as.data.frame(cbind(Label_PM$FUorRInterval, Label_PM$Recurrence, ClinicalFeatures))
colnames(Data_ML)[1] <- "time"
colnames(Data_ML)[2] <- "status"
task <- as_task_surv(Data_ML, id = "right_censored",
                    time = "time", event = "status", type = "right")
print(task)
set.seed(0728)
learner <- lrn("surv.coxph")
learner$param_set
learner
rr <- resample(task_TV, learner, rsmp("cv", folds = 3))
rr$aggregate(msr("surv.cindex"))
#0.7248366
learner$train(task,row_ids = TrainingIndex)
learner$predict(task,row_ids = TrainingIndex)$score(msr("surv.cindex"))
#0.7961538
learner$predict(task,row_ids = TestingIndex)$score(msr("surv.cindex"))
#0.7096774

#Rsf Cli
set.seed(0203)
learner_RF_Best <- lrn("surv.rfsrc", ntree = 1000, mtry = 3, nodesize = 8, importance = "permute")
learner_RF_Best$train(task,row_ids = TrainingIndex)
print(learner_RF_Best$predict(task,row_ids = TrainingIndex)$score(msr("surv.cindex")))
#0.8730769
print(learner_RF_Best$predict(task,row_ids = TestingIndex)$score(msr("surv.cindex")))
#0.7419355

#Coxph Prot
Data_ML <- as.data.frame(cbind(Label_PM$FUorRInterval, Label_PM$Recurrence, Ratio_Matrix_R1_impseqrob_combat_PM[,Protein_Index_R2_retain]))
colnames(Data_ML)[1] <- "time"
colnames(Data_ML)[2] <- "status"
task <- as_task_surv(Data_ML, id = "right_censored",
                     time = "time", event = "status", type = "right")
print(task)
learner_LL <- lrn("surv.glmnet", alpha = 1, lambda = 0.1116327)
learner_LL$train(task,row_ids = TrainingIndex)
print(learner_LL$predict(task, row_ids = TrainingIndex)$score(msr("surv.cindex")))
#0.9576923
print(learner_LL$predict(task, row_ids = TestingIndex)$score(msr("surv.cindex")))
#0.8172043
learner_LL$selected_features()

#Rsf Prot
Data_ML_RF <- as.data.frame(cbind(Label_PM$FUorRInterval, Label_PM$Recurrence, Ratio_Matrix_R1_impseqrob_combat_PM[,Protein_Importance_RF]))
colnames(Data_ML_RF)[1] <- "time"
colnames(Data_ML_RF)[2] <- "status"
task_RF <- as_task_surv(Data_ML_RF, id = "right_censored",
                       time = "time", event = "status", type = "right")
print(task_RF)
Data_ML_RF_TV <- as.data.frame(cbind(Label_PM$FUorRInterval[TrainingIndex], Label_PM$Recurrence[TrainingIndex], Ratio_Matrix_R1_impseqrob_combat_PM[TrainingIndex, Protein_Importance_RF]))
colnames(Data_ML_RF_TV)[1] <- "time"
colnames(Data_ML_RF_TV)[2] <- "status"
task_RF_TV <- as_task_surv(Data_ML_RF_TV, id = "right_censored",
                       time = "time", event = "status", type = "right")
print(task_RF_TV)
set.seed(5)
learner_RF_TV <- lrn("surv.rfsrc", ntree = 1000, mtry = 4, nodesize = 6)
rr <- resample(task_RF_TV, learner_RF_TV, rsmp("cv", folds = 3))
rr$aggregate(msr("surv.cindex"))
#0.9686486
learner_RF_Best <- lrn("surv.rfsrc", ntree = 1000, mtry = 4, nodesize = 6, importance = "permute")
learner_RF_Best$train(task_RF, row_ids = TrainingIndex)
print(learner_RF_Best$predict(task_RF, row_ids = TrainingIndex)$score(msr("surv.cindex")))
#0.9961538
print(learner_RF_Best$predict(task_RF, row_ids = TestingIndex)$score(msr("surv.cindex")))
#0.8494624

#Rsf All
Data_ML_RF2 <- as.data.frame(cbind(Label_PM$FUorRInterval, Label_PM$Recurrence, Ratio_Matrix_R1_impseqrob_combat_PM[,Protein_Importance_RF2]))
colnames(Data_ML_RF2)[1] <- "time"
colnames(Data_ML_RF2)[2] <- "status"
task_RF2 <- as_task_surv(Data_ML_RF2, id = "right_censored",
                       time = "time", event = "status", type = "right")
print(task_RF2)
Data_ML_RF_TV2 <- as.data.frame(cbind(Label_PM$FUorRInterval[TrainingIndex], Label_PM$Recurrence[TrainingIndex], Ratio_Matrix_R1_impseqrob_combat_PM[TrainingIndex, Protein_Importance_RF2]))
colnames(Data_ML_RF_TV2)[1] <- "time"
colnames(Data_ML_RF_TV2)[2] <- "status"
task_RF_TV2 <- as_task_surv(Data_ML_RF_TV2, id = "right_censored",
                          time = "time", event = "status", type = "right")
print(task_RF_TV2)
set.seed(5)
learner_RF_TV2 <- lrn("surv.rfsrc", ntree = 1000, mtry = 5, nodesize = 7)
rr2 <- resample(task_RF_TV2, learner_RF_TV2, rsmp("cv", folds = 3))
rr2$aggregate(msr("surv.cindex"))
#0.9686486
learner_RF_Best2 <- lrn("surv.rfsrc", ntree = 1000, mtry = 5, nodesize = 7)
learner_RF_Best2$train(task_RF2, row_ids = TrainingIndex)
print(learner_RF_Best2$predict(task_RF2, row_ids = TrainingIndex)$score(msr("surv.cindex")))
#0.9961538
print(learner_RF_Best2$predict(task_RF2, row_ids = TestingIndex)$score(msr("surv.cindex")))
#0.827957
