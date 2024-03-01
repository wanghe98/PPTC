#Hyperparameter Tuning & Feature Selection

library(mlr3)
library(mlr3proba)
library(survival)
library(mlr3extralearners)
library(mlr3learners)
library(mlr3verse)
library(mlr3viz)
library(mlr3tuning)
library(ggplot2)

#Rsf
Data_ML_TV <- as.data.frame(cbind(Label_PM$FUorRInterval[TrainingIndex], Label_PM$Recurrence[TrainingIndex], ClinicalFeatures[TrainingIndex,]))
colnames(Data_ML_TV)[1] <- "time"
colnames(Data_ML_TV)[2] <- "status"
task_TV <- as_task_surv(Data_ML_TV, id = "right_censored",
                        time = "time", event = "status", type = "right")
print(task_TV)
set.seed(41)
learner_RF <- lrn("surv.rfsrc", ntree = 1000, mtry = 3)
search_space <- ps(
  nodesize = p_int(lower = 1, upper = 15)
)
cv_inner <- rsmp("cv", folds = 3)
measure <- msr("surv.cindex")
evals <- trm("evals", n_evals = 15)
instance <- TuningInstanceSingleCrit$new(
  task = task_TV,
  learner = learner_RF,
  resampling = cv_inner,
  measure = measure,
  search_space = search_space,
  terminator = evals
)
tuner <- tnr("grid_search", resolution = 15)
tuner$optimize(instance)
#nodesize = 8
#0.7500753

#Coxph
Data_ML_TV <- as.data.frame(cbind(Label_PM$FUorRInterval[TrainingIndex], Label_PM$Recurrence[TrainingIndex], Ratio_Matrix_R1_impseqrob_combat_PM[TrainingIndex, Protein_Index_R2_retain]))
colnames(Data_ML_TV)[1] <- "time"
colnames(Data_ML_TV)[2] <- "status"

task_TV <- as_task_surv(Data_ML_TV, id = "right_censored",
                        time = "time", event = "status", type = "right")
print(task_TV)
learner <- lrn("surv.glmnet", maxit = 5000000, alpha = 1)
learner$param_set
learner
set.seed(1)
search_space <- ps(
  lambda = p_dbl(lower = 0.11, upper = 0.15)
)
cv_inner <- rsmp("cv", folds = 3)
measure <- msr("surv.cindex")
evals <- trm("evals", n_evals = 50)
instance <- TuningInstanceSingleCrit$new(
  task = task_TV,
  learner = learner,
  resampling = cv_inner,
  measure = measure,
  search_space = search_space,
  terminator = evals
)
tuner <- tnr("grid_search", resolution = 50)
tuner$optimize(instance)
(lambda <- instance$result_learner_param_vals$lambda)
#0.1116327
(alpha <- instance$result_learner_param_vals$alpha)
#1
instance$result_y
#0.82177

#Rsf
learner_RF <- lrn("surv.rfsrc", ntree = 1000, mtry = 39)
set.seed(1)
search_space <- ps(
  nodesize = p_int(lower = 1, upper = 15)
)
search_space
cv_inner <- rsmp("cv", folds = 3)
measure <- msr("surv.cindex")
evals <- trm("evals", n_evals = 15)
instance <- TuningInstanceSingleCrit$new(
  task = task_TV,
  learner = learner_RF,
  resampling = cv_inner,
  measure = measure,
  search_space = search_space,
  terminator = evals
)
tuner <- tnr("grid_search", resolution = 15)
tuner$optimize(instance)
#nodesize = 6

#Rsf
Data_ML_TV <- as.data.frame(cbind(Label_PM$FUorRInterval[TrainingIndex], Label_PM$Recurrence[TrainingIndex], Ratio_Matrix_R1_impseqrob_combat_PM[TrainingIndex,], ClinicalFeatures[TrainingIndex,]))
colnames(Data_ML_TV)[1] <- "time"
colnames(Data_ML_TV)[2] <- "status"
task_TV <- as_task_surv(Data_ML_TV, id = "right_censored",
                        time = "time", event = "status", type = "right")
print(task_TV)
learner_RF <- lrn("surv.rfsrc", ntree = 1000, mtry = 39)
set.seed(3)
search_space <- ps(
  nodesize = p_int(lower = 1, upper = 15)
)
search_space
cv_inner <- rsmp("cv", folds = 3)
measure <- msr("surv.cindex")
evals <- trm("evals", n_evals = 15)
instance <- TuningInstanceSingleCrit$new(
  task = task_TV,
  learner = learner_RF,
  resampling = cv_inner,
  measure = measure,
  search_space = search_space,
  terminator = evals
)
tuner <- tnr("grid_search", resolution = 15)
tuner$optimize(instance)
#nodesize = 7

#Importance
Importance_RF <- matrix(0,50,100)
for (i in 1:100) {
  set.seed(i)
  learner_RF <- lrn("surv.rfsrc", ntree = 1000, mtry = 39, nodesize = 6, importance = "permute")
  learner_RF$train(task_TV, row_ids = 1:50)
  print(learner_RF$predict(task_TV, row_ids = 1:50)$score(msr("surv.cindex")))
  Importance_RF[,i] <- as.vector(rownames(as.data.frame(learner_RF$importance()))[1:50])
}
rm(i)
table(as.vector(as.matrix(Importance_RF)))

#Importance
Importance_RF2 <- matrix(0,50,100)
for (i in 1:100) {
  set.seed(i)
  learner_RF <- lrn("surv.rfsrc", ntree = 1000, mtry = 39, nodesize = 7, importance = "permute")
  learner_RF$train(task_TV,row_ids = 1:50)
  print(learner_RF$predict(task_TV,row_ids = 1:50)$score(msr("surv.cindex")))
  Importance_RF2[,i] <- as.vector(rownames(as.data.frame(learner_RF$importance()))[1:50])
}
rm(i)
table(as.vector(as.matrix(Importance_RF2)))
