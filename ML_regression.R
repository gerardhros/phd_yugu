##through random forest and xgboost to do the regression individually
#########################XGBOOST###############################
library(data.table);library(ggplot2);library(caret);
library(mlr3);library(mltools);library(xgboost);library(mlr3);library(mlr3learners)
library(mlr3measures); library(mlr3tuning);library(paradox);library(DALEX);library(ggplot2)
library(mlr3hyperband);library(emoa)
# Load data
d1 <- fread("C:/Users/gu021/OneDrive - Wageningen University & Research/1 WUR_work/2 meta/paper1_lr_origional.csv")
length(unique(d1$title))
# Rename columns
cols_old <- colnames(d1)
cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
setnames(d1,cols_old,cols_new)
d1[,alfe:= al+fe]
d1[,title2 := .GRP,by='title']
d1[,logq:= log(qmax)]
d2 <- d1[-c(68:73),][title2!=7][,c("ph_h2o","som","clay","soil_type_fao2","totalp", "logq")]
cols_soilty <- c('soil_type_fao2')
d2[,c(cols_soilty) := lapply(.SD,as.factor),.SDcols = cols_soilty][]
d3 <- mltools::one_hot(d2, cols = c("soil_type_fao2"))

#--------
target <- "logq"
tune_method <- "hyperband"#methods for best parameter

# if model is trained before, then change this in TRUE
runMLmodel <- FALSE
dt <- as.data.table(d3)
set.seed(123)
fr.test <- 0.30
rows.test <- sample(1:nrow(dt), size = fr.test * nrow(dt), replace = FALSE)
rows.train <- which(! 1:nrow(dt) %in% rows.test)

dt.train <- dt[rows.train, ]
dt.test <- dt[rows.test, ]

if(!runMLmodel){
  # Set the task
  task.train <- TaskRegr$new(id = target, backend = dt.train, target = target)
  
  # Set the learner
  learner <- lrn("regr.xgboost")
  
  # Set the parameters of the learner
  learner$param_set$values <- list(
    eval_metric = "rmsle",
    verbose = 0,
    nthread = 8,
    nrounds=1000,
    early_stopping_rounds = 100,
    early_stopping_set = "train"
  )
  
  ps.tune <- ParamSet$new(list(
    ParamInt$new("nrounds", lower = 100, upper = 1600, default = 750, tags = "budget"),
    ParamInt$new("max_depth", lower = 5, upper = 20, default = 12),
    ParamDbl$new("min_child_weight", lower = 1, upper = 5, default =1),
    ParamDbl$new("subsample", lower = 0.5, upper = 1, default = 0.5),
    ParamDbl$new("colsample_bytree", lower = 0.5, upper = 1, default = 0.5),
    ParamDbl$new("eta", lower = 2^-9, upper = 2^-3, default = 2^-5),
    ParamDbl$new("gamma", lower = 0, upper = 5, default = 0)
  ))
  
  # Set the measures
  measures <- list(msr("regr.rmse"), msr("time_train"))
  
  # Set the resampling
  resampling <- rsmp("cv", folds = 3L)
  
  # Set the tuning
  terminator <- trm("run_time", secs = 60 * 60)
  
  # Tune the model ---------------------------------------------------------娴兼ê瀵插Ο鈥崇€?
  
  if (tune_method == "random") {
    
    tuner <- mlr3tuning::tnr("random_search")
    
    at <- AutoTuner$new(
      learner = learner,
      resampling = resampling,
      measures = measures,
      tune_ps = ps.tune,
      terminator = terminator,
      tuner = tuner
    )
    at$train(task = task.train)
    
    model <- at$model$learner
    
  } else if (tune_method == "hyperband") {
    
    
    tuner = tnr("hyperband", eta = 2L)
    
    inst <- TuningInstanceSingleCrit$new(
      task = task.train,
      learner = learner,
      resampling = resampling,
      measure = msr("regr.rmse"),
      search_space = ps.tune,
      terminator = terminator
    )
    
    result <- tuner$optimize(inst)
    
    # Set the best hyperparameters
    learner$param_set$values <- inst$result_learner_param_vals
    
    # Train the model
    model <- learner$train(task = task.train)
  }
  
  # save the model for later use
  saveRDS(model,'ml_model.rds')
  # xgb.save(model,'ml_model.rds')
} else{
  #model <- xgb.load('ml_model.rds')
  model <- readRDS('ml_model.rds')
}

# Interpret the model -----------------------------------------------------

custom_predict <- function(model, newdata) {
  x <- predict(model, newdata = newdata)
  return(x)
}
# Create the explainer
explainer.train <- explain(model, 
                           data = as.data.frame(dt.train[, .SD, .SDcols = !c(target)]), 
                           y = dt.train$logq, label = "model_train")
explainer.test <- explain(model, 
                          data = as.data.frame(dt.test[, .SD, .SDcols = !c(target)]), 
                          y = dt.test$logq, label = "model_test")


explainer.all <- explain(model, 
                         data = as.data.frame(dt[, .SD, .SDcols = !c(target)]), 
                         y = dt$logq, label = "model_all")



#--------------
performance1 <- model_performance(explainer.train)
performance2 <- model_performance(explainer.test)
performance3<- model_performance(explainer.all)

library(caret)
#
trainpred <- predict(model, newdata = dt.train)
dt.train[, predicted := trainpred]
defaultSummary(data.frame(obs = dt.train$logq, pred = trainpred))
#
trainpred <- predict(model, newdata = dt.test)
dt.test[, predicted := trainpred]
defaultSummary(data.frame(obs = dt.test$logq, pred = trainpred))


plotfun <- function(){
  #combine residuals1 and 2 in one dataframe
  train <- data.table(observed = dt.train$logq,
                      predicted = dt.train$predicted,
                      type = "training")
  test <- data.table(observed = dt.test$logq,
                     predicted = dt.test$predicted,
                     type = "test")
  combined <- rbind(train, test)
  #create label
  dt.label <- data.table(
    R2 = c(caret::R2(train$observed, train$predicted), caret::R2(test$observed, test$predicted)),
    RMSE = c(caret::RMSE(train$observed, train$predicted) ,caret::RMSE(test$observed, test$predicted)),
    type = c("training", "test"))
  
  dt.label[, RMSE := round(RMSE, 2)]
  dt.label[, R2 := round(R2, 2)]
  dt.label[, label := paste("R2 =", R2, "\nRMSE = ", RMSE)]
  #split
  train.label <- dt.label[type == "training"]
  test.label <- dt.label[type == "test"]
  
  #plot
  ggplot(combined, aes(x = observed, y = predicted, col = type)) + geom_point() +
    scale_color_manual(values  = c(`training` = "#e31a1c", `test` = "#1a9641")) +
    theme_bw() + labs(col = "Dataset") +
    geom_smooth(method = 'lm', se = FALSE, linetype = "solid") +
    theme(axis.text = element_text(family="A",size = 12,colour ="black" ),#
          axis.title.x = element_text(family="A",size = 12,colour ="black",face = "bold"),#
          axis.title.y  = element_text(family="A",size = 12,colour ="black",face = "bold"))+#
    theme(legend.position = c(0.85,0.15),#
          legend.text = element_text(family="A", size = 9,face = "bold"),#
          legend.title = element_text(family="A",size = 9,face = "bold"),#
          legend.key =  element_rect(fill = "transparent", colour = NA))+#
    #xlim(0,1.5)+ylim(0,1.5)+
    ggtitle("Log(Qmax)~(SOM+clay+ph+totalp+soil type)")+
    theme(plot.title = element_text(hjust = 0.5,family="A",size = 12 )) +
    geom_abline(intercept = 0,slope=1, colour='black',linetype=2, size=0.5)
  # geom_text(data = train.label, aes(x = 0.2, y = 1.3), label = train.label$label, 
  #           inherit.aes = FALSE, col = "#e31a1c",family="A", size = 4,face = "bold") +
  # geom_text(data = test.label, aes(x = 0.2, y = 1.1), label = test.label$label, 
  #           inherit.aes = FALSE, col = "#1a9641",family="A", size = 4,face = "bold")
}