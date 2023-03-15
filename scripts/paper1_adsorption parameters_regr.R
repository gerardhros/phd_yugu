# preparation regression between Qmax/Kl and soil properties
 # clear environment
 rm(list = ls())
 #########################XGBOOST###############################
 library(data.table);library(ggplot2);library(caret);library(mlr3hyperband);library(emoa)
 library(mlr3);library(mltools);library(xgboost);library(mlr3learners)
 library(mlr3measures); library(mlr3tuning);library(paradox);library(DALEX)
 
 # Load data
 d1 <- fread("data/paper1_lr_origional.csv")
 length(unique(d1$title))
 
 # Rename columns
 cols_old <- colnames(d1)
 cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
 setnames(d1,cols_old,cols_new)
 
 # remove strange data from specific study
 d1[,title2 := .GRP,by='title']
 
 # make sleected columns with numbers also numeric
 cols.num <- c("ph_h2o","som","clay","caco3","al","fe","qmax","klnew", "totalp")
 d1[,c(cols.num) := lapply(.SD,as.numeric),.SDcols = cols.num]
 
 # paper in title 2==7 is natural soil and no human intervention;titl2==80 is dong fang's unpublished paper
 # selected variables for regression
 d2 <- d1[title2!=7& title2!=80][,c("ph_h2o","som","clay","soil_type_fao2","totalp","qmax","klnew")]
 d2[,ln_qmax:= log(qmax)][,ln_k:= log(klnew)][,ln_som:=log(som)][,ln_clay:=log(clay)][,ln_totalp:=log(totalp)]
 
 #add new categorical variable for soil type to decrease the number of soil types
 d2[soil_type_fao2=='Calcisols'|soil_type_fao2=='Luvisols', soilty:= 'alfisols']
 d2[soil_type_fao2=='Anthrosols'|soil_type_fao2=='Gleysols'|soil_type_fao2=='Inceptisols', soilty:= 'inceptisols']
 d2[soil_type_fao2=='Chernozems'|soil_type_fao2=='Phaeozems', soilty:= 'mollisols']
 d2[soil_type_fao2=='Ferralosols'|soil_type_fao2=='Ferrosols', soilty:= 'oxisols']
 d2[soil_type_fao2=='Solonchaks'|soil_type_fao2=='Solonets', soilty:= 'saline']
 d2[soil_type_fao2=='Fluvisols', soilty:= 'fluvisols']
 

 # which columns are numeric
 cols <- colnames(d2[ , .SD, .SDcols = is.numeric])
 
 # save mean and sd per variable
 d2.mean <- d2[,lapply(.SD,mean,na.rm=T),.SDcols = cols]
 d2.sd <- d2[,lapply(.SD,sd,na.rm=T),.SDcols = cols] 
 
 # scale
 d2[,c(cols) := lapply(.SD,scale),.SDcols = cols]
 
 # remove the variable not needed anymore
 cols <- c("som","clay","totalp","qmax","klnew","soil_type_fao2")
 d2 <- as.data.table(na.omit(d2[,c(cols):=NULL]))
 
# make a linear regression model
 
 #one-hot code for soil type
 cols_soilty <- c('soilty')
 d2[,c(cols_soilty) := lapply(.SD,as.factor),.SDcols = cols_soilty]
 d2 <- mltools::one_hot(d2, cols = c("soilty"))
 
 #since always with error, i then transfer 'soil type' into 'factor' again but doesn't work
 d2[, c('soilty_alfisols', 'soilty_fluvisols', 'soilty_inceptisols', 'soilty_mollisols', 'soilty_oxisols', 'soilty_saline') :=
      .(as.factor(soilty_alfisols), as.factor(soilty_fluvisols), as.factor(soilty_inceptisols),
        as.factor(soilty_mollisols), as.factor(soilty_oxisols), as.factor(soilty_saline))]

 #split in train and test set
 set.seed(123)
 fr.test <- 0.30
 rows.test <- sample(1:nrow(d2), size = fr.test * nrow(d2), replace = FALSE)
 rows.train <- which(! 1:nrow(d2) %in% rows.test)
 
 # make two separate data.tables with training and testing data
 dt.train <- d2[rows.train,]
 dt.test <- d2[rows.test, ]
 
# make a linear model on training data for only total P

 glm_qmax <-lm(ln_qmax~ph_h2o+ln_clay*ln_som+ln_totalp+I((ln_totalp)^2)+soilty_alfisols+
             soilty_inceptisols+soilty_mollisols+soilty_oxisols+soilty_fluvisols,#+soilty_saline,
           data = dt.train)

# building xgboost model---------
 
 #make local copy and split dataset, selected relevant variables
 d3 <- copy(d2[,-c('ln_k')])
 dt.train.xgb <- d3[rows.train,]
 dt.test.xgb <- d3[rows.test,]
 
 # set tuning methods
 target <- "ln_qmax"
 tune_method <- "hyperband"#methods for best parameter

 ###########
 # Set the task
 task.train <- TaskRegr$new(id = target, backend = dt.train.xgb, target = target)
 
 # Set the learner
 learner <-  mlr3::lrn("regr.xgboost")
 
 # Set the parameters of the learner
 learner$param_set$values <- list(
   #eval_metric = "rmsle",
   verbose = 0,
   nthread = 8,
   early_stopping_rounds = 12,
   early_stopping_set = "train"
  )
 
  ps.tune <- paradox::ParamSet$new(list(
   ParamInt$new("nrounds", lower = 100, upper = 1600, default = 750, tags = "budget"),
   ParamInt$new("max_depth", lower = 5, upper = 20, default = 12),
   ParamDbl$new("min_child_weight", lower = 1, upper = 5, default = 1),
   ParamDbl$new("subsample", lower = 0.5, upper = 1, default = 0.5),
   ParamDbl$new("colsample_bytree", lower = 0.5, upper = 1, default = 0.5),
   ParamDbl$new("eta", lower = 2^-9, upper = 2^-3, default = 2^-5),
   ParamDbl$new("gamma", lower = 0, upper = 5, default = 0)
  ))
 
 # Set the measures
 measures <- list(msr("regr.rmse"), msr("time_train"))
 
 # Set the resampling
 resampling <- mlr3::rsmp("cv", folds = 3L)
 
 # Set the tuning
 terminator <- mlr3tuning::trm("run_time", secs = 2 * 60 * 60)
 
 # Tune the model ---------------------------------------------------------
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
 
 
 # Interpret the model -----------------------------------------------------
 
 custom_predict <- function(model, newdata) {
   x <- predict(model, newdata = newdata)
   return(x)
 }
 
 
# Create the explainers for both models
 explainer.train.xgb <- explain(model, 
                                data = as.data.frame(dt.train.xgb[, .SD, .SDcols = !c(target)]), 
                                y = dt.train.xgb$ln_qmax, 
                                label = "model_train_xgb")
 explainer.test.xgb <- explain(model, 
                               data = as.data.frame(dt.test.xgb[, .SD, .SDcols = !c(target)]), 
                               y = dt.test.xgb$ln_qmax, 
                               label = "model_test_xgb")
 
 cols.lm <- c('ph_h2o','ln_som','ln_totalp','ln_clay','soilty_alfisols',
              'soilty_fluvisols','soilty_inceptisols','soilty_mollisols',
              'soilty_oxisols')
 explainer.train.lm <- explain(model = glm_qmax, 
                               data = as.data.frame(dt.train[, .SD, .SDcols = cols.lm]), 
                               y = dt.train$ln_qmax, 
                               label = "model_train_glm")
 explainer.test.lm <- explain(model =glm_qmax, 
                              data = as.data.frame(dt.test[, .SD, .SDcols = cols.lm]), 
                              y = dt.test$ln_qmax, 
                              label = "model_test_glm")
 
 # make VIP plot on training sets
 imp.xgb <- ingredients::feature_importance(explainer.train.xgb, 
                                            loss_function = loss_root_mean_square, 
                                            type = "difference")
 imp.lm <- ingredients::feature_importance(explainer.train.lm, 
                                           loss_function = loss_root_mean_square, 
                                           type = "difference")
 
 plot.vip <-plot(imp.xgb, imp.lm)
 # make a residual plot on test
 res.xgb.test <- auditor::model_residual(explainer.test.xgb)
 res.xgb.train <- auditor::model_residual(explainer.train.xgb)
 res.lm.test <- auditor::model_residual(explainer.test.lm)
 res.lm.train <- auditor::model_residual(explainer.train.lm)  
 
 auditor::score_r2(explainer.test.lm)
 auditor::score_r2(explainer.test.xgb)
 
 plot.res <- plot(res.xgb.train,res.xgb.test)
 ggsave(plot = plot.res,filename = 'products/plot_res.png',width = 13.16, height = 8.90, units='cm')
 
 library(caret)
 #
 trainpred.xgb <- predict(model, newdata = dt.train.xgb)
 dt.train.xgb[, predicted := trainpred.xgb]
 defaultSummary(data.frame(obs = dt.train.xgb$ln_qmax, pred = trainpred.xgb))
 #
 testpred.xgb <- predict(model, newdata = dt.test.xgb)
 dt.test.xgb[, predicted := testpred.xgb]
 defaultSummary(data.frame(obs = dt.test.xgb$ln_qmax, pred = testpred.xgb))
 #
 trainpred.glm <- predict(glm_qmax, newdata = dt.train)
 dt.train[, predicted := trainpred.glm]
 defaultSummary(data.frame(obs = dt.train$ln_qmax, pred = trainpred.glm))
 #
 testpred.glm <- predict(glm_qmax, newdata = dt.test)
 dt.test[, predicted := testpred.glm]
 defaultSummary(data.frame(obs = dt.test$ln_qmax, pred = testpred.glm))
 
#visualization 

 plotfun <- function(){
  #combine residuals1 and 2 in one dataframe
  train <- data.table(observed = dt.train$ln_q,
                      predicted = dt.train$predicted,
                      type = "training")
  
  test <- data.table(observed = dt.test$ln_q,
                     predicted = dt.test$predicted,
                     type = "test")
  
  combined <- rbind(train, test)
  
  #create label
  dt.label <- data.table(
    R2 = c(caret::R2(train$observed, train$predicted), caret::R2(test$observed, test$predicted)),
    RMSE = c(caret::RMSE(train$observed, train$predicted) ,caret::RMSE(test$observed, test$predicted)),
    type = c("training", "test"))
  
  dt.label[, RMSE := round(RMSE, 2)]
  dt.label[, R2 := round(R2, 2)]#http://127.0.0.1:12767/graphics/plot_zoom_png?width=1200&height=900
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

    ggtitle("Log(Qmax)~(SOM+clay+ph+totalp+soil type)")+
    theme(plot.title = element_text(hjust = 0.5,family="A",size = 12 )) +
    geom_abline(intercept = 0,slope=1, colour='black',linetype=2, size=0.5)
  # geom_text(data = train.label, aes(x = 0.2, y = 1.3), label = train.label$label, 
  #           inherit.aes = FALSE, col = "#e31a1c",family="A", size = 4,face = "bold") +
  # geom_text(data = test.label, aes(x = 0.2, y = 1.1), label = test.label$label, 
  #           inherit.aes = FALSE, col = "#1a9641",family="A", size = 4,face = "bold")
}

