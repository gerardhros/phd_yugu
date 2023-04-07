# preparation regression between Qmax/Kl and soil properties
 # clear environment
 rm(list = ls())
 #########################XGBOOST###############################
  library(data.table);library(ggplot2);library(caret);library(mlr3hyperband);library(emoa)
  library(mlr3);library(mltools);library(xgboost);library(mlr3learners)
  library(mlr3measures); library(mlr3tuning);library(paradox);library(DALEX)
 
 # Load data
  d1 <- fread("data/experimentalData/paper1_lr_origional.csv")
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
  d2 <- d1[title2!=7& title2!=80][,c("ph_h2o","som","clay","soil_type_fao2","al","fe","totalp","qmax","klnew")]

  d2[,ln_qmax:= log(qmax)][,ln_k:= log(klnew)][,ln_qk:= log(qmax/klnew)][,ln_som:=log(som)][,ln_clay:=log(clay)][,ln_totalp:=log(totalp)]
  d2[, ln_al:= log(al)][,ln_fe:= log(fe)]
  
 #add new categorical variable for soil type to decrease the number of soil types
  
  d2[soil_type_fao2=='Calcisols'|soil_type_fao2=='Luvisols', soilty:= 'alfisols']
  d2[soil_type_fao2=='Anthrosols'|soil_type_fao2=='Gleysols'|soil_type_fao2=='Inceptisols', soilty:= 'inceptisols']
  d2[soil_type_fao2=='Chernozems'|soil_type_fao2=='Phaeozems', soilty:= 'mollisols']
  d2[soil_type_fao2=='Ferralosols'|soil_type_fao2=='Ferrosols', soilty:= 'oxisols']
  d2[soil_type_fao2=='Solonchaks'|soil_type_fao2=='Solonets', soilty:= 'saline']
  d2[soil_type_fao2=='Fluvisols', soilty:= 'fluvisols']
 
   # which columns are numeric
  
  cols <- colnames(d2[ , .SD, .SDcols = is.numeric])
  d2.mean <- d2[,lapply(.SD,mean,na.rm=T),.SDcols = cols]
  d2.sd <- d2[,lapply(.SD,sd,na.rm=T),.SDcols = cols] 

  # save mean and sd per variable
  
  d2.mean <- d2[,lapply(.SD,mean,na.rm=T),.SDcols = cols]
  d2.sd <- d2[,lapply(.SD,sd,na.rm=T),.SDcols = cols] 
 
  # scale
  
  d2[,c(cols) := lapply(.SD,scale),.SDcols = cols]
 
  # remove the variable not needed anymore
  
  cols <- c("som","clay","totalp","qmax","klnew","soil_type_fao2","al","fe")
  d2 <- as.data.table(na.omit(d2[,c(cols):=NULL]))
  
  # perform correlation test for all numeric variables
  
  library(ggcorrplot)
  
  # Calculate correlation matrix and corresponding P values
  corr_mat <- cor(d2[, -'soilty'])
  p_mat <- cor_pmat(d2[, -'soilty'], method = "pearson")
  
  # Create correlation plot 
  ggcorrplot(corr_mat, 
             p.mat = p_mat)
  
  # Perform chi-square tests for the correlation between soil type and all numeric variables
  
  chisq <- data.frame(variable = character(), p_value = numeric(), stringsAsFactors = FALSE)
  
  for (var in colnames(d2)) {
    if (is.numeric(d2[[var]])) {
      chisq_result <- chisq.test(table(d2$soilty, d2[[var]]))
      chisq <- rbind(chisq, data.frame(variable = var, p_value = chisq_result$p.value))
    } else {
      warning(sprintf("Skipping non-numeric variable '%s'", var))
    }
  }
  
  # View the results
  chisq
  

# make a linear regression model
 
  # one-hot code for soil type
  
  cols_soilty <- c('soilty')
  d2[,c(cols_soilty) := lapply(.SD,as.factor),.SDcols = cols_soilty]
  d2 <- mltools::one_hot(d2, cols = c("soilty"))

  # split in train and test set
  
  set.seed(123)
  fr.test <- 0.30
  rows.test <- sample(1:nrow(d2), size = fr.test * nrow(d2), replace = FALSE)
  rows.train <- which(! 1:nrow(d2) %in% rows.test)
 
  # make two separate data.tables with training and testing data
  
  dt.train <- d2[rows.train,]
  dt.test <- d2[rows.test,]
 
  # make a linear model on training data with Fe&Al
 
  glm_qmax <-lm(formula = ln_qmax ~ ln_totalp + ph_h2o:ln_totalp + ln_som:ln_al + 
                  ln_al, data = dt.train)
  
  glm_k <-lm(formula = ln_k ~ ln_totalp + ln_al + ph_h2o:ln_clay + ph_h2o:ln_al + 
                               ph_h2o:ln_fe + ln_clay:ln_al + ln_clay:ln_fe + ln_al:ln_fe, 
                         data = dt.train)
  
  glm_qk <-lm(formula = ln_qk ~ ln_clay:ln_som + ph_h2o:ln_fe + ln_totalp, data = dt.train)
 
  # make a linear model on training data without Fe&Al
  
  # ph_h2o*ln_som+ph_h2o*ln_clay+ph_h2o*ln_al+ph_h2o*ln_fe+ph_h2o*ln_totalp+ln_som*ln_clay+ ln_som*ln_al
  # +ln_som*ln_fe+ln_clay*ln_al+ln_clay*ln_fe+ln_al*ln_fe
  # 
  # glm_qmax <-lm(ln_qmax~ph_h2o+ln_clay*ln_som+ln_totalp+I((ln_totalp)^2)+soilty_alfisols+
  #                soilty_inceptisols+soilty_mollisols+soilty_oxisols+soilty_fluvisols,#+soilty_saline,
  #              data = dt.train)
  # glm_k <-lm(ln_k~ph_h2o+ln_clay*ln_som+ln_totalp+I((ln_totalp)^2)+soilty_alfisols+
  #             soilty_inceptisols+soilty_mollisols+soilty_oxisols+soilty_fluvisols,#+soilty_saline,
  #           data = dt.train)
  # glm_qk <-lm(ln_qk~ph_h2o+ln_clay*ln_som+ln_totalp) #combine residuals1 and 2 in one dataframe

# building xgboost model---------
 
 #make local copy and split dataset, selected relevant variables
 d3 <- copy(d2[,-c('ln_qmax','ln_qk')])
 #d3 <- copy(d2[,-c(10:14)][,-c('ln_k','ln_qk')])
 dt.train.xgb <- d3[rows.train,]
 dt.test.xgb <- d3[rows.test,]
 
 # set tuning methods
 target <- "ln_k"
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
  tuner = mlr3tuning::tnr("hyperband", eta = 2L)
 
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
                                y = dt.train.xgb$ln_k, 
                                label = "model_train_xgb")
  explainer.test.xgb <- explain(model, 
                               data = as.data.frame(dt.test.xgb[, .SD, .SDcols = !c(target)]), 
                               y = dt.test.xgb$ln_k, 
                               label = "model_test_xgb")
 
  cols.lm <- c('ph_h2o','ln_som','ln_totalp','ln_clay','soilty_alfisols','ln_al','ln_fe',
              'soilty_fluvisols','soilty_inceptisols','soilty_mollisols',
              'soilty_oxisols')
 
  explainer.train.lm <- explain(model = glm_k, 
                               data = as.data.frame(dt.train[, .SD, .SDcols = cols.lm]), 
                               y = dt.train$ln_k, 
                               label = "model_train_glm")
  explainer.test.lm <- explain(model =glm_k, 
                              data = as.data.frame(dt.test[, .SD, .SDcols = cols.lm]), 
                              y = dt.test$ln_k, 
                              label = "model_test_glm")
 
 # make VIP plot on training sets
  imp.xgb <- ingredients::feature_importance(explainer.train.xgb, 
                                            loss_function = loss_root_mean_square, 
                                            type = "difference")
  imp.lm <- ingredients::feature_importance(explainer.train.lm, 
                                           loss_function = loss_root_mean_square, 
                                           type = "difference")
 
  
  plot.vip <-ggplot_imp(imp.xgb, imp.lm)
  ggsave(plot = plot.vip,filename = 'products/lnk/withFeAl/plot2_vip_bar.png',width = 13.16, height = 8.90, units='cm')
  
  
 # make a residual plot on test
  res.xgb.test <- auditor::model_residual(explainer.test.xgb)
  res.xgb.train <- auditor::model_residual(explainer.train.xgb)
  res.lm.test <- auditor::model_residual(explainer.test.lm)
  res.lm.train <- auditor::model_residual(explainer.train.lm)  
 
  auditor::score_r2(explainer.test.lm)
  auditor::score_r2(explainer.test.xgb)
 
  plot.res <- ggplot_hist(res.xgb.train,res.xgb.test,res.lm.train,res.lm.test)
  ggsave(plot = plot.res,filename = 'products/lnk/withFeAl/plot_res.png',width = 13.16, height = 8.90, units='cm')
 
 # make a 1-to-1 plot of both regressions
  
  plot.mp <- ggplot_onetoone(res.xgb.train,res.xgb.test,res.lm.train,res.lm.test)
  plot.mp <- ggplot_onetoone(res.xgb.test,res.lm.test)
  
  ggsave(plot = plot.mp,filename = 'products/lnk/withFeAl/plot3_1-to-1.png',width = 13.16, height = 8.90, units='cm')
  
  require(patchwork)
  pcombi <- plot.vip + plot.mp
  ggsave(plot = pcombi,filename = 'products/lnk/withFeAl/plot2_combi.png',width = 13.16, height = 8.90, units='cm')
  
 # plot ALE plots
  
  ale.clay.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_clay')
  ale.clay.lm <- ingredients::accumulated_dependency(explainer.train.lm, 'ln_clay')
  ale.ph.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ph_h2o')
  ale.ph.lm <- ingredients::accumulated_dependency(explainer.train.lm, 'ph_h2o')
  ale.som.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_som')
  ale.som.lm <- ingredients::accumulated_dependency(explainer.train.lm, 'ln_som')
  ale.totalp.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_totalp')
  ale.totalp.lm <- ingredients::accumulated_dependency(explainer.train.lm, 'ln_totalp')
  
  # for regression with Fe+Al
  
  ale.feox.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_fe')
  ale.feox.lm <- ingredients::accumulated_dependency(explainer.train.lm, 'ln_fe')
  ale.alox.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_al')
  ale.alox.lm <- ingredients::accumulated_dependency(explainer.train.lm, 'ln_al')
  plot.ale <- ggplot_ale(ale.clay.xgb,ale.clay.lm,ale.ph.xgb,ale.ph.lm,
                         ale.som.xgb,ale.som.lm,ale.totalp.xgb,ale.totalp.lm,
                         ale.feox.xgb,ale.feox.lm,ale.alox.xgb,ale.alox.lm)
  ggsave(plot = plot.ale,filename = 'products/lnk/withFeAl/plot2_ale.png',width = 13.16, height = 8.90, units='cm')
  
  # calculate R2 and RMSE for all models
  library(caret)
  #
  trainpred.xgb <- predict(model, newdata = dt.train.xgb)
  dt.train.xgb[, predicted := trainpred.xgb]
  defaultSummary(data.frame(obs = dt.train.xgb$ln_k, pred = trainpred.xgb))
  #
  testpred.xgb <- predict(model, newdata = dt.test.xgb)
  dt.test.xgb[, predicted := testpred.xgb]
  defaultSummary(data.frame(obs = dt.test.xgb$ln_k, pred = testpred.xgb))
  #
  trainpred.glm <- predict(glm_k, newdata = dt.train)
  dt.train[, predicted := trainpred.glm]
  defaultSummary(data.frame(obs = dt.train$ln_k, pred = trainpred.glm))
  #
  testpred.glm <- predict(glm_k, newdata = dt.test)
  dt.test[, predicted := testpred.glm]
  defaultSummary(data.frame(obs = dt.test$ln_k, pred = testpred.glm))
  
# plot 1 - 1 for both regression models
  
  # ##XGBOOST
    plotfun <- function(){
      #combine residuals1 and 2 in one dataframe
      train <- data.table(observed = dt.train.xgb$ln_k,
                        predicted = dt.train.xgb$predicted,
                        type = "training")
      test <- data.table(observed = dt.test.xgb$ln_k,
                       predicted = dt.test.xgb$predicted,
                       type = "test")
      combined <- rbind(train, test)
      #create label
      dt.label <- data.table(
        R2 = c(caret::R2(train$observed, train$predicted), caret::R2(test$observed, test$predicted)),
        RMSE = c(caret::RMSE(train$observed, train$predicted) ,caret::RMSE(test$observed, test$predicted)),
        type = c("training", "test"))
      dt.label[, RMSE := round(RMSE, 2)]
      dt.label[, R2 := round(R2, 2)]#
      dt.label[, label := paste("R2 =", R2, "\nRMSE = ", RMSE)]
      #split
      train.label <- dt.label[type == "training"]
      test.label <- dt.label[type == "test"]
      #plot
      ggplot(combined, aes(x = observed, y = predicted, col = type)) + geom_point() +
          scale_color_manual(values  = c(`training` = "#e31a1c", `test` = "#2c7bb6")) +
          theme_bw() + labs(col = "Dataset") +
          geom_smooth(method = 'lm', se = FALSE, linetype = "solid") +
            theme(axis.text = element_text(family="A",size = 12,colour ="black" ),#
                  axis.title.x = element_text(family="A",size = 12,colour ="black",face = "bold"),#
                  axis.title.y  = element_text(family="A",size = 12,colour ="black",face = "bold"))+#
            theme(legend.position = c(0.85,0.15),#
                  legend.text = element_text(family="A", size = 9,face = "bold"),#
                  legend.title = element_text(family="A",size = 9,face = "bold"),#
                  legend.key =  element_rect(fill = "transparent", colour = NA))+#
          ggtitle("Log(kl)~with Fe&AL[XGBOOST]")+
            theme(plot.title = element_text(hjust = 0.5,family="A",size = 12 )) +
        xlim(-2,3)+ylim(-2,3)+
          geom_abline(intercept = 0,slope=1, colour='black',linetype=2, size=0.5)+
          geom_text(data = train.label, aes(x = -0.1, y = 2), label = train.label$label,
                inherit.aes = FALSE, col = "#e31a1c",family="A", size = 4,face = "bold") +
          geom_text(data = test.label, aes(x = -0.1, y = 1.5), label = test.label$label,
                inherit.aes = FALSE, col = "#2c7bb6",family="A", size = 4,face = "bold")
  }
    ggsave("products/lnk/withFeAl/plotEX_xgb.png",width = 13.16, height = 13.16, units='cm', dpi = 800)
  
  # ##GLM 
    plotfun <- function(){
       #combine residuals1 and 2 in one dataframe
       train <- data.table(observed = dt.train$ln_k,
                           predicted = dt.train$predicted,
                           type = "training")
       test <- data.table(observed = dt.test$ln_k,
                          predicted = dt.test$predicted,
                          type = "test")
       combined <- rbind(train, test)
       #create label
       dt.label <- data.table(
       R2 = c(caret::R2(train$observed, train$predicted), caret::R2(test$observed, test$predicted)),
       RMSE = c(caret::RMSE(train$observed, train$predicted) ,caret::RMSE(test$observed, test$predicted)),
       type = c("training", "test"))
       dt.label[, RMSE := round(RMSE, 2)]
       dt.label[, R2 := round(R2, 2)]#
       dt.label[, label := paste("R2 =", R2, "\nRMSE = ", RMSE)]
       #split
       train.label <- dt.label[type == "training"]
       test.label <- dt.label[type == "test"]
       #plot
       ggplot(combined, aes(x = observed, y = predicted, col = type)) + geom_point() +
          scale_color_manual(values  = c(`training` = "#e31a1c", `test` = "#2c7bb6")) +
          theme_bw() + labs(col = "Dataset") +
          geom_smooth(method = 'lm', se = FALSE, linetype = "solid") +
          theme(axis.text = element_text(family="A",size = 12,colour ="black" ),#
                axis.title.x = element_text(family="A",size = 12,colour ="black",face = "bold"),#
                axis.title.y  = element_text(family="A",size = 12,colour ="black",face = "bold"))+#
          theme(legend.position = c(0.85,0.15),#
                legend.text = element_text(family="A", size = 9,face = "bold"),#
                legend.title = element_text(family="A",size = 9,face = "bold"),#
                legend.key =  element_rect(fill = "transparent", colour = NA))+#
          ggtitle("Log(Kl)~with Fe&AL [Linear]")+
         xlim(-2,3)+ylim(-2,3)+
          theme(plot.title = element_text(hjust = 0.5,family="A",size = 12 )) +
          geom_abline(intercept = 0,slope=1, colour='black',linetype=2, size=0.5)+
          geom_text(data = train.label, aes(x = -0.1, y = 2), label = train.label$label,
                    inherit.aes = FALSE, col = "#e31a1c",family="A", size = 4,face = "bold") +
          geom_text(data = test.label, aes(x = -0.1, y = 1.5), label = test.label$label,
                    inherit.aes = FALSE, col = "#2c7bb6",family="A", size = 4,face = "bold")
  }
  
    ggsave("products/lnk/withFeAl/plotEX_glm.png",width = 13.16, height = 13.16, units='cm', dpi = 800)
  

 
