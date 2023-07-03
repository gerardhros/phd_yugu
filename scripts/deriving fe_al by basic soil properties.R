 ##explore the relationship between Fe&Al and SOM & CLAY at global level
 
 ##data: quzhou + qiyang inventory data + switherlands inventory data

 # clear environment
   rm(list = ls())
 #########################XGBOOST###############################
   library(data.table);library(ggplot2);library(caret);library(mlr3hyperband);library(emoa)
   library(mlr3);library(mltools);library(xgboost);library(mlr3learners)
   library(mlr3measures); library(mlr3tuning);library(paradox);library(DALEX)
   #load data
   
   d1 <- fread("data/experimentalData/data for deriving FeAl.csv")
   
   # Rename columns 
   
   cols_old <- colnames(d1)
   cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
   setnames(d1,cols_old,cols_new)
   
   # transfer fe&al from mg/kg to mmol/kg
   d1[,fe_mmol:= fe/56][,al_mmol:= al/27][,id:= NULL][,id_other:=NULL][,source:= NULL][,pox:=NULL][,polsen:= NULL]
   
   # multiple variable outlier test

   kmeans.result <- kmeans(as.data.table(na.omit(d1[,-c(4:5)])),centers = 5)
   
   centers <- kmeans.result$centers[kmeans.result$cluster, ]
   distances<-sqrt(rowSums((as.data.table(na.omit(d1[,-c(4:5)]))-centers)^2))
   out<-order(distances,decreasing=TRUE)[1:200]
   
   ##manually delete datapoints where Fe/Alox in mmol is higher than 400, need reference, 400 is from the boxplot of Fe/Alox in mmol
   
   d2 <- d1[fe_mmol < 400][som < 200][clay < 60][ al_mmol < 400][,ln_fe:= log(fe/56)][,ln_al:= log(al/27)][,al:= NULL][,fe:= NULL]#[,fe_mmol:= NULL][,al_mmol:= NULL]
   
   # Loop over variables
   for (var in c("fe_mmol", "al_mmol", "ph", "clay","som")) {
     
     # Set up plot grid
     par(mfrow = c(1, 2))
     
     # Histogram
     hist(d2[[var]], main = paste("Histogram of", var))
     
     # Boxplot
     boxplot(d2[[var]], main = paste("Boxplot of", var))
     
     # Normal Q-Q plot
     #qqnorm(d2[[var]], main = paste("Normal Q-Q plot of", var))
     
     # Save combined plot
     combined_plot <- recordPlot()
     
     # Export combined plot as png
     png(paste(var, "combined_plot.png", sep = "_"))
     replayPlot(combined_plot)
     dev.off()
     
   }
   
   variables <- c("som", "ph", "clay", "fe_mmol", "al_mmol")
   
   for (variable in variables) {
     plot_title <- switch(variable,
                          "som" = "SOM[mg/kg]",
                          "ph" = "pH",
                          "clay" = "Clay content[%]",
                          "fe_mmol" = "Feox[mmol/kg]",
                          "al_mmol" = "Alox[mmol/kg]")
     
     plot <- ggplot(d2) +
       geom_histogram(aes(x = !!rlang::sym(variable))) +
       theme_bw() +
       labs(x = plot_title, y = "Frequency")
     
     print(plot)
   }
   
   
   
   
   
   
   # make nice boxplot
   
   ggplot(d2)+
     geom_histogram(aes(x=som))+
     theme_bw()+
     labs(x="SOM[mg/kg]",y="Frequency")
   
   ggplot(d2)+
     geom_histogram(aes(x=ph))+
     theme_bw()+
     labs(x="pH",y="Frequency")
   
   ggplot(d2)+
     geom_histogram(aes(x=clay))+
     theme_bw()+
     labs(x="Clay content[%]",y="Frequency")
   
   ggplot(d2)+
     geom_histogram(aes(x=fe_mmol))+
     theme_bw()+
     labs(x="Feox[mmol/kg]",y="Frequency")
     
   ggplot(d2)+
     geom_histogram(aes(x=al_mmol))+
     theme_bw()+
     labs(x="Alox[mmol/kg]",y="Frequency")
   
   
   
   # split in train and test set
   
   set.seed(123)
   fr.test <- 0.30
   rows.test <- sample(1:nrow(d2), size = fr.test * nrow(d2), replace = FALSE)
   rows.train <- which(! 1:nrow(d2) %in% rows.test)
   
   # make two separate data.tables with training and testing data
   
   d3 <- copy(d2[,-'ln_al'])
   
   dt.train.xgb <- d3[rows.train,]
   dt.test.xgb <- d3[rows.test,]
   
   # set tuning methods
   target <- "ln_fe"
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
   
   # Save the model
   saveRDS(model, file = "dev/model_logfetest.rds")
   
   # Interpret the model -----------------------------------------------------
   
   custom_predict <- function(model, newdata) {
     x <- predict(model, newdata = newdata)
     return(x)
   }
 
   # Create the explainers for both models
   explainer.train.xgb <- explain(model, 
                                  data = as.data.frame(dt.train.xgb[, .SD, .SDcols = !c(target)]), 
                                  y = dt.train.xgb$ln_fe, 
                                  label = "model_train_xgb")
   explainer.test.xgb <- explain(model, 
                                 data = as.data.frame(dt.test.xgb[, .SD, .SDcols = !c(target)]), 
                                 y = dt.test.xgb$ln_fe, 
                                 label = "model_test_xgb")
   
   
   # calculate R2 and RMSE for all models
   library(caret)
   #
   trainpred.xgb <- predict(model, newdata = dt.train.xgb)
   dt.train.xgb[, predicted := trainpred.xgb]
   defaultSummary(data.frame(obs = dt.train.xgb$ln_fe, pred = trainpred.xgb))
   #
   testpred.xgb <- predict(model, newdata = dt.test.xgb)
   dt.test.xgb[, predicted := testpred.xgb]
   defaultSummary(data.frame(obs = dt.test.xgb$ln_fe, pred = testpred.xgb))
   
   
   # plot 1 - 1 for both regression models
   
   # ##XGBOOST
   plotfun <- function(){
     #combine residuals1 and 2 in one dataframe
     train <- data.table(observed = dt.train.xgb$ln_fe,
                         predicted = dt.train.xgb$predicted,
                         type = "training")
     test <- data.table(observed = dt.test.xgb$ln_fe,
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
       ggtitle("Log(Fe) with all variables")+
       theme(plot.title = element_text(hjust = 0.5,family="A",size = 12 )) +
       xlim(2,6)+ylim(2,6)+
       geom_abline(intercept = 0,slope=1, colour='black',linetype=2, size=0.5)+
       geom_text(data = train.label, aes(x =2.5, y = 5.5), label = train.label$label,
                 inherit.aes = FALSE, col = "#e31a1c",family="A", size = 4,face = "bold") +
       geom_text(data = test.label, aes(x = 2.5, y = 5), label = test.label$label,
                 inherit.aes = FALSE, col = "#2c7bb6",family="A", size = 4,face = "bold")
   }
   
   