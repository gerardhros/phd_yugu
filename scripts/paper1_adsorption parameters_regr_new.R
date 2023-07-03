# clear environment
 rm(list = ls())
#########################XGBOOST###############################
 library(data.table);library(ggplot2);library(caret);library(mlr3hyperband);library(emoa)
 library(mlr3);library(mltools);library(xgboost);library(mlr3learners)
 library(mlr3measures); library(mlr3tuning);library(paradox);library(DALEX)

# load data for clay from HWSD for the missing values for database

 clay_hwsd <- fread("data/experimentalData/meta_clay_missing.csv")
 
# load data for meta database
 
 metadata0510 <- fread("data/experimentalData/paper1_origional_new.csv")
 
# input missing clay data from HWSD to the meta database
 
 meta <- merge( clay_hwsd,metadata0510, by = "ID_new", all= T)
 meta[,clay := fifelse(is.na(clay), clay_missing1, clay)][, clay:= ifelse(clay == 0, NA, clay)]->d1

# Rename columns 
 
 cols_old <- colnames(d1)
 cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
 setnames(d1,cols_old,cols_new)

# calculation for SOC to SOM; log transform 
 
 cols.num <- c("ph","som","clay","ptotal","qmax","kl", "caco3","alox",
               "feox","polsen")
 d1[,c(cols.num) := lapply(.SD,as.numeric),.SDcols = cols.num] 
 d1[,som:= fifelse(is.na(som), soc*1.724,som)][`landuse type` != 'sediment'][`kl_unit` == 'L/mg']-> d2
 

 
# outlier test and data visualization
 
 # Loop over variables
 for (var in c("ph", "som", "clay", "ptotal","qmax" )) {
   
   # Set up plot grid
   par(mfrow = c(1, 3))
   
   # Histogram
   hist(d5[[var]], main = paste("Histogram of", var))
   
   # Boxplot
   boxplot(d5[[var]], main = paste("Boxplot of", var))
   
   # Normal Q-Q plot
   qqnorm(d5[[var]], main = paste("Normal Q-Q plot of", var))
   
   # Save combined plot
   combined_plot <- recordPlot()
   
   # Export combined plot as png
   png(paste(var, "combined_plot.png", sep = "_"))
   replayPlot(combined_plot)
   dev.off()
   
 }
 
 # delete Jingguang Zhang & Guoping Wang since adsorption experi is in different temperature instead of daily temperature
 # delete Qiujun Wang since 2 points with Quite high Qmax (10000mg/kg), which is almost 80 times higher than Qmax from other adsorption equation
 # why delete clay > 60??
  d3 <- d2[author != 'Jingguang Zhang' & author != 'Qiujun Wang' & author != 'Wang, Guo-Ping'][qmax <= 3000][som <= 100][ptotal <= 4][clay <=60][!is.na(ph)]

 # basic calculation 
 d4 <- d3[,ln_qmax:= log(qmax)][,ln_kl:= log(kl)][,ln_som:= log(som)][,ln_clay:=log(clay)][,ln_totalp:=log(ptotal)][,ln_caco3:= log(caco3)][,ln_fe:= log(feox/56)][,ln_al:= log(alox/27)]
 
 # selected related columns

 # d5 <- d4[,c("soilty_usda","ph","ln_som","ln_clay", "ln_totalp","ln_qmax","ln_kl","ln_caco3","ln_fe","ln_al")]
 # d5 <- d4[,c("ph","ln_som","ln_clay", "ln_totalp","ln_qmax","ln_kl","soilty_usda")]
 d5 <- d4[,c("ph","ln_som","ln_clay", "ln_totalp","ln_qmax","soilty_usda")]
 
 # code for categorical variable
 
 cols_soilty <- c('soilty_usda')
 
 d5[,c(cols_soilty) := lapply(.SD,as.factor),.SDcols = cols_soilty]
 d6 <- mltools::one_hot(d5, cols = c("soilty_usda"))

 # split in train and test set by random sampling
 
 set.seed(123)
 fr.test <- 0.30
 rows.test <- sample(1:nrow(d5), size = fr.test * nrow(d5), replace = FALSE)
 rows.train <- which(! 1:nrow(d5) %in% rows.test)
 
 # make two separate data.tables with training and testing data

 
 dt.train.xgb <- d6[rows.train,]
 dt.test.xgb <- d6[rows.test,]

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
 #saveRDS(model, file = "dev/0518model.rds")
 
 # Interpret the model -----------------------------------------------------
 
 custom_predict <- function(model, newdata) {
   x <- predict(model, newdata = newdata)
   return(x)
 }
 
 # Load the saved model
 #model0 <- readRDS(file = "dev/0518model_nofealca.rds")
 # Create the explainers for both models
 explainer.train.xgb <- explain(model, 
                                data = as.data.frame(dt.train.xgb[, .SD, .SDcols = !c(target)]), 
                                y = dt.train.xgb$ln_qmax, 
                                label = "model_train_xgb")
 explainer.test.xgb <- explain(model, 
                               data = as.data.frame(dt.test.xgb[, .SD, .SDcols = !c(target)]), 
                               y = dt.test.xgb$ln_qmax, 
                               label = "model_test_xgb")
 

 
 # make VIP plot on training sets
 imp.xgb <- ingredients::feature_importance(explainer.train.xgb, 
                                            loss_function = loss_root_mean_square, 
                                            type = "difference")

 
 plot.vip <-ggplot_imp(imp.xgb)
 #ggsave(plot = plot.vip,filename = 'products/lnk/withoutFeAl/0410plot2_vip_bar.png',width = 13.16, height = 8.90, units='cm')
 
 
 # make a residual plot on test
 res.xgb.test <- auditor::model_residual(explainer.test.xgb)
 
 res.xgb.train <- auditor::model_residual(explainer.train.xgb)

 

 auditor::score_r2(explainer.test.xgb)
 
 plot.res <- ggplot_hist(res.xgb.train,res.xgb.test)
 #ggsave(plot = plot.res,filename = 'products/lnk/withoutFeAl/0410plot_res.png',width = 13.16, height = 8.90, units='cm')
 
 # make a 1-to-1 plot of both regressions
 
 plot.mp <- ggplot_onetoone(res.xgb.train,res.xgb.test)

 
 #ggsave(plot = plot.mp,filename = 'products/lnk/withoutFeAl/0410plot3_1-to-1.png',width = 13.16, height = 8.90, units='cm')
 
 require(patchwork)
 pcombi <- plot.vip + plot.mp
 #ggsave(plot = pcombi,filename = 'products/lnk/withoutFeAl/0410plot2_combi.png',width = 13.16, height = 8.90, units='cm')
 
 # plot ALE plots
 
 ale.clay.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_clay')

 ale.ph.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ph')
 
 ale.som.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_som')

 ale.totalp.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_totalp')
 
 # ale.feox.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_fe')
 # 
 # ale.alox.xgb <- ingredients::accumulated_dependency(explainer.train.xgb, 'ln_al')

 plot.ale <- ggplot_ale(ale.clay.xgb,ale.ph.xgb,
                        ale.som.xgb,ale.totalp.xgb)
 ggsave(plot = plot.ale,filename = 'products/lnk/withoutFeAl/0410plot2_ale.png',width = 13.16, height = 8.90, units='cm')
 
 # calculate R2 and RMSE for all models
 library(caret)
 #
 trainpred.xgb <- predict(model, newdata = dt.train.xgb)
 dt.train.xgb[, predicted := trainpred.xgb]
 defaultSummary(data.frame(obs = dt.train.xgb$ln_qmax, pred = trainpred.xgb))
 #
 testpred.xgb <- predict(model, newdata = dt.test.xgb)
 dt.test.xgb[, predicted := testpred.xgb]
 defaultSummary(data.frame(obs = dt.test.xgb$ln_qmax, pred = testpred.xgb))

 
 # plot 1 - 1 for both regression models
 
 # ##XGBOOST
 plotfun <- function(){
   #combine residuals1 and 2 in one dataframe
   train <- data.table(observed = dt.train.xgb$ln_qmax,
                       predicted = dt.train.xgb$predicted,
                       type = "training")
   test <- data.table(observed = dt.test.xgb$ln_qmax,
                      predicted = dt.test.xgb$predicted,
                      type = "test")
   combined <- rbind(train, test)
   #create label
   dt.label <- data.table(
     R2 = c(caret::R2(train$observed, train$predicted), caret::R2(test$observed, test$predicted)),
     RMSE = c(caret::RMSE(train$observed, train$predicted) ,caret::RMSE(test$observed, test$predicted)),
     MAE = c(caret::MAE(train$observed, train$predicted) ,caret::MAE(test$observed, test$predicted)),
     type = c("training", "test"))
   dt.label[, RMSE := round(RMSE, 2)]
   dt.label[, R2 := round(R2, 2)]#
   dt.label[, MAE := round(MAE, 2)]#
   dt.label[, label := paste("R2 =", R2, "\nRMSE = ", RMSE, "\nMAE = ", MAE)]
   #split
   train.label <- dt.label[type == "training"]
   test.label <- dt.label[type == "test"]
   #plot
   ggplot(combined, aes(x = observed, y = predicted, col = type)) + geom_point() +
     scale_color_manual(values  = c(`training` = "#e31a1c", `test` = "#2c7bb6")) +
     theme_bw() + labs(col = "Dataset") +
     #geom_smooth(method = 'lm', se = FALSE, linetype = "solid") +
     theme(axis.text = element_text(family="A",size = 12,colour ="black" ),#
           axis.title.x = element_text(family="A",size = 12,colour ="black",face = "bold"),#
           axis.title.y  = element_text(family="A",size = 12,colour ="black",face = "bold"))+#
     theme(legend.position = c(0.85,0.15),#
           legend.text = element_text(family="A", size = 9,face = "bold"),#
           legend.title = element_text(family="A",size = 9,face = "bold"),#
           legend.key =  element_rect(fill = "transparent", colour = NA))+#
     ggtitle("Log(Qmax)~pH+SOM+clay+totalp+soil type")+
     theme(plot.title = element_text(hjust = 0.5,family="A",size = 12 )) +
     xlim(3.5,8)+ylim(3.5,8)+
     geom_abline(intercept = 0,slope=1, colour='black',linetype=2, size=0.5)+
      geom_text(data = train.label, aes(x = 4, y = 7.8), label = train.label$label,
                inherit.aes = FALSE, col = "#e31a1c",family="A", size = 4,face = "bold") +
      geom_text(data = test.label, aes(x = 4, y = 7.1), label = test.label$label,
                inherit.aes = FALSE, col = "#2c7bb6",family="A", size = 4,face = "bold")
   }
 