
 # Convert the importance order results to a data table
 imp.xgb.dt <- data.table(Feature = imp.xgb$variable, 
                      Importance = imp.xgb$dropout_loss)
 imp.lm.dt <- data.table(Feature = imp.lm$variable, 
                          Importance = imp.lm$dropout_loss)
 
 
 # Set the order of the factor levels for plotting
 order_levels <- c("ln_som","ln_totalp","ln_clay","ph_h2o","ln_fe","ln_al","soilty_fluvisols",
                   "soilty_saline","soilty_alfisols","soilty_inceptisols",
                   "soilty_oxisols","soilty_mollisols")
 
 # Sort the data tables by importance in decreasing order
 imp.xgb.dt <- imp.xgb.dt[order(-Importance)]
 imp.lm.dt <- imp.lm.dt[order(-Importance)]
 
 # Calculate importance as a percentage of the total importance
 imp.xgb.dt$Importance <- imp.xgb.dt$Importance / sum(imp.xgb.dt$Importance) * 100
 imp.lm.dt$Importance <- imp.lm.dt$Importance / sum(imp.lm.dt$Importance) * 100
 
 # Plot the bar chart
 ggplot(imp.xgb.dt, aes(x = Importance, y = reorder(Feature, Importance))) +
   geom_bar(stat = "identity") +
   theme_bw() +
   labs(x = "Relative Importance (%)", y = "Feature") +
   ggtitle("XGBoost Feature Importance")
 
 ggplot(imp.lm.dt, aes(x = Importance, y = reorder(Feature, Importance))) +
   geom_bar(stat = "identity") +
   theme_bw() +
   labs(x = "Relative Importance (%)", y = "Feature") +
   ggtitle("Linear Regression Feature Importance")
 
 # make the ALE plots for different variables
 ggsave('products/lnk/withFeAl/plotVIP_xgb.png',width = 13.16, height = 10.90, units='cm')
 
 # try to make a plot for preicting Qmax with 2 important variables
 
 yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata))
 ALE_2vari <- ALEPlot::ALEPlot(dt.train.xgb[,-'ln_k'], X.model = explainer.train.xgb , 
                  pred.fun = yhat, J = c("ph_h2o","ln_totalp"), NA.plot = TRUE) 

 ggsave( 'products/qmax_kl/withFeAl/plotALE_2vari.png',width = 13.16, height = 10.90, units='cm')
 