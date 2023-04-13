# clear environment
  rm(list = ls())
#########################XGBOOST###############################
  library(data.table);library(ggplot2);library(ggpmisc)

  d1 <- fread("data/experimentalData/POLSEN_psd.csv")
  # Rename columns
  cols_old <- colnames(d1)
  cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
  setnames(d1,cols_old,cols_new)
  
  d1[,psd:= (`pox`/(0.5*(`feox`+`alox`))*100)]
  d1[,mpsc:=0.5*(`alox`+`alox`) ]

  # different classifications for the range of MPSC(maximum P sorption capacity)
  d1$class2 <- ifelse(d1$mpsc <= 30, "<30", ifelse(d1$mpsc <= 50, "30-50", 
                                                    ifelse(d1$mpsc <= 100, "50-100",
                                                           ifelse(d1$mpsc <= 200, "100-200", ">200"
                                                                  ))))
  d1$class2 <- factor(d1$class2,levels = c("<30","30-50","50-100","100-200",">200"))
 
  ggplot(d1[`psd`<=100][`pols`<=600],aes(x=`pols`,y=psd))+
    geom_point(aes(fill=mpsc))+
    scale_color_gradient(low = "#ffffb2",high = "#b10026")+
   scale_fill_manual(values  = c('<30'="#2c7bb6", '30-50'="#abd9e9", '50-100'= "#ffffbf",'100-200'="#fdae61",'>200'="#d7191c"))+
   labs( x = "P Olsen [mg/kg]", y = "PSD[%]")+
   guides(fill=guide_legend(title = "MaxPSC /n[mmol/kg]"))+
   #stat_poly_line(aes(fill= class2,linetype=class2),se = FALSE) +
   geom_smooth(method = 'lm', se = FALSE,aes(fill= class2,linetype=class2) )+
   stat_poly_eq(use_label(c("eq", "R2"))) +
   theme_bw()
  
 ggsave("dev/pols_psd.png")

 