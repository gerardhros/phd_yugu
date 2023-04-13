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

  d1$class2 <- ifelse(d1$mpsc <= 30, "<30", ifelse(d1$mpsc <= 50, "30-50", 
                                                   ifelse(d1$mpsc <= 100, "50-100",
                                                          ifelse(d1$mpsc <= 200, "100-200", ">200"
                                                          ))))
  d1$class2 <- as.factor(factor(d1$class2,levels = c("<30","30-50","50-100","100-200",">200")))
  
  
  color <- c("<30"="#fef0d9", "30-50"="#fdcc8a", "50-100"= "#fc8d59","100-200"="#e34a33",">200"="#b30000")
  ggplot(d1[`psd`<=100][`pols`<=600][`author`!="MURPHY"],aes(x=`pols`,y=psd, fill=class2))+
    geom_point(aes(color=class2))+
    geom_smooth(method = 'lm', se = FALSE,aes(linetype=class2, colour=class2) )+
    
    scale_color_manual(values =color )+
    labs( x = "P Olsen [mg/kg]", y = "PSD[%]")+
    guides(fill=guide_legend(title = "MaxPSC[mmol/kg]"))+
    stat_poly_eq(use_label(c("eq", "R2"))) +
    theme(legend.position = c(0.85,0.85))+
    theme_bw()
  
  ggsave("dev/psd_pol.png",width = 18, height = 13.16, units='cm', dpi = 800)

 