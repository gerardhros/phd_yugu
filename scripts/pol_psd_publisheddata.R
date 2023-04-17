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

  d1$class2 <- ifelse(d1$mpsc <= 10, "<10", ifelse(d1$mpsc <= 30, "10-30", 
                                                   ifelse(d1$mpsc <= 50, "30-50",
                                                          ifelse(d1$mpsc <= 80, "50-80", 
                                                                 ifelse(d1$mpsc<= 150,"80-150",
                                                                        ifelse(d1$mpsc<= 300,"150-300", ">300"))))))
  d1$class2 <- as.factor(factor(d1$class2,levels = c("<10","10-30","30-50","50-80","80-150","150-300",">300")))
  
  
  color <- c("<10"="#ffffd4", "10-30"="#fee391", "30-50"= "#fec44f",
             "50-80"="#fe9929","80-150"="#ec7014", "150-300"="#cc4c02", ">300" = "#8c2d04" )
  ggplot(d1[`psd`<=100][`pols`<=600][`author`!="MURPHY"][note!="caco3[34%]" & note!="caco3[14%]" & note!="caco3[24%]"][`mpsc`<=300],
         aes(x=`pols`,y=psd, fill=class2))+
    geom_point(aes(color=class2))+
    geom_smooth(method = 'lm', se = FALSE,aes(linetype=class2, colour=class2) )+
    
    scale_color_manual(values =color )+
    labs( x = "P Olsen [mg/kg]", y = "PSD[%]")+
    guides(fill=guide_legend(title = "MaxPSC[mmol/kg]"))+
    stat_poly_eq(use_label(c("eq", "R2"))) +
    theme(legend.position = c(0.85,0.85))+
    theme_bw()
  
  #test for the multilinear regression between PSD & polsen plus mpsc
  summary(lm(formula = psd ~ pols * mpsc + I(pols^2) + I(mpsc^2),
             data = d1[psd <= 100][pols <= 600][author != "MURPHY"][note != "caco3[34%]" &
                                                                      note!="caco3[14%]" & note!="caco3[24%]"][`mpsc`<=300]))
                                                                       
  ggsave("dev/psd_pol.png",width = 18, height = 13.16, units='cm', dpi = 800)
  
  # plot basic soil properies
  d2 <- d1[psd <= 100][pols <= 600][author != "MURPHY"][note != "caco3[34%]" &
                                                          note!="caco3[14%]" & note!="caco3[24%]"][`mpsc`<=300]
 
  ggplot(d2)+
    geom_histogram(aes(pox))+
    theme_bw()+
    xlab("P Oxalate[mmol/kg]")+ylab("Count")
  
####### predict PSD based on POLSEN and fertilization dataset
  
  ## load fertilization dataset 
  fert <- fread("data/experimentalData/fertilizationData_cn.csv")
  
  #select related variables, x, y, olsenP and SOM
  
  fert_cn <- fert[,c("latitute","longtitude","som","pols","province")]
 
  
  ## select data for Fe and Al for each province from published dataset for paper1
  
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
  
  d2 <- d1[`fe`!= " " & `al`!=" "]
  
  d2[, feox_av:= mean(fe),by= 'province']
  d2[, alox_av:= mean(al),by= 'province']
  d3 <- d2[,c("lontitudee","latituden","province","alox_av","feox_av","alox_av","id","author","title")]
  #slect one value with related row forming into a new tataable
  d4 <- d3[, .(lontitudee = lontitudee[1], 
               latituden = latituden[1], 
               province = province[1], 
               alox_av = alox_av[1], 
               feox_av = feox_av[1], 
               id = id[1], author = author[1], 
               title = title[1]), by = `province`][,mpsc:= 0.5*(alox_av+feox_av)]
 ## add feox and alox from published dataset to fertilization dataset
  fert_cn$feox_av <- d4$feox_av[match(fert_cn$province,d4$province)]
  fert_cn$alox_av <- d4$alox_av[match(fert_cn$province,d4$province)]

  fert_cn[, psd_meas:= -0.002*pols^2+0.001*((feox_av+alox_av))+0.59*pols-0.27*(feox_av+alox_av)]
  
  ggplot(fert_cn[psd_meas > 0,])+
    geom_histogram(aes(psd_meas))+
    xlab("Predicted PSD[%]")+
    ylab("Number of points")+
    ggtitle("Predicted PSD by Pols and MPSC")+
    theme_bw()
  ggsave("dev/psd_predicted.png",width = 13.16, height = 13.16, units='cm', dpi = 800)
  