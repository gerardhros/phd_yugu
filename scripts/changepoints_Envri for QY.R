 
#load packages
packages<- c("data.table", "ggplot2")
 lapply(packages, library, character.only = TRUE)
 # library(ggpubr);library(dplyr);library(Rmisc);library(lattice);library(plyr);library(gcookbook);library(rsq)
 # library(splines);library(drc);library(nlme)
 colors1 <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f")
 shape1 <- c(20,15,16,17,18,19)
 windowsFonts(A=windowsFont("Times New Roman"),
             B=windowsFont("Arial"))
 #load data
 dt <- fread("Origionaldata for Qiyang lte3.csv")
 cols_old <- colnames(dt)
 cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
 setnames(dt,cols_old,cols_new)
 dt[,som:= soc*2*0.1]
 dt[,bulk_density := 1000 / (0.625 + 0.025 * som + 0.0015 *43.86) ]
 ###Yi er al., 2016 Pedotransfer Functions for Estimating Soil Bulk Density: A Case Study in the Three-River Headwater Region of Qinghai Province, China
 #dt[,bulk_density := (1.398-0.0047*43.86-0.042*soc/10)*1000 ]
 #conversion factor from mg/kg to kg/ha
 dt[,conv_factor_kgha:= bulk_density*0.2*10000*10^(-6)]
 #2.25 is parameter from mg/kg to kg/ha and calculate different P pools(currently not use)
 #constants
 molarmassfe <- 55.845
 molarmassal <- 26.981539
 molarmassp <- 30.97376
 dt[, `p_ox (mmol/kg)` := `oxalatep`/molarmassp]
 dt[, `al_ox (mmol/kg)` := `oxalateal`/molarmassal]
 dt[, `fe_ox (mmol/kg)` := `oxalatefe`/molarmassfe]
 dt[,`olsen_p(mmol/kg)`:=`olsenp`/molarmassp]
 dt[,`total_p(mmol/kg)`:=`totalp`*1000/molarmassp]
 dt[,oxalatep_mmolkg:=oxalatep/molarmassp]#pox for measured data
 #calc qmax
 dt[, `qmax (mmol/kg)` := 0.5 * (`fe_ox (mmol/kg)` + `al_ox (mmol/kg)`)]###if can change into another value
 #assume maximum adsorption capacity is a constant
 dt[,  `FeAl (mmol/kg)` := mean((`al_ox (mmol/kg)`+`fe_ox (mmol/kg)`), na.rm = T)]
 dt[, `psi_ox_measure` := oxalatep_mmolkg / `FeAl (mmol/kg)` ]#measurement

 #prepare data regression
 #delete NA values for oxalatep
 #delete the rows with cacl2p>5.2 i assume that point was outlier
 d1 <- dt[treatments%in%c("CK","NP","NPK","NPKM","HNPKM","M")][oxalatep!='NA'][year!='1990'][cacl2p <=5.2]
 ##regression for CacL2 vs olsen
 lm1 <- lm(cacl2p~olsenp,data=d1)
 os1 <- segmented(lm1,~olsenp, psi = 50)
 #davies.test(os,"psi_ox_measure",k=1)
 slope(os1)
 plot(os1)
 summary(os1)

  fig1 <- ggplot(d1,aes(x=olsenp,y=cacl2p))+
  geom_point(aes(shape= treatments,color= treatments), 
             dt[treatments%in%c("CK","NP","NPK","NPKM","HNPKM","M")][oxalatep!='NA'][year!='1990'],size=2, color="#e31a1c")+
  geom_smooth(method="lm",aes(x=olsenp,y=cacl2p),
              formula=y ~ bs(x,degree=1,knots=c(47)),
              data=d1,
              se=F,color='black')+
  scale_shape_manual(values = shape1,name="Treatments", 
                     breaks = c("CK","NP","NPK","NPKM","M","HNPKM"),
                     labels=c("CK","NP","NPK","NPKM","M","HNPKM"))+
  geom_point(aes(x=olsenp,y=cacl2p,shape=treatments),color='#d3d3d3',size=2,
             data=dt[treatments%in%c("CK","NP","NPK","NPKM","HNPKM","M")][oxalatep!='NA'][cacl2p > 5.2])+
  stat_cor(aes(x=olsenp,y=cacl2p, label =c( ..rr.label..)),
           method = 'pearson', label.x = 0.01, label.y = 6,family="A", size =4,colour="black")+
  scale_x_continuous(breaks = seq(0,250,25),limits=c(0,250))+
  scale_y_continuous(breaks = seq(0,6,1),limits=c(0,6))+
  labs(x=expression(Olsen~P~(mg~kg^-1)),y=expression(CaCl[2]~P~(mg~kg^-1)))+
  theme(axis.text = element_text(family="A",size = 9,colour ="black" ),
        axis.title.x = element_text(family="A",size = 10,colour ="black",face = "bold"),
        axis.title.y  = element_text(family="A",size = 10,colour ="black",face = "bold"))+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(fill = NA))+
  theme(legend.key=element_blank(),legend.position = "none")

##regression for CaCl2 P vs PSI and plot
  lm2 <- lm(cacl2p~`psi_ox_measure`,data=d1)
  os2 <- segmented(lm2,~`psi_ox_measure`, psi = 0.123)
  #davies.test(os,"psi_ox_measure",k=1)
  slope(os2)
  plot(os2)
  summary(os2)##change point 0.141, but both coefficients are insignificant 

  fig2 <- ggplot(d1, aes(x=`psi_ox_measure`,y=cacl2p))+
  geom_point(aes(shape= treatments,color= treatments), data=dt[treatments%in%c("CK","NP","NPK","NPKM","HNPKM","M")][oxalatep!='NA'][year!='1990'],size=2, color="#e31a1c")+
  geom_smooth(method="lm",aes(x=`psi_ox_measure`,y=cacl2p),
              formula=y ~ bs(x,degree=1,knots=c(0.141)),
              data=d1,
              se=F,color='black')+
  scale_shape_manual(values = shape1,name="Treatments", 
                     breaks = c("CK","NP","NPK","NPKM","M","HNPKM"),
                     labels=c("CK","NP","NPK","NPKM","M","HNPKM"))+
  geom_point(aes(x=`psi_ox_measure`,y=cacl2p,shape=treatments),color='#d3d3d3',size=2,
             data=dt[treatments%in%c("CK","NP","NPK","NPKM","HNPKM","M")][oxalatep!='NA'][cacl2p > 5.2])+
  stat_cor(aes(x=`psi_ox_measure`,y=cacl2p, label =c( ..rr.label..)),
           method = 'pearson', label.x = 0.01, label.y = 6,family="A", size =4,colour="black")+
  scale_x_continuous(breaks = seq(0,0.4,0.1),limits=c(0,0.4))+
  scale_y_continuous(breaks = seq(0,6,1),limits=c(0,6))+
  labs(x=expression(P~Saturation~Index~(PSI)),y=expression(CaCl[2]~P~(mg~kg^-1)))+
  theme(axis.text = element_text(family="A",size = 9,colour ="black" ),
        axis.title.x = element_text(family="A",size = 10,colour ="black",face = "bold"),
        axis.title.y  = element_text(family="A",size = 10,colour ="black",face = "bold"))+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(fill = NA))+
  theme(legend.key=element_blank(),legend.position = "bottom")

  fig4 <- ggarrange( fig2,fig1,ncol = 2,nrow = 1)
  ggsave(filename = "C:/Users/gu021/OneDrive - Wageningen University & Research/1 WUR_work/2 longtermExperi/plots_p2/cacl2-soilp-2.png",
         units = "in",width=6.91,height = 5.09, dpi = 3700) 

