#------------calculation for PSD
# windowsFonts(A=windowsFont("Times New Roman"),
#              B=windowsFont("Arial"))
#------
library(data.table);library(ggplot2);library(dplyr)
#load data from Nawara
ry_eu1<- fread("data/experimentalData/20220906Nawara_data.csv")
cols_old <- colnames(ry_eu1)
cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
setnames(ry_eu1,cols_old,cols_new)
ry_eu1[, `p_ols(mmol/kg)`:= p_olsen/30.973762]
ry_eu1[,psd:= `p-ox (mmol/kg)`/(0.5*(`al-ox (mmol/kg)` +`fe-ox (mmol/kg)`))]
ry_eu1[,psd3:= `p-ox (mmol/kg)`/(0.8*(`al-ox (mmol/kg)` +`fe-ox (mmol/kg)`))]
ry_eu1[,psol:= `p_ols(mmol/kg)`/(0.5*(`al-ox (mmol/kg)` +`fe-ox (mmol/kg)`))]
#ry_eu1[,ratio:= `p_olsen`/`fe-ox (mmol/kg)`]

ry_eu <- ry_eu1[gewassen=='wheat'|gewassen=='maize']#avoid change all the names in the following steps
ry_eu[, max_yield_year:= max(yield), by =c("year","site") ]
ry_eu[, relative_yield_eu:= yield/max_yield_year*100]

#-----
#load data for Qiyang
dt <- fread("data/experimentalData/Origionaldata for Qiyang lte3.csv")
cols_old <- colnames(dt)
cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
setnames(dt,cols_old,cols_new)
dt[,som:= soc*2*0.1]
dt[,bulk_density := 1000 / (0.625 + 0.025 * som + 0.0015 *43.86) ]
#constants
molarmassfe <- 55.845
molarmassal <- 26.981539
molarmassp <- 30.973762
#calculation
dt[, `p_ox(mmol/kg)` := `oxalatep`/molarmassp]
dt[, `al_ox(mmol/kg)` := `oxalateal`/molarmassal]
dt[, `fe_ox(mmol/kg)` := `oxalatefe`/molarmassfe]
dt[, `p_ols(mmol/kg)` := `olsenp`/molarmassp]
dt[,`psi`:= `p_ox(mmol/kg)`/(`al_ox(mmol/kg)`+`fe_ox(mmol/kg)`)]
dt[,`psd`:= `p_ox(mmol/kg)`/(0.5 * (`al_ox(mmol/kg)`+`fe_ox(mmol/kg)`))]
dt[,`psd2`:= `p_ox(mmol/kg)`/(0.75 * (`al_ox(mmol/kg)`+`fe_ox(mmol/kg)`))]
dt[,`psd3`:= `p_ox(mmol/kg)`/(0.8 * (`al_ox(mmol/kg)`+`fe_ox(mmol/kg)`))]
#calculate relative yield for Qiyang
#selected treatments and delete NA values
dt_1 <- dt[treatments %in% 
             c("CK","NP","NPK","NPKM","HNPKM","M"),] [maizegrain!= 'NA',][wheatgrain!= 'NA',]
dt_1[, max_maize_year:= max(maizegrain), by= year]
dt_1[, max_wheat_year:= max(wheatgrain), by= year]
dt_1[, ry_maize_year:= maizegrain/max_maize_year]
dt_1[, ry_wheat_year:= wheatgrain/max_wheat_year]
dt_1[, ry_mean_year:= (ry_maize_year+ry_wheat_year)/2*100]#mean ry for maize and wheat by year
ry_qy <- filter(dt_1, cacl2p!= 'NA') 

##regression for plot
#--------
library(segmented);library(caret)
#environmental risk
##olsenp
#eu
l1_eu <- lm(p_cacl2~p_olsen,data=ry_eu)
os1_eu <- segmented(l1_eu,~p_olsen,psi = 30)#meothod 1
os11_eu <- glm(p_cacl2 ~ p_olsen + I(pmax(0,p_olsen - 39.2)), data=ry_eu)#method 2
#qy
l1_qy <- lm(cacl2p~olsenp,data=ry_qy)#delete 6 9 18 individually or combined
os1_qy <- segmented(l1_qy,~olsenp,psi =15) #meothod 1
os11_qy <- glm(cacl2p ~ olsenp + I(pmax(0,olsenp - 23.2)), data=ry_qy)#method 2

##PSD
#eu psd
l3_eu <- lm(p_cacl2~ psd,data=ry_eu)
os3_eu <- segmented(l3_eu,~ psd,psi = 0.1)#meothod 1
os33_eu <- glm(p_cacl2 ~ psd + I(pmax(0,psd - 0.33)), data=ry_eu)#method 2

l3_eu3 <- lm(p_cacl2~ psd3,data=ry_eu)
os3_eu3 <- segmented(l3_eu3,~ psd3,psi = 0.1)#meothod 1
os33_eu3 <- glm(p_cacl2 ~ psd3 + I(pmax(0,psd3 - 0.21)), data=ry_eu)#method 2

#qy
l3_qy1 <- lm(cacl2p~psd,data=ry_qy)
os3_qy1 <- segmented(l3_qy1,~psd, psi = 0.1)#meothod 1
os33_qy1 <- glm(cacl2p ~ psd + I(pmax(0,psd - 0.32)), data=ry_qy)#method 2

l3_qy2 <- lm(cacl2p~psd2,data=ry_qy)
os3_qy2 <- segmented(l3_qy2,~psd2, psi = 0.1)#meothod 1
os33_qy2 <- glm(cacl2p ~ psd2 + I(pmax(0,psd2 - 0.21)), data=ry_qy)#method 2

l3_qy3 <- lm(cacl2p~psd3,data=ry_qy)
os3_qy3 <- segmented(l3_qy3,~psd3, psi =0.20)#meothod 1
os33_qy3<- glm(cacl2p ~ psd3 + I(pmax(0,psd3 - 0.20)), data=ry_qy)#method 2

caret::R2(ry_eu$p_cacl2, predict(os11_eu,data = ry_eu))
caret::R2(ry_qy$cacl2p, predict(os33_qy3,data = ry_qy))
caret::R2(ry_qy$cacl2p, predict(os33_qy1,data = ry_qy))
caret::R2(predict(os11_eu,data = ry_eu),ry_eu$p_cacl2)
caret::RMSE(predict(os33_qy3,data = ry_qy),ry_qy$cacl2p)



##relative yield
ry_olsenp_qy <-  nls(ry_mean_year~SSasymp(olsenp,Asym,Rs0,lrc), 
                     data = ry_qy[ph>5][year!=1990])
ry_psd_qy <-  nls(ry_mean_year~SSasymp(psd,Asym,R0,lrc), 
                   data = ry_qy[ph>5][year!=1990])
ry_psd3_qy <-  nls(ry_mean_year~SSasymp(psd3,Asym,R0,lrc), 
                  data = ry_qy[ph>5][year!=1990])

ry_psd_eu <-  nls(relative_yield_eu~SSasymp(psd,Asym,R0,lrc), data = ry_eu)
ry_psd3_eu <-  nls(relative_yield_eu~SSasymp(psd3,Asym,R0,lrc), data = ry_eu)
ry_olsenp_eu <-  nls(relative_yield_eu~SSasymp(p_olsen,Asym,R0,lrc),  data = ry_eu)
ry_oxalatep_eu <-  nls(relative_yield_eu~SSasymp(p_ox,Asym,R0,lrc), data = ry_eu)


# standard error of prediction
library(caret)
caret::R2(ry_qy[ph>5][year!=1990]$ry_mean_year, predict(ry_olsenp_qy,data = ry_qy[ph>5][year!=1990]))
caret::R2(ry_qy[ph>5][year!=1990]$ry_mean_year, predict(ry_psd3_qy,data = ry_qy[ph>5][year!=1990]))
caret::R2(ry_qy[ph>5][year!=1990]$ry_mean_year, predict(ry_psd_qy,data = ry_qy[ph>5][year!=1990]))
caret::R2(ry_eu$relative_yield, predict(ry_olsenp_eu,data = ry_eu))
caret::RMSE( predict(ry_olsenp_qy,data = ry_qy[ph>5][year!=1990]),ry_qy[ph>5][year!=1990]$ry_mean_year)
caret::RMSE( predict(ry_psd_qy,data = ry_qy[ph>5][year!=1990]),ry_qy[ph>5][year!=1990]$ry_mean_year)
caret::RMSE(predict(ry_olsenp_eu,data = ry_eu),ry_eu$relative_yield)
caret::RMSE(predict(ry_psd_eu,data = ry_eu),ry_eu$relative_yield)

caret::RMSE(predict(ry_psd3_qy,data = ry_qy[ph>5][year!=1990]),ry_qy[ph>5][year!=1990]$ry_mean_year)

dx = data.table(place=c('olsen_qy','psd_qy','psd3_qy',
                       'olsen_eu','psd_eu','psd_eu3'),
                Asym=c( summary(ry_olsenp_qy)$coefficients[1,1],
                        summary(ry_psd_qy)$coefficients[1,1],
                        summary(ry_psd3_qy)$coefficients[1,1],
                        summary(ry_olsenp_eu)$coefficients[1,1],
                       summary(ry_psd_eu)$coefficients[1,1],
                       summary(ry_psd3_eu)$coefficients[1,1]),
                R0= c( summary(ry_olsenp_qy)$coefficients[2,1],
                     summary(ry_psd_qy)$coefficients[2,1],
                      summary(ry_psd3_qy)$coefficients[2,1],
                     summary(ry_olsenp_eu)$coefficients[2,1],
                      summary(ry_psd_eu)$coefficients[2,1],
                     summary(ry_psd3_eu)$coefficients[2,1]),
                lrc= c(summary(ry_olsenp_qy)$coefficients[3,1],
                      summary(ry_psd_qy)$coefficients[3,1],
                       summary(ry_psd3_qy)$coefficients[3,1],
                       summary(ry_olsenp_eu)$coefficients[3,1],
                       summary(ry_psd_eu)$coefficients[3,1],
                      summary(ry_psd3_eu)$coefficients[3,1]))

dx[, cp_80:= log((Asym*0.80-Asym)/(R0-Asym))/(-exp(lrc))]
dx[, cp_85:= log((Asym*0.85-Asym)/(R0-Asym))/(-exp(lrc))]
dx[, cp_90:= log((Asym*0.90-Asym)/(R0-Asym))/(-exp(lrc))]
dx[, cp_95:= log((Asym*0.95-Asym)/(R0-Asym))/(-exp(lrc))]
#environemntal thresholds
##add thresholds for all soil P indicators to dataframe dx (except data combination)
dx[, thre_envri:= c(NA,os3_qy1$psi[2],os3_qy3$psi[2],
                    os1_eu$psi[2] ,os3_eu$psi[2],
                    os3_eu3$psi[2])]
dx[,thre_cacl2:= c(NA,  predict(os33_qy1, newdata = data.frame(psd =c(os3_qy1$psi[2]))),
                   predict(os33_qy3, newdata = data.frame(psd3=c(os3_qy3$psi[2]))),
                   predict(os11_eu, newdata = data.frame(p_olsen=c(os1_eu$psi[2]))),
                   predict(os33_eu, newdata = data.frame(psd=c(os3_eu$psi[2]))),
                   predict(os33_eu3, newdata = data.frame(psd3=c(os3_eu3$psi[2]))))]

#plot with environmental risk
library(splines)
#Envri_psd <- 
ggplot()+
  geom_point(aes(x=psd,y=p_cacl2),size=2,colour='#1f78b4',data=ry_eu)+
  geom_point(aes(x=psd,y=cacl2p),size=2,colour='#d7191c',data=ry_qy )+#0.5
  # geom_point(aes(x=psd2,y=cacl2p),color= '#d7191c', size=2,
  #            data = ry_qy)+#
  geom_point(aes(x=psd3,y=cacl2p),color= '#fdae61', size=2,
             data = ry_qy)+#alpha = 0.75
  geom_smooth(method="lm",aes(x=psd,y=p_cacl2),
              formula=y ~ bs(x,degree=1,knots=c(0.33)),
              data=ry_eu,se=F,color='#1f78b4')+
   geom_smooth(method="lm",aes(x=psd,y=cacl2p),
               formula=y ~ bs(x,degree=1,knots=c(0.32)),
              data=ry_qy,se=F,color='#d7191c')+
  geom_smooth(method="lm",aes(x=psd3,y=cacl2p),
               formula=y ~ bs(x,degree=1,knots=c(0.19)),
              data=ry_qy,se=F,color='#fdae61')+
  # geom_smooth(method="lm",aes(x=psd3,y=cacl2p),
  #             formula=y ~ bs(x,degree=1,knots=c(0.28)),
  #             data=ry_qy,se=F,color='#99d594')+
  scale_y_continuous(breaks = seq(0,20,4),limits=c(0,20))+
  theme(legend.text = element_text(family="Times", size = 14))+
  theme(legend.title = element_text(family="Times",size = 14))+
  guides(col=guide_legend(ncol = 1))+
  theme(legend.key = element_blank())+
  labs(x=expression(P~Saturation~Degree~(PSD)),y=expression(CaCl[2]~P~(mg~kg^-1)))+
  theme(axis.text = element_text(family="A",size = 9,colour ="black" ),
        axis.title.x = element_text(family="A",size = 10,colour ="black",face = "bold"),
        axis.title.y  = element_text(family="A",size = 10,colour ="black",face = "bold"))+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(fill = NA))+
  theme(legend.position = "none")


RY_psd <- 
ggplot()+
  geom_point(aes(x=psd,y=ry_mean_year),size=2,colour='#d7191c', data= ry_qy[ph>5])+#1
  geom_point(aes(x=psd3,y=ry_mean_year),size=2,colour='#fdae61', data= ry_qy[ph>5])+# for alpha = 0.75
  geom_point(aes(x=psd,y=ry_mean_year),size=2,colour='#d3d3d3', data= ry_qy[ph<=5])+
  geom_point(aes(x=psd,y=relative_yield_eu),size=2,colour='#1f78b4',data = ry_eu)+
  geom_smooth(aes(x=psd,y=relative_yield_eu),method="nls", 
              formula =y~SSasymp(x,Asym,R0,lrc),color="#1f78b4",se=F,
              data = ry_eu) +
  geom_smooth(aes(x=psd,y=ry_mean_year),method="nls", 
              formula =y~SSasymp(x,Asym,R0,lrc),color="#d7191c",se=F,
              data = ry_qy[ph>5]) +
  geom_smooth(aes(x=psd3,y=ry_mean_year),method="nls", 
               formula =y~SSasymp(x,Asym,R0,lrc),color="#fdae61",se=F,
               data = ry_qy[ph>5]) +
  theme(legend.text = element_text(family="Times", size = 14))+
  theme(legend.title = element_text(family="Times",size = 14))+
  guides(col=guide_legend(ncol = 1))+
  theme(legend.key = element_blank())+
  labs(x=expression(P~Saturation~Degree~(PSD)),y=expression(Relative~Crop~Yield))+
  theme(axis.text = element_text(family="A",size = 9,colour ="black" ),
        axis.title.x = element_text(family="A",size = 10,colour ="black",face = "bold"),
        axis.title.y  = element_text(family="A",size = 10,colour ="black",face = "bold"))+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(fill = NA))+
  geom_hline(yintercept = 100,colour='red',linetype=2, size=0.5)


RY_polsen <-
 ggplot()+
  geom_point(aes(x=olsenp,y=ry_mean_year),size=2,colour='#e31a1c', data= ry_qy)+
  geom_point(aes(x=olsenp,y=ry_mean_year),size=2,colour='#d3d3d3', data= ry_qy[ph<=5])+
  geom_point(aes(x=p_olsen,y=relative_yield_eu),size=2,colour='#1f78b4',data = ry_eu)+
  geom_smooth(aes(x=p_olsen,y=relative_yield_eu),method="nls", 
              formula =y~SSasymp(x,Asym,R0,lrc),color="#1f78b4",se=F,
              data = ry_eu) +
  geom_smooth(aes(x=olsenp,y=ry_mean_year),method="nls", 
              formula =y~SSasymp(x,Asym,R0,lrc),color="#e31a1c",se=F,
              data =ry_qy[ph>5][year!=1990]) +
  theme(legend.text = element_text(family="Times", size = 14))+
  theme(legend.title = element_text(family="Times",size = 14))+
  guides(col=guide_legend(ncol = 1))+
  theme(legend.key = element_blank())+
  labs(x=expression(Olsen~P~(mg~kg^-1)),y=expression(Relative~crop~yield))+
  theme(axis.text = element_text(family="A",size = 9,colour ="black" ),
        axis.title.x = element_text(family="A",size = 10,colour ="black",face = "bold"),
        axis.title.y  = element_text(family="A",size = 10,colour ="black",face = "bold"))+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(fill = NA))+
  geom_hline(yintercept = 100,colour='red',linetype=2, size=0.5)

Envri_polsen <- 
ggplot()+
   geom_point(aes(x=p_olsen,y=p_cacl2),size=2,colour='#1f78b4',data= ry_eu)+
   geom_point(aes(x=olsenp,y=cacl2p),color= '#e31a1c',size=2,
              data = ry_qy)+
   # geom_point(aes(x=olsenp,y=cacl2p),color= '#d3d3d3',size=2,
   #           data = ry_qy[c(6,9,18),])+ 
   geom_smooth(method="lm",aes(x=p_olsen,y=p_cacl2),
               formula=y ~ bs(x,degree=1,knots=c(39.22)),
               data=ry_eu,se=F,color='#1f78b4')+
   # geom_smooth(method="lm",aes(x=olsenp,y=cacl2p),
   #             formula=y ~ bs(x,degree=1,knots=c(23.2)),
   #             data=ry_qy[-c(6,9,18),],se=F,color='#e31a1c')+
   theme(legend.position = c(0.01,5))+
   scale_y_continuous(breaks = seq(0,20,4),limits=c(0,20))+
   theme(legend.text = element_text(family="Times", size = 14))+
   theme(legend.title = element_text(family="Times",size = 14))+
   guides(col=guide_legend(ncol = 1))+
   theme(legend.key = element_blank())+
   labs(x=expression(Olsen~P~(mg~kg^-1)),y=expression(CaCl[2]~P~(mg~kg^-1)))+
   theme(axis.text = element_text(family="A",size = 9,colour ="black" ),
         axis.title.x = element_text(family="A",size = 10,colour ="black",face = "bold"),
         axis.title.y  = element_text(family="A",size = 10,colour ="black",face = "bold"))+
   theme(panel.grid.major = element_line(colour = NA),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA),
         panel.border = element_rect(fill = NA))+
   theme(legend.position = "none")


 
##spatial analysis
#------
# library(ggpubr) 
# ggarrange(Envri_psd,RY_psd,ncol = 2,nrow = 1) 
 
 #ggsave(filename = "F:/Onedrive/OneDrive - Wageningen University & Research/1 WUR_work/2 PSI paper/plots/indi_psd.png",units = "in", dpi = 800)
##spatial analysis
#------

 library(sf);library(raster);library(terra);library(gstat);library(stars)
 ##load data for Qiyang shape
 QY <- st_read("data/spatialInfor/Qiyang.shp") 
 QY_map <- QY[,c(2,5)]
 QY_sf <-st_transform(QY_map, "epsg:32649")#
 sum(st_area(QY_sf))#area
 st_crs(QY_sf)#get/set Coordinate reference system
 QYBorder <- st_cast(QY_sf, "MULTILINESTRING")# 
 
 # Land use raster for idw
 QYtif <- terra::rast('data/spatialInfor/ld2015.tif') 
 QYtif <- terra::project(QYtif, "epsg:32649")
 
 
 reclass_table <- matrix(data=cbind(  c(1,11,18,26,34,46),
                                      c(11,18,26,34,46,100),
                                      c(1,2,3,3,4,5)), ncol=3)# evidence for re_group 5 groups
 # 1 = paddy
 # 2 = upland
 # 3 = Forest + grassland
 # 4 = waterfield
 # 5 = urban and infrastructure + nature

 QY_recl <- classify(QYtif, reclass_table, filename="data/spatialInfor/QY_reclass.tif",
                     overwrite=T)
 #frequency distribution for different types 
 freq(QY_recl)
 QY_rast <- st_as_stars(QY_recl) 
 
 
 ##QY_data is origional dataset, and QY_data3 is dataset after finish outlier test
 QY_data <- fread("data/experimentalData/qiyang-inventory_coodinate.csv")
 cols_old <- colnames(QY_data)
 cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))# doesn't work for read.csv
 setnames(QY_data,cols_old,cols_new)
 QY_data[, psd:= (oxalatep/31)/(0.5*((oxalateal/27)+(oxalatefe/56)))]
 QY_data[, psd3:= ((oxalatep/31)/(0.8*((oxalateal/27)+(oxalatefe/56))))]

 ##outlier test
 QY_data2 <- QY_data[, .(lat,lon,landuse_code,cacl2p,olsenp,oxalatep,oxalatefe,oxalateal,totalp,som,psd,psd3,ph,claycontent)]
 QY_data3 <- QY_data2[-c(176,143,111,112,80),][oxalatep<491]##why delete oxalatep>491 
 
 qy <- QY_data3[,c("lat","lon","landuse_code",
                   "claycontent","som","ph","cacl2p","olsenp","oxalatep",
                   "oxalatefe","oxalateal","totalp","psd","psd3")]
 qy <- st_as_sf(qy, coords= c("lat","lon"), crs = 4214)
 st_crs(qy)
 qy_sf<-st_transform(qy , "epsg:32649")

# spatial interpolation for psd regional alpha = 0.5
 idw.weighing <- lapply(X = 1:3, FUN = function(x){
   #inverse distance weigh
   distanceweigh <-idw(formula = psd ~ 1,#key parameters
                       locations = qy_sf[qy_sf$landuse_code == x,],
                       newdata =QY_rast ,
                       idp = 2)
   #form stars to terra
   distanceweigh <- terra::rast(distanceweigh)
   distanceweigh$var1.var_lyr.1 <- NULL
   names(distanceweigh) <- "psd"
   #resample to landuse
   distanceweigh <- terra::resample(distanceweigh, terra::rast(QY_rast))
 })
 names(idw.weighing) <- c("paddy", "upland", "other")

 
 #bind land use map to list
 idw.weighing[["landuse"]] <- terra::rast(QY_rast)
 #stack
 idw.stacked <- terra::rast(idw.weighing)
 #to data.table
 idw.stacked.dt <- terra::as.data.frame(idw.stacked, xy = TRUE) |> as.data.table()
 idw.stacked.dt[landuse == 1, psd := paddy]
 idw.stacked.dt[landuse == 2, psd := upland]
 #we remove show NA for forests and nautral areas since we are only interested in agriculture
 idw.stacked.dt[landuse == 3, psd := NA]
 #to raster
 idw.final <- terra::rast(idw.stacked.dt)
 
#visualisation steps
 visualisation.dt <- terra::as.data.frame(idw.final,xy =TRUE)
 visualisation.dt <- as.data.table(visualisation.dt)
 #general spatial distribution
 visualisation.dt[, type_qy:= ifelse(psd< 0.1, '< 0.1',ifelse(psd< 0.15,'0.1-0.15',ifelse(psd<0.2,'0.15-0.2','> 0.2')))]
 
 ##alpha=0.5 for regional PSD estimation and 0.8 for TR-agri-QY BASED
 visualisation.dt[, type_qy:= ifelse(psd< as.numeric(dx[3,5]), 'Low',ifelse(psd<as.numeric(dx[3,7]),'Moderate','Sufficient'))]#qy_psd3 (only change CV)
 ##alpha=0.5 for regional PSD estimation and 0.8 for TR-envri-QY BASED
 visualisation.dt[, type_qy:= ifelse(psd< (as.numeric(dx[3,9])*0.5), 'Low risk',ifelse(psd<as.numeric(dx[3,9]),'Moderate risk','High risk'))]#qy_psd3 0.8 (only change CV)
 ##alpha=0.5 for regional PSD estimation and 0.5 for TR-agri-EU BASED
 visualisation.dt[, type_qy:= ifelse(psd< as.numeric(dx[5,5]), 'Low',ifelse(psd<as.numeric(dx[5,7]),'Moderate','Sufficient'))]#eu psd
 ##alpha=0.5 for regional PSD estimation and 0.5 for TR-envri-EU BASED
 visualisation.dt[, type_qy:= ifelse(psd<(as.numeric(dx[5,9])*0.5), 'Low risk',ifelse(psd<as.numeric(dx[5,9]),'Moderate risk','High risk'))]#eu
 
## calculate the percentage of area where soil P is lower than thressholds
 farea1  = nrow(visualisation.dt[type_qy=='Low'])*1000*1000 # area in m2
 farea2  = nrow(visualisation.dt[type_qy=='Moderate'])*1000*1000 # area in m2
 farea3  = nrow(visualisation.dt[type_qy=='Sufficient'])*1000*1000 # area in m2
 farea1/(farea1+farea2+farea3)
 farea2/(farea1+farea2+farea3)
 farea3/(farea1+farea2+farea3)
 
 farea4  = nrow(visualisation.dt[type_qy=='Low risk'])*1000*1000 # area in m2
 farea5  = nrow(visualisation.dt[type_qy=='Moderate risk'])*1000*1000 # area in m2
 farea6  = nrow(visualisation.dt[type_qy=='High risk'])*1000*1000 # area in m2
 farea4/(farea4+farea5+farea6)
 farea5/(farea4+farea5+farea6)
 farea6/(farea4+farea5+farea6)
 
 visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("< 0.1", "0.1-0.15", "0.15-0.2","> 0.2"))
 #visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("Low", "Moderate", "Sufficient"))
 #visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("Low risk", "Moderate risk", "High risk"))
 
 #geom_tile
 ggplot() + geom_tile(data = visualisation.dt, aes(x = x, y = y, fill = visualisation.dt$type_qy)) +
   coord_sf(crs = st_crs(QYBorder)) +
   scale_fill_manual(values  = c('< 0.1'="#2c7bb6", '0.1-0.15'="#abdda4", '0.15-0.2'= "#fdae61",'> 0.2'="#d7191c"), na.value = 'transparent')+
   #scale_fill_manual(values  = c('Low'="#2c7bb6",  'Moderate'= "#abdda4",'Sufficient'="#fdae61"), na.value = 'transparent')+
   #scale_fill_manual(values  = c('Low risk'="#abdda4",  'Moderate risk'= "#fdae61",'High risk'="#d7191c"), na.value = 'transparent')+
   #scale_fill_viridis_c() + 
   geom_sf(data = QYBorder) +
   #guides(fill=guide_legend(title = "PSD Status\n\u03b1[1]==0.5; \u03b1[2]==0.5\n[EU base]"))+
   guides(fill=guide_legend(title = "PSD Status\n\u03b1=0.5"))+
   theme(axis.text = element_text(family="A",size = 13,face = "bold",colour ='black'),
         axis.title.x = element_blank(),
         axis.title.y  = element_blank())+
   theme(panel.grid.major = element_line(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA),#
         panel.border = element_rect(fill = NA))+
   theme(legend.position = c(0.2,0.15),
         legend.text = element_text(family="A", size = 11,face = "bold",colour ='black'),
         legend.title = element_text(family="A",size = 11,face = "bold",colour ='black'),
         legend.key =  element_rect(fill = "transparent", colour = NA),
         legend.background = element_blank())
 
########################################alpha =0.8 for regional PSD estimation
# spatial interpolation for psd regional alpha = 0.8
 idw.weighing <- lapply(X = 1:3, FUN = function(x){
   #inverse distance weigh
   distanceweigh <-idw(formula = psd3 ~ 1,#key parameters
                       locations = qy_sf[qy_sf$landuse_code == x,],
                       newdata =QY_rast ,
                       idp = 2)
   #form stars to terra
   distanceweigh <- terra::rast(distanceweigh)
   distanceweigh$var1.var_lyr.1 <- NULL
   names(distanceweigh) <- "psd"
   #resample to landuse
   distanceweigh <- terra::resample(distanceweigh, terra::rast(QY_rast))
 })
 names(idw.weighing) <- c("paddy", "upland", "other")
 
 
 #bind land use map to list
 idw.weighing[["landuse"]] <- terra::rast(QY_rast)
 #stack
 idw.stacked <- terra::rast(idw.weighing)
 #to data.table
 idw.stacked.dt <- terra::as.data.frame(idw.stacked, xy = TRUE) |> as.data.table()
 idw.stacked.dt[landuse == 1, psd := paddy]
 idw.stacked.dt[landuse == 2, psd := upland]
 #we remove show NA for forests and nautral areas since we are only interested in agriculture
 idw.stacked.dt[landuse == 3, psd := NA]
 #to raster
 idw.final <- terra::rast(idw.stacked.dt)
 
#visualisation steps
 visualisation.dt <- terra::as.data.frame(idw.final,xy =TRUE)
 visualisation.dt <- as.data.table(visualisation.dt)
 
 ##alpha=0.8 for regional PSD estimation and 0.8 for TR-agri-QY BASED
 visualisation.dt[, type_qy:= ifelse(psd< as.numeric(dx[3,5]), 'Low',ifelse(psd<as.numeric(dx[3,5]),'Moderate','Sufficient'))]#qy_psd3 (only change CV)
 ##alpha=0.8 for regional PSD estimation and 0.8 for TR-envri-QY BASED
 visualisation.dt[, type_qy:= ifelse(psd< (as.numeric(dx[3,9])*0.5), 'Low risk',ifelse(psd<as.numeric(dx[3,9]),'Moderate risk','High risk'))]#qy_psd3 0.8 (only change CV)
 ##alpha=0.8 for regional PSD estimation and 0.5 for TR-agri-EU BASED
 visualisation.dt[, type_qy:= ifelse(psd< as.numeric(dx[5,5]), 'Low',ifelse(psd<as.numeric(dx[5,7]),'Moderate','Sufficient'))]#eu psd
 ##alpha=0.8 for regional PSD estimation and 0.5 for TR-envri-EU BASED
 visualisation.dt[, type_qy:= ifelse(psd<(as.numeric(dx[5,9])*0.5), 'Low risk',ifelse(psd<as.numeric(dx[5,9]),'Moderate risk','High risk'))]#eu
 

 ## calculate the percentage of area where soil P is lower than thressholds
 farea1  = nrow(visualisation.dt[type_qy=='Low'])*1000*1000 # area in m2
 farea2  = nrow(visualisation.dt[type_qy=='Moderate'])*1000*1000 # area in m2
 farea3  = nrow(visualisation.dt[type_qy=='Sufficient'])*1000*1000 # area in m2
 farea1/(farea1+farea2+farea3)
 farea2/(farea1+farea2+farea3)
 farea3/(farea1+farea2+farea3)
 
 farea4  = nrow(visualisation.dt[type_qy=='Low risk'])*1000*1000 # area in m2
 farea5  = nrow(visualisation.dt[type_qy=='Moderate risk'])*1000*1000 # area in m2
 farea6  = nrow(visualisation.dt[type_qy=='High risk'])*1000*1000 # area in m2
 farea4/(farea4+farea5+farea6)
 farea5/(farea4+farea5+farea6)
 farea6/(farea4+farea5+farea6)
 # 
 # dy = data.table(STP=c('qyR0.5_TR0.8','qyR0.8_TR0.8','euR0.5_TR0.5','euR0.8_TR0.5'),
 #                 ClassA_exL=c( farea1/(farea1+farea2+farea3)),
 #                 ClassA_L= c( farea2/(farea1+farea2+farea3)),
 #                 ClassA_suff= c(farea3/(farea1+farea2+farea3)),
 #                 ClassE_L=c( farea4/(farea4+farea5+farea6)),
 #                 ClassE_M= c( farea5/(farea4+farea5+farea6)),
 #                 ClassE_H= c(farea6/(farea4+farea5+farea6)))
 
 
 #visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("< 0.1", "0.1-0.15", "0.15-0.2","> 0.2"))
 #visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("Low", "Moderate", "Sufficient"))
 visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("Low risk", "Moderate risk", "High risk"))

  ggplot() + geom_tile(data = visualisation.dt, aes(x = x, y = y, fill = visualisation.dt$type_qy)) +
     coord_sf(crs = st_crs(QYBorder)) +
     #scale_fill_manual(values  = c('< 0.1'="#2c7bb6", '0.1-0.15'="#abdda4", '0.15-0.2'= "#fdae61",'> 0.2'="#d7191c"), na.value = 'transparent')+
     #scale_fill_manual(values  = c('Low'="#2c7bb6",  'Moderate'= "#abdda4",'Sufficient'="#fdae61"), na.value = 'transparent')+
     scale_fill_manual(values  = c('Low risk'="#abdda4",  'Moderate risk'= "#fdae61",'High risk'="#d7191c"), na.value = 'transparent')+
     #scale_fill_viridis_c() + 
     geom_sf(data = QYBorder) +
     guides(fill=guide_legend(title = "PSD Status\n\u03b1[1]==0.8; \u03b1[2]==0.8\n[QY base]"))+
     theme(axis.text = element_text(family="A",size = 13,face = "bold",colour ='black'),
           axis.title.x = element_blank(),
           axis.title.y  = element_blank())+
     theme(panel.grid.major = element_line(),
           panel.grid.minor = element_blank(),
           panel.background = element_rect(fill = "transparent", colour = NA),#
           panel.border = element_rect(fill = NA))+
     theme(legend.position = c(0.2,0.15),
           legend.text = element_text(family="A", size = 11,face = "bold",colour ='black'),
           legend.title = element_text(family="A",size = 11,face = "bold",colour ='black'),
           legend.key =  element_rect(fill = "transparent", colour = NA),
           legend.background = element_blank())   
############################# POLSEN
 
 idw.weighing <- lapply(X = 1:3, FUN = function(x){
   #inverse distance weigh
   distanceweigh <-idw(formula = olsenp ~ 1,
                       locations = qy_sf[qy_sf$landuse_code == x,],
                       newdata =QY_rast ,
                       idp = 2)
   
   #form stars to terra
   distanceweigh <- terra::rast(distanceweigh)
   distanceweigh$var1.var_lyr.1 <- NULL
   names(distanceweigh) <- "olsenp"
   
   #resample to landuse
   distanceweigh <- terra::resample(distanceweigh, terra::rast(QY_rast))
 })
 names(idw.weighing) <- c("paddy", "upland", "other")
 
 
 #bind land use map to list
 idw.weighing[["landuse"]] <- terra::rast(QY_rast)
 
 #stack
 idw.stacked <- terra::rast(idw.weighing)
 
 #to data.table
 idw.stacked.dt <- terra::as.data.frame(idw.stacked, xy = TRUE) |> as.data.table()
 idw.stacked.dt[landuse == 1, olsenp := paddy]
 idw.stacked.dt[landuse == 2, olsenp := upland]
 #we remove show NA for forests and nautral areas since we are only interested in agriculture
 idw.stacked.dt[landuse == 3, olsenp := NA]
 #to raster
 idw.final <- terra::rast(idw.stacked.dt)
 
 #visualisation steps
 visualisation.dt <- terra::as.data.frame(idw.final, xy =TRUE)
 visualisation.dt <- as.data.table(visualisation.dt)
 visualisation.dt[, type_qy:= ifelse(olsenp< 15, '< 15',ifelse(olsenp< 35,'15-35',ifelse(olsenp<70,'35-70','> 70')))]
 
 visualisation.dt[, type_qy:= ifelse(olsenp< as.numeric(dx[1,5]), 'Low',ifelse(olsenp<as.numeric(dx[1,7]),'Moderate','Sufficient'))]#qy
 visualisation.dt[, type_qy:= ifelse(olsenp< as.numeric(dx[4,5]), 'Low',ifelse(olsenp<as.numeric(dx[4,7]),'Moderate','Sufficient'))]#eu
 visualisation.dt[, type_qy:= ifelse(olsenp< 250, 'Low risk',ifelse(olsenp<300,'Moderate risk','High risk'))]#qy
 visualisation.dt[, type_qy:= ifelse(olsenp< (as.numeric(dx[4,9])*0.5), 'Low risk',ifelse(olsenp<as.numeric(dx[4,9]),'Moderate risk','High risk'))]#eu

 ## calculate the percentage of area where soil P is lower than thressholds
 farea1  = nrow(visualisation.dt[type_qy=='Low'])*1000*1000 # area in m2
 farea2  = nrow(visualisation.dt[type_qy=='Moderate'])*1000*1000 # area in m2
 farea3  = nrow(visualisation.dt[type_qy=='Sufficient'])*1000*1000 # area in m2
 farea1/(farea1+farea2+farea3)
 farea2/(farea1+farea2+farea3)
 farea3/(farea1+farea2+farea3)
 
 farea4  = nrow(visualisation.dt[type_qy=='Low risk'])*1000*1000 # area in m2
 farea5  = nrow(visualisation.dt[type_qy=='Moderate risk'])*1000*1000 # area in m2
 farea6  = nrow(visualisation.dt[type_qy=='High risk'])*1000*1000 # area in m2
 farea4/(farea4+farea5+farea6)
 farea5/(farea4+farea5+farea6)
 farea6/(farea4+farea5+farea6)
 
 
 visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("< 15", "15-35", "35-70","> 70"))
 #visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("Low", "Moderate", "Sufficient"))
 #visualisation.dt$type_qy <- factor(visualisation.dt$type_qy, levels = c("Low risk", "Moderate risk", "High risk"))

 #geom_tile
 ggplot() + geom_tile(data = visualisation.dt, aes(x = x, y = y, fill = visualisation.dt$type_qy)) +
   coord_sf(crs = st_crs(4326)) +
   scale_fill_manual(values  = c('< 15'="#2c7bb6", '15-35'="#abdda4", '35-70'= "#fdae61",'> 70'="#d7191c"), na.value = 'transparent')+
   #scale_fill_manual(values  = c('Low'="#2c7bb6",  'Moderate'= "#abdda4",'Sufficient'="#fdae61"), na.value = 'transparent')+
   #scale_fill_manual(values  = c('Low risk'="#abdda4",  'Moderate risk'= "#fdae61",'High risk'="#d7191c"), na.value = 'transparent')+
   geom_sf(data = QYBorder) +
   guides(fill=guide_legend(title = "Olsen P Status\n [mg/kg]"))+
   theme(axis.text = element_text(family="A",size = 13,face = "bold",colour ='black'),
         axis.title.x = element_blank(),
         axis.title.y  = element_blank())+
   theme(panel.grid.major = element_line(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA),#
         panel.border = element_rect(fill = NA))+
   theme(legend.position = c(0.2,0.15),
         legend.text = element_text(family="A", size = 11,face = "bold",colour ='black'),
         legend.title = element_text(family="A",size = 11,face = "bold",colour ='black'),
         legend.key =  element_rect(fill = "transparent", colour = NA),
         legend.background = element_blank())  
 
###############################alpha value calculation
 library(forcats)
 alpha<- fread("data/experimentalData/Alpha_data.csv")
 cols_old <- colnames(alpha)
 cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
 setnames(alpha,cols_old,cols_new)
 alpha[,mean_alpha_em:= mean(alpha_em)]
 alpha$climate2 <- fct_relevel(alpha$climate2, "Tropical", "Subtropical", "Temperate")
 
 ggplot(aes(x = alpha_em, fill= climate2), data=alpha) + 
   geom_histogram(aes(y =..density..),
                  breaks = seq(0.1, 1.8, by = 0.1),colour= "black") + 
   stat_function(fun = dnorm, args = list(mean = mean(alpha$alpha_em), sd=sd(alpha$alpha_em)))+
   geom_vline(data = alpha, aes(xintercept = mean_alpha_em, color = 'balck'), linetype = "dashed", color= "black")+
   geom_rug(aes(fill=climate2)) + facet_wrap(~climate2) +
   theme_bw() +
   scale_fill_manual(values = c("#a50026", "#f46d43", "#d9ef8b"))+
   #labs(title = "Alpha Values Distribution by Climate", x = "\u03b1", y = "Density")+
   labs( x = "\u03b1", y = "Density")+
   guides(fill=guide_legend(title = "Climate"))+
   theme(axis.text = element_text(family="A",size = 9,colour ="black" ),
         axis.title.x = element_text(family="A",size = 10,colour ="black",face = "bold"),
         axis.title.y  = element_text(family="A",size = 10,colour ="black",face = "bold"))+
   theme(panel.grid.major = element_line(colour = NA),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA),
         panel.border = element_rect(fill = NA))+
   theme(legend.position = c(0.9,0.8),
         legend.text = element_text(family="A", size = 7,colour ="black"),
         legend.title = element_text(family="A",size = 8,face = "bold",colour ="black"),
         legend.key =  element_rect(fill = "transparent", colour = NA),
         legend.background = element_blank())
 
 
 