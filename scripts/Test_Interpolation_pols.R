## test for POlsen interpolation for the whole China but too slow if make the grid like 100*100 

library(sf);library(raster);library(terra);library(gstat);library(stars)
 ##load data for CN shape
  
  CN <- st_read("data/spatialInfor/simplied_china_country.shp")
  CN_sf <-st_transform(CN, "epsg:32649")#
  CNBorder <- st_as_stars(st_cast(CN_sf, "MULTILINESTRING"))# 

 ## load fertilization dataset 
  fert <- fread("data/experimentalData/fertilizationData_china.csv")
  
  #select related variables, x, y, olsenP and SOM
  
  fert_cn <- fert[,c("latitute","longtitude","som","pols")]
  
  #project and set coordnate systems
  fert_cn <- st_as_sf(fert_cn, coords= c("latitute","longtitude"), crs = 4214)
  
  fert_cn_sf<-st_transform(fert_cn , "epsg:32649")
  fert_cn_sf <- na.omit(fert_cn_sf)
  
 ## make the grid for CN border
  
  cngrid <- st_make_grid(CNBorder,n=c(100,100))

 ## idw interpolation for whole CN
  
  pols <- idw(formula = pols ~ 1,#key parameters
                     locations = fert_cn_sf,
                      newdata =cngrid,
                      idp = 2)
  #save interplotion results since too slow
  save(pols,file = "dev/idw1.Rdata")
  
 ##crop CN region and make sure same the CNBroder
  
  pols <- st_intersection(pols, st_union(CNBorder))
 
 ##visualization

  #convert raster to data.frame
  pols_rast <- rast(pols)


pols_raster <- raster(pols_rast)
pols_spdf <- as(pols_raster, "SpatialPixelsDataFrame")
pols_df <- as.data.frame(pols_spdf)


ggplot()+
  geom_raster(data = pols_df , aes(x = x, y = y,fill = var1.pred_lyr.1))+
  geom_sf(data = QYBorder, color="black", size=0.5)+
  labs(y = 'lon', x = 'lat')+theme_bw()+
  scale_fill_gradient2(low ="#fee8c8",high="#e34a33",mid = "#fdbb84", space ="Lab", midpoint = 14 )
