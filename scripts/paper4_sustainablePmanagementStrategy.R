
 ##QY_data is origional dataset, and QY_data3 is dataset after finish outlier test
 QY_data <- fread("data/experimentalData/qiyang-inventory_coodinate.csv")
 cols_old <- colnames(QY_data)
 cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))# doesn't work for read.csv
 setnames(QY_data,cols_old,cols_new)
 
 ##calculation (firstly all the parameters are the simplest one)
 
 ##---------calculate integral average Ptotal/P-Av ratio  going from current to the critical P status------
 # Define the integration function using numerical integration, x is accumulated P surplus
 calculate_integral_average <- function(totalp, olsenp, from, to) {
   # Define the integrand function
   integrand <- function(x) {
     ratio_func(totalp(x), olsenp(x))
   }
   
   # Perform numerical integration
   integral <- integrate(integrand, from = from(), to = to())
   
   # Return the average ratio
   return(integral$value / (to() - from()))
 }
 # Usage example
 totalp <- function(x) {
   
   # Define the function for Total Phosphorus (P) as a function of x
   # Replace with your own function or data points
   # For example: totalp <- 2 * x^2 + 3 * x + 1
   totalp <- 15.83*x^0.65/(1.4*10)### assume 1.4*10 is the parameter from kg/ha to mg/kg
   #the relationship between total P and accumulated P surplus was derived in paper1
 }
 oslenp <- function(x) {
   # Define the function for Olsen P as a function of x
   # Replace with your own function or data points
   # For example: olsenp <- 0.5 * x^2 + 2 * x + 4
   olsenp <- 493*x/(1.4*10)### assume 1.4*10 is the parameter from kg/ha to mg/kg
   # the relationship between olsen P and accumulated P surplus was derived in paper1
 }
 from <- 0
 to <-3200 # 0 and 3200 are also from paper 1(when soil reactive P pools were saturated)
 # Calculate the integral average ratio for a specific starting point
 average_ratio <- calculate_integral_average(totalp, olsenp, from, to)
 
 ##----- 
 
 myield_tar <- 5096 #  maximum maize yield kg/ha from paper1
 wyield_tar <- 1272 #  maximum wheat yield kg/ha from paper1
 mgrainP_con <- 2.18 # maize grain P concentraion g/kg, ref(Qichao phd thesis)
 mstrawP_con <- 1.52 # maize straw P concentraion g/kg, ref(Qichao phd thesis)
 mstraw_grainRatio <- 1.1 #ratio for maize straw over by grain 
 wgrainP_con <- 3.5 # wheat grain P concentraion g/kg, ref(Qichao phd thesis)
 wstrawP_con <- 0.8 # wheat straw P concentraion g/kg, ref(Qichao phd thesis)
 wstraw_grainRatio <- 1.3 # #ratio for wheat straw over by grain 
 soildepth <- 0.2 # assume thickness is 0.2m
 bulkDensity <- 1250   # assume a bulk density near 1250 kg m-3
 pol_tar <- 28   # assume the target Polsen is 28 mg/kg for both wheat and maize
 psd_tar <- 0.2  # assume the traget PSD is 0.2
 T_tar <- 35  # assume target time is 35 years
 
 ##-----calculate for required amount based on P olsen
 # Pin,req,mt = Pup,target  + (ρ x Tsoil x (POlsentarget – POlsencurrent) x Ptotal/POlsen x 0.01)/Ttarget
 QY_data[,priOls_mid:= (wyield_tar*wgrainP_con+ wyield_tar*wstrawP_con*wstraw_grainRatio+myield_tar*mgrainP_con+myield_tar*mstrawP_con*mstraw_grainRatio)/1000+
           (bulkDensity*soildepth*(pol_tar-olsenp)*(totalp*1000/olsenp/31)*0.01)/T_tar]
 
 # same with POLSEN but be careful with 0.5*[Feox+Alox]....
 QY_data[, psd:= (oxalatep/31)/(0.5*((oxalateal/27)+(oxalatefe/56)))]
 QY_data[,priPSD_mid:= (wyield_tar*wgrainP_con+ wyield_tar*wstrawP_con*wstraw_grainRatio+myield_tar*mgrainP_con+myield_tar*mstrawP_con*mstraw_grainRatio)/1000+
           (soildepth*bulkDensity*(psd_tar-psd)*((oxalateal/27)+(oxalatefe/56))*((totalp*1000/oxalatep)/31)*31*0.01)/T_tar]
 
 ## plot spatial P required amount
 library(sf);library(raster);library(terra);library(gstat);library(stars)
 
 ##load data and set coordinate system for the county 
 QY <- st_read("data/spatialInfor/Qiyang.shp") 
 QY_map <- QY[,c(2,5)]
 QY_sf <-st_transform(QY_map, "epsg:32649")#
 QYBorder <- st_cast(QY_sf, "MULTILINESTRING")# 

 #raster data for landuse map,maybe useful later
 # load land use raster data 
 
 # QYtif <- raster('data/spatialInfor/ld2015.tif') 
 # 
 # # specify project coordinate system 
 # 
 # QYtif <-projectRaster(QYtif, crs = "epsg:32649")
 
 #selected data for IDW interploation 
 
 qy <- QY_data[,c("lat","lon","priPSD_mid","priOls_mid")]
 
 # specify coordinates
 
 coordinates(qy) <- ~lon+lat
 
 # specify coordinates reference system
 
 proj4string(qy) <- CRS("EPSG:32649")

 # idw interpolation for required P amount
 
 pols <-gstat::idw(formula = priPSD_mid ~ 1,#key parameters
                     locations = qy,
                     newdata = QY_sf,
                     idp = 1)
 
 # visualize P required amount for the county
 
 ggplot() + geom_tile(data = pols, aes(x = x, y = y, fill= var1.pred)) +
   coord_sf(crs = st_crs(QY_sf)) 
