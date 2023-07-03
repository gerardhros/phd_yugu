
 ##QY_data is origional dataset, and QY_data3 is dataset after finish outlier test
 QY_data <- fread("data/experimentalData/qiyang-inventory_coodinate.csv")
 cols_old <- colnames(QY_data)
 cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))# doesn't work for read.csv
 setnames(QY_data,cols_old,cols_new)
 
 ##calculation (firstly all the parameters are the simplest one)
 ## assume bulk density is same in the whole county 1.1 kg/m3
 ## assume thickness is 0.2m
 ## assume wheat-maize rotation is the only crop system for Qiyang county
 ## assume the max maize and crop yield from Qiyang long term term exp are the target yield for traget p uptake
 ## 1272 kg/ha for wheat and 5096 kg/ha for maize
 ## grain and straw P concentration for wheat: 3.5 and 0.8 g/kg, 1.3 for the ratio straw: grain
 ## grain and straw P concentration for maize: 2.18 and 1.52 g/kg, 1.1 for the ratio straw: grain
 ## assume the target Polsen is 28 mg/kg for both wheat and maize
 ## assume the traget PSD is 0.2
 
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
   totalp <- 15.83*x^0.65/(1.4*10)
 }
 
 oslenp <- function(x) {
   # Define the function for Olsen P as a function of x
   # Replace with your own function or data points
   # For example: olsenp <- 0.5 * x^2 + 2 * x + 4
   olsenp <- 493*x/(1.4*10)### assume 1.4*10 is the parameter from kg/ha to mg/kg
 }
 
 from <- 0
 
 to <-3200
 
 # Calculate the integral average ratio for a specific starting point
 average_ratio <- calculate_integral_average(totalp, olsenp, from, to)
 
 
 ##---calculate for required amount based on P olsen
 
 QY_data[, psd:= (oxalatep/31)/(0.5*((oxalateal/27)+(oxalatefe/56)))]
 
 QY_data[,priOls_mid:= (1272*3.5+ 1272*1.3*0.8+5096*2.18+75096*1.52*1.1)/1000000+
           (1.1*0.2*(28-olsenp)*average_ratio*0.01)/35]
 

 ## same with PSD but be carful with 0.5
 
 