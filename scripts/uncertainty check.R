##PSD psd(alpha=0.5)  psd3(alpha=0.8)
library(segmented)
library(data.table)

ry_eu <- fread('D:/ESA/02 phd projects/03 yu gu/01 data/ry_eu.csv')
#Pcacl2 vs PSD
l3_eu <- lm(p_cacl2~ psd,data=ry_eu)
os3_eu <- segmented(l3_eu,~ psd,psi = 0.1)#meothod 1
os33_eu <- glm(p_cacl2 ~ psd + I(pmax(0,psd - 0.33)), data=ry_eu)#method 2

confint(os3_eu) #delta CI for the 1st variable 
confint(os3_eu, "psd", method="score") #also method="g"
cint2 <- confint(os3_eu, "psd", method="g")

l3_eu3 <- lm(p_cacl2~ psd3,data=ry_eu)
os3_eu3 <- segmented(l3_eu3,~ psd3,psi = 0.1)#meothod 1
os33_eu3 <- glm(p_cacl2 ~ psd3 + I(pmax(0,psd3 - 0.21)), data=ry_eu)#method 2
##relative crop yield vs psd
ry_psd_eu <-  nls(relative_yield_eu~SSasymp(psd,Asym,R0,lrc), data = ry_eu)
ry_psd3_eu <-  nls(relative_yield_eu~SSasymp(psd3,Asym,R0,lrc), data = ry_eu)

pred <- predict(os33_eu,newdata = ry_eu)
plot(ry_eu$p_cacl2~ry_eu$psd)
points(pred~ry_eu$psd,pch=16,col='blue')
abline(v = cint[2],lty=2)
abline(v = cint[3],lty=2)
