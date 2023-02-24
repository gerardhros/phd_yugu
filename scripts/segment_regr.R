require(data.table)
require(segmented)
setwd('D:/ESA/02 phd projects/03 yu gu')

#load data from Nawara
ry_eu1<- fread("01 data/20220906Nawara_data.csv")
cols_old <- colnames(ry_eu1)
cols_new <- tolower(unlist(tstrsplit(cols_old,'\\[',keep=1)))
setnames(ry_eu1,cols_old,cols_new)
ry_eu1[,psd:= `p-ox (mmol/kg)`/(0.5*(`al-ox (mmol/kg)` +`fe-ox (mmol/kg)`))]
ry_eu <- ry_eu1[gewassen=='wheat'|gewassen=='maize']#avoid change all the names in the following steps
ry_eu[, max_yield_year:= max(yield), by =c("year","site") ]
ry_eu[, relative_yield_eu:= yield/max_yield_year*100]
#load data for Qiyang(selected related data not origional)
ry_qy <- fread("01 data/segmented_regre.csv")
l3_qy3 <- lm(cacl2p~psd3,data=ry_qy)
os3_qy3 <- segmented(l3_qy3,~psd3, psi =0.28)#meothod 1
os33_qy3<- glm(cacl2p ~ psd3 + I(pmax(0,psd3 - 0.28)), data=ry_qy)#method 2

#plot
require(ggplot2)
ggplot()+
  geom_point(aes(x=psd,y=p_cacl2),size=2,colour='#1f78b4',data=ry_eu)+
  geom_point(aes(x=psd,y=cacl2p),size=2,colour='#d7191c',data=ry_qy )+#0.5
  # geom_point(aes(x=psd2,y=cacl2p),color= '#d7191c', size=2,
  #            data = ry_qy)+#
  geom_point(aes(x=psd3,y=cacl2p),color= '#fdae61', size=2,
             data = ry_qy)+
  geom_smooth(method="lm",aes(x=psd,y=p_cacl2),
              formula=y ~ bs(x,degree=1,knots=c(0.33)),
              data=ry_eu,se=F,color='#1f78b4')+
  geom_smooth(method="lm",aes(x=psd,y=cacl2p),
              formula=y ~ bs(x,degree=1,knots=c(0.32)),
              data=ry_qy,se=F,color='#d7191c')+
  geom_smooth(method="lm",aes(x=psd3,y=cacl2p),
              formula=y ~ bs(x,degree=1,knots=c(0.28)),
              data=ry_qy,se=F,color='#fdae61')+
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
