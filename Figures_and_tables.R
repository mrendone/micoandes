
###########################################
#                                         #
#           Figures                       #
#                                         #
###########################################

library(dplyr)
library(readxl)


Data10kmFiltrado13_12_20 <- read_excel("./insumos/Data10kmFiltrado13_12_20.xlsx")
Data10kmFiltrado13_12_20<- Data10kmFiltrado13_12_20 %>% rename(id="IDceldas10km.IDcelda10km")
mayores_a_100m<-Data10kmFiltrado13_12_20 %>% filter(Elev10kmcell>100)

library(ggpubr)
library(cowplot)
library(ggplot2)
############
# Figure 1 #
############


library(raster)
library(viridis)
library(ggspatial)


ColomabiaNew<-raster::bind(andes, Colombia) 

andes<-st_read("./insumos/Andes.shp")
andes<-as(andes, "Spatial")

Colombia<-st_read("./insumos/Colombia.shp")
Colombia<-as(Colombia, "Spatial")

spat_mico<-st_read("./insumos/Spat_mico.shp")
spat_mico<-as(spat_mico,"Spatial")


Colombia_sf<-st_as_sf(ColomabiaNew)
cordillera_sf<-st_as_sf(andes)


spat_mico<-st_read("./insumos/Spat_mico.shp")
spat_mico<-as(spat_mico,"Spatial")


bb10 <- bbox(andes) #This command retrieves the bounding box of the spatial object andes and stores it in the variable bb10
cs10 <- c(1/111.325,1/111.325)*10  # This line defines the cell size for a grid.
cc10 <- bb10[, 1] + (cs10/2)  # calculate the cell center offset.
cd10 <- ceiling(diff(t(bb10))/cs10)  # calculate the number of cells per direction
grd10 <- GridTopology(cellcentre.offset=cc10, cellsize=cs10, cells.dim=cd10)#creates a grid topology with the cell center offset, cell size, and number of cells per direction.
grd10
sp_grd10 <- SpatialGridDataFrame(grd10,
                                 data=data.frame(id=1:prod(cd10)),
                                 proj4string=CRS(proj4string(andes))) #This line creates a SpatialGridDataFrame from the grid topology and some associated data.



spPoly<-as(sp_grd10, "SpatialPolygonsDataFrame")


spPoly@data<- full_join(spPoly@data, Data10kmFiltrado13_12_20, by="id")

andes_sf<-st_as_sf(spPoly)
andes_sf<-na.omit(andes_sf)
Colombia_sf<-st_as_sf(ColomabiaNew)
cordillera_sf<-st_as_sf(andes)

andes_sf<-andes_sf %>% rename(Elevation=Elev10kmcell, AM_ratio=AM_ratiocell, EM_ratio=EM_ratiocell,
                              ErM_ratio=ErM_ratiocell, OM_ratio=OM_ratiocell, WanNm_ratio=WanNm_ratiocell,
                              Nfix_ratio=Nfix_ratiocell, Shannon_Index=shannon10kmRegistros)



Elev_C<-ggplot(data=Colombia_sf)+
  geom_sf()+
  geom_sf(data = cordillera_sf)+
  geom_sf(data = andes_sf,aes(fill = Elevation, color = Elevation), lwd=0) +
  scale_fill_viridis_c(option = "plasma")+
  annotation_scale(location = "bl", width_hint = 0.4)+
  theme_bw()





##########
#Figure 2#
##########

#command is reading a CSV file named “tabShaEnvirRaref_50.csv” that is a rarefiated data to 50 samples per cell

tabShaEnvirRaref_50<-read.csv("./insumos/tabShaEnvirRaref_50.csv")

#Here it is extracting the row names from the dataframe tabShaEnvirRaref_50 and storing them in a variable called index

index<-row.names(tabShaEnvirRaref_50)
tabShaEnvirRaref_50$id<-as.numeric(index)

#Here it is selecting columns 1 and 53 from the dataframe Data10kmFiltrado13_12_20 and storing them in a variable called nitro

nitro<-Data10kmFiltrado13_12_20[,c(1,53)]

#This command is performing an “inner join” of the dataframes tabShaEnvirRaref_50 and nitro based on the id column. The result is stored back in tabShaEnvirRaref_50.

tabShaEnvirRaref_50<-inner_join(tabShaEnvirRaref_50, nitro, by="id")

#Shannon graph and shannon rarefied graphs  vs elevation and Carbon stock

shannon10kmRegistrosGraphElev<-ggplot(data=mayores_a_100m, aes(Elev10kmcell, shannon10kmRegistros))+
  geom_point()+
  geom_smooth(method = "lm", color = "#1f77b4")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text= element_text(size=12),axis.title.x= element_text(color="black", size=15), axis.title.y= element_text(color="black", size=15))+
  labs(x="Elevation", y= "Diversity of types of mycorrhizal associations", size=12)+
  stat_cor(method = "pearson",cor.coef.name=c("r"))+
  ggtitle("Original Data")

RAFshannonRegistrosGraphElev<-ggplot(tabShaEnvirRaref_50, aes(x = Elev10kmcell, y = mean_sha)) +
  geom_point()+
  geom_smooth(method = "lm", color = "#1f77b4")+
  geom_segment(aes(xend = Elev10kmcell, yend = s_p5), linetype = "dotdash", color = "grey") +
  geom_segment(aes(xend = Elev10kmcell, yend = s_p95), linetype = "dotdash", color = "grey") +
  
  geom_point(color = "black", size = 0.75)+ stat_cor(method = "pearson",cor.coef.name=c("r")) + labs(x="Elevation", y= "Diversity of types of mycorrhizal associations")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text= element_text(size=12), axis.title.x= element_text(color="black", size=15), axis.title.y= element_text(color="black", size=15))+
  ggtitle("Rarified Data")

shannon10kmRegistrosGraphOCSTHA<-ggplot(data=mayores_a_100m, aes(OCSTHA, shannon10kmRegistros))+
  geom_point()+
  geom_smooth(method = "lm", color = "#1f77b4")+
  theme_classic()+
  theme(axis.text= element_text(size=12), axis.title.x= element_text(color="black", size=15), axis.title.y= element_text(color="black", size=15))+
  labs(x="Soil organic carbon stock (ton/ha)", y= "Diversity of types of mycorrhizal associations", size=12) + stat_cor(method = "pearson",cor.coef.name=c("r"))

RAFshannonRegistrosGraphOCSTHA<-ggplot(tabShaEnvirRaref_50, aes(x = OCSTHA, y = mean_sha)) +
  geom_point()+
  geom_smooth(method = "lm", color = "#1f77b4")+
  geom_segment(aes(xend = OCSTHA, yend = s_p5), linetype = "dotdash", color = "grey") +
  geom_segment(aes(xend = OCSTHA, yend = s_p95), linetype = "dotdash", color = "grey") +
  
  geom_point(color = "black", size = 0.75)+ stat_cor(method = "pearson",cor.coef.name=c("r")) + labs(x="Soil organic carbon stock (ton/ha)", y= "Diversity of types of mycorrhizal associations")+
  theme_classic()+
  theme(axis.text= element_text(size=12), axis.title.x= element_text(color="black", size=15), axis.title.y= element_text(color="black", size=15))

# Fit lm of shannon10kmRegistros ~ Nitrogeno

modelo <- lm(shannon10kmRegistros ~ Nitrogeno, data = mayores_a_100m)

# Predict values for x
predicciones <- predict(modelo, newdata = mayores_a_100m)

# Get confidence intervals
conf_int <- predict(modelo, newdata = mayores_a_100m, interval = "confidence")

# Create a new data set with only the points where the predictions are >= 0
datos_filtrados <- subset(mayores_a_100m, predicciones >= 0)

# Filter predictions and confidence intervals for filtered data
predicciones_filtradas <- predicciones[predicciones >= 0]
conf_int_filtrado <- conf_int[predicciones >= 0, ]

# Add the filtered predictions and confidence intervals to the filtered data set
datos_filtrados <- cbind(datos_filtrados, predicciones = predicciones_filtradas, conf_int_filtrado)

# Plot the data, linear regression line and confidence intervals
shannon10kmRegistrosGraphNitro <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Nitrogeno, y = shannon10kmRegistros)) +
  geom_line(data = datos_filtrados, aes(x = Nitrogeno, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +  
  geom_ribbon(data = datos_filtrados, aes(x = Nitrogeno, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Total Nitrogeno in soil (cg/kg)", y = "Diversity of types of mycorrhizal associations")

# Add Pearson correlation coefficient and p value
shannon10kmRegistrosGraphNitro <- shannon10kmRegistrosGraphNitro +
  geom_text(data = NULL, aes(label = paste("r =", round(cor(mayores_a_100m$Nitrogeno, mayores_a_100m$shannon10kmRegistros), 2), ", p < 2.2e- 16 ")), x = min(mayores_a_100m$Nitrogeno), y = max(mayores_a_100m$shannon10kmRegistros), hjust = 0, vjust = 1, size = 4, color = "black")

RAFshannonRegistrosGraphNitrogeno<-ggplot(tabShaEnvirRaref_50, aes(x = Nitrogeno, y = mean_sha)) +
  geom_point()+
  geom_smooth(method = "lm", color = "#1f77b4")+
  geom_segment(aes(xend = Nitrogeno, yend = s_p5), linetype = "dotdash", color = "grey") +
  geom_segment(aes(xend = Nitrogeno, yend = s_p95), linetype = "dotdash", color = "grey") +
  
  geom_point(color = "black", size = 0.75)+ stat_cor(method = "pearson",cor.coef.name=c("r")) + labs(x="Total Nitrogeno in soil (cg/kg)", y= "Diversity of types of mycorrhizal associations")+
  theme_classic()+
  theme(axis.text= element_text(size=12), axis.title.x= element_text(color="black", size=15), axis.title.y= element_text(color="black", size=15))

#preview of figure 2
plot_grid(shannon10kmRegistrosGraphElev,RAFshannonRegistrosGraphElev,shannon10kmRegistrosGraphOCSTHA,RAFshannonRegistrosGraphOCSTHA,shannon10kmRegistrosGraphNitro,RAFshannonRegistrosGraphNitrogeno, labels=c("a","b","c","d","e","f"), ncol=2,label_size = 15)


#Save figure 2 in a variable 
fig2<-plot_grid(shannon10kmRegistrosGraphElev,RAFshannonRegistrosGraphElev,shannon10kmRegistrosGraphOCSTHA,RAFshannonRegistrosGraphOCSTHA,shannon10kmRegistrosGraphNitro,RAFshannonRegistrosGraphNitrogeno, labels=c("a","b","c","d","e","f"), ncol=2,label_size = 20)


##########
#Figure 3#
##########

Data10kmFiltrado_4000 <-Data10kmFiltrado13_12_20%>%filter(Elev10kmcell<=4000)
Data10kmFiltrado_200 <-Data10kmFiltrado13_12_20%>%filter(Elev10kmcell>=200)

#lm; mycorrhizal types ratio ~ Elev10kmcell

AM_ratiocell10kmGLM<-glm(AM_ratiocell~ Elev10kmcell,
                         data = mayores_a_100m, family=gaussian(link = "identity")) 

EM_ratiocell10kmGLM<-glm(EM_ratiocell~ Elev10kmcell,
                         data = mayores_a_100m, family=gaussian(link = "identity")) 

ErM_ratiocell10kmGLM<-glm(ErM_ratiocell~ Elev10kmcell,
                          data = Data10kmFiltrado_200, family=gaussian(link = "identity")) 

OM_ratiocell10kmGLM<-glm(OM_ratiocell~ Elev10kmcell,
                         data = mayores_a_100m, family=gaussian(link = "identity")) 

WanM_ratiocell10kmGLM<-glm(WanNm_ratiocell~ Elev10kmcell,
                           data = mayores_a_100m, family=gaussian(link = "identity")) 

Nfix_ratiocell10kmGLM<-glm(Nfix_ratiocell~ Elev10kmcell,
                           data = Data10kmFiltrado_4000, family=gaussian(link = "identity")) 

#Graph of mycorrhizal ratios 

Nfix_ratiocell10kmGraph<-ggplot(data=Data10kmFiltrado_4000, aes(Elev10kmcell, Nfix_ratiocell))+
  geom_point()+
  xlim(min(mayores_a_100m$Elev10kmcell), max(4500))+
  geom_smooth(method = "lm", color = "#1f77b4")+
  theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title.x= element_text(color="black", size=13), axis.title.y= element_text(color="black", size=13))+
  labs(x="Elevation", y= "Proportion of Nfixing") + 
  stat_cor(method = "pearson",cor.coef.name=c("r"))



EM_ratiocell10kmGraph<-ggplot(data=mayores_a_100m, aes(Elev10kmcell, EM_ratiocell))+
  geom_point()+
  geom_smooth(method = "glm", color = "#1f77b4")+
  theme_classic()+
  theme(axis.text = element_text(size = 12),axis.title.x= element_text(color="black", size=13), axis.title.y= element_text(color="black", size=13))+
  labs(x="Elevation", y= "Proportion of Em") + stat_cor(method = "pearson",cor.coef.name=c("r"))


AM_ratiocell10kmGraph<-ggplot(data=mayores_a_100m, aes(Elev10kmcell, AM_ratiocell))+
  geom_point()+
  geom_smooth(method = "lm", color = "#1f77b4")+
  theme_classic()+
  theme(axis.text = element_text(size = 12),axis.title.x= element_text(color="black", size=13), axis.title.y= element_text(color="black", size=13))+
  labs(x="Elevation", y= "Proportion of Am") + stat_cor(method = "pearson",cor.coef.name=c("r"),label.x =1000 , label.y = 0.25)

OM_ratiocell10kmGraph<-ggplot(data=mayores_a_100m, aes(Elev10kmcell, OM_ratiocell))+
  geom_point()+
  geom_smooth(method = "glm", color = "#1f77b4")+
  theme_classic()+
  theme(axis.text = element_text(size = 12),axis.title.x= element_text(color="black", size=13), axis.title.y= element_text(color="black", size=13))+
  labs(x="Elevation", y= "Proportion of OM") + stat_cor(method = "pearson",cor.coef.name=c("r"))


#AQUÍ TUVE QUE HACER EL CORTE DE DATOS PARA QUE NO QUEDE BAJO CERO EN ErM#

# Filter the data for those with prediction values >= 0

mayores_a_100m_filtered <- subset(mayores_a_100m, predict(lm(ErM_ratiocell ~ Elev10kmcell, mayores_a_100m)) >= 0)

# Plot the filtered data and add the linear regression line
# fit lm

modelo <- lm(ErM_ratiocell ~ Elev10kmcell, data = mayores_a_100m)

# Predict values for x

predicciones <- predict(modelo, newdata = mayores_a_100m)

# Get confidence intervals
conf_int <- predict(modelo, newdata = mayores_a_100m, interval = "confidence")

# Create a new data set with only the points where the predictions are >= 0
datos_filtrados <- subset(mayores_a_100m, predicciones >= 0)

# Filter predictions and confidence intervals for filtered data

predicciones_filtradas <- predicciones[predicciones >= 0]
conf_int_filtrado <- conf_int[predicciones >= 0, ]

# Add the filtered predictions and confidence intervals to the filtered data set
datos_filtrados <- cbind(datos_filtrados, predicciones = predicciones_filtradas, conf_int_filtrado)

# Plot the data, linear regression line and confidence intervals

ErM_ratiocell10kmGraph <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Elev10kmcell, y = ErM_ratiocell)) + 
  geom_line(data = datos_filtrados, aes(x = Elev10kmcell, y = predicciones),  color = "#1f77b4", linetype = "solid", linewidth = 1.3)+
  geom_ribbon(data = datos_filtrados, aes(x = Elev10kmcell, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  geom_text(data = NULL, aes(label = paste("r =", round(cor(mayores_a_100m$Elev10kmcell, mayores_a_100m$shannon10kmRegistros), 2), ", p < 2.2e-16")), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 4, color = "black") + labs(x = "Elevation", y = "Proportion of ErM")

print(ErM_ratiocell10kmGraph)





WanNm_ratiocell10kmGraph<-ggplot(data=mayores_a_100m, aes(Elev10kmcell, WanNm_ratiocell))+
  geom_point()+
  geom_smooth(method = "lm", color = "#1f77b4")+
  theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title.x= element_text(color="black", size=13), axis.title.y= element_text(color="black", size=13))+
  labs(x="Elevation", y= "Proportion of WanNm") + stat_cor(method = "pearson",cor.coef.name=c("r"))

#Create figure 3

plot_grid(AM_ratiocell10kmGraph, EM_ratiocell10kmGraph, ErM_ratiocell10kmGraph,  OM_ratiocell10kmGraph, 
          WanNm_ratiocell10kmGraph, Nfix_ratiocell10kmGraph, labels=c("a","b","c","d","e","f"))

#Save figure 3 in a variable 


grap3<-plot_grid(AM_ratiocell10kmGraph, EM_ratiocell10kmGraph, ErM_ratiocell10kmGraph,  OM_ratiocell10kmGraph, 
                 WanNm_ratiocell10kmGraph, Nfix_ratiocell10kmGraph, labels=c("a","b","c","d","e","f"))


#######################################
#  Figure 5 and supplementary tables  #
#######################################

library(MASS)
library(car)

# Identify Final models with stepAIC

AM_CarbonoAIC<-stepAIC(glm(AM_ratiocell~ BIO1 + BIO12 + PHIHOX + OCSTHA, 
                           data =mayores_a_100m, gaussian(link = "identity")))


ModeloFinal_AM_Carbono<- glm(formula = AM_ratiocell ~ BIO1 + BIO12 + OCSTHA, family = gaussian(link = "identity"), 
                             data = mayores_a_100m)

AM_NitrogenoAIC<-stepAIC(glm(AM_ratiocell~ BIO1 + BIO12 + PHIHOX + Nitrogeno, 
                             data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_AM_Nitrogeno<- glm(formula = AM_ratiocell ~ BIO1 + BIO12 + Nitrogeno, family = gaussian(link = "identity"), 
                               data = mayores_a_100m)

EM_CarbonoAIC<-stepAIC(glm(EM_ratiocell~ BIO1 + BIO12  + PHIHOX + OCSTHA, 
                           data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_EM_Carbono<-glm(formula = EM_ratiocell ~ BIO1 + BIO12, family = gaussian(link = "identity"), 
                            data = mayores_a_100m)

EM_Nitrogeno<-stepAIC(glm(EM_ratiocell~ BIO1 + BIO12  + PHIHOX + Nitrogeno, 
                          data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_EM_Nitrogeno<-glm(formula = EM_ratiocell ~ BIO1 + BIO12, family = gaussian(link = "identity"), 
                              data = mayores_a_100m)

ErM_Carbono<-stepAIC(glm(ErM_ratiocell~ BIO1 + BIO12 + PHIHOX + OCSTHA, 
                         data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_ErM_Carbono<- glm(formula = ErM_ratiocell ~ BIO1 + BIO12 + OCSTHA, family = gaussian(link = "identity"), 
                              data = mayores_a_100m)

ErM_Nitrogeno<-stepAIC(glm(ErM_ratiocell~ BIO1 + BIO12 + PHIHOX + Nitrogeno, 
                           data =mayores_a_100m, gaussian(link = "identity")))


ModeloFinal_ErM_Nitrogeno<- glm(formula = ErM_ratiocell ~ BIO1 + BIO12 + PHIHOX + Nitrogeno, 
                                family = gaussian(link = "identity"), data = mayores_a_100m)

OM_Carbono<-stepAIC(glm(OM_ratiocell~ BIO1 + BIO12 + PHIHOX + OCSTHA, 
                        data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_OM_Carbono<- glm(formula = OM_ratiocell ~ BIO12 + OCSTHA, family = gaussian(link = "identity"), 
                             data = mayores_a_100m)

OM_Nitrogeno<-stepAIC(glm(OM_ratiocell~ BIO1 + BIO12  + PHIHOX + Nitrogeno, 
                          data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_OM_Nitrogeno<-glm(formula = OM_ratiocell ~ BIO1 + Nitrogeno, family = gaussian(link = "identity"), 
                              data = mayores_a_100m)

Nfix_Carbono<-stepAIC(glm(Nfix_ratiocell~ BIO1 + BIO12  + PHIHOX + OCSTHA, 
                          data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_Nfix_Carbono<-glm(formula = Nfix_ratiocell ~ BIO1 + OCSTHA, family = gaussian(link = "identity"), 
                              data = mayores_a_100m)

Nfix_Nitrogeno<-stepAIC(glm(Nfix_ratiocell~ BIO1 + BIO12  + PHIHOX + Nitrogeno, 
                            data =mayores_a_100m, gaussian(link = "identity")))


ModeloFinal_Nfix_Nitrogeno<-glm(formula = Nfix_ratiocell ~ BIO1 + Nitrogeno, family = gaussian(link = "identity"), 
                                data = mayores_a_100m)


WanNM_Nitrogeno<-stepAIC(glm(WanNm_ratiocell~ BIO1 + BIO12 + PHIHOX + Nitrogeno, 
                             data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_WanNM_Nitrogeno<-glm(formula = WanNm_ratiocell ~ BIO1 + Nitrogeno, family = gaussian(link = "identity"), 
                                 data = mayores_a_100m)


WanNM_Carbono<-stepAIC(glm(WanNm_ratiocell~ BIO1 + BIO12  + PHIHOX + OCSTHA, 
                           data =mayores_a_100m, gaussian(link = "identity")))

ModeloFinal_WanNM_Carbono<- glm(formula = WanNm_ratiocell ~ BIO1+ PHIHOX , family = gaussian(link = "identity"), 
                                data = mayores_a_100m)



#Create supplementary table of models 

tabla_glm_nitro_myco<-tab_model(ModeloFinal_AM_Nitrogeno,
                                ModeloFinal_EM_Nitrogeno,
                                ModeloFinal_ErM_Nitrogeno,
                                ModeloFinal_OM_Nitrogeno,
                                ModeloFinal_WanNM_Nitrogeno,
                                ModeloFinal_Nfix_Nitrogeno,
                                auto.label = TRUE)


tabla_glm_carbono_myco<-tab_model(ModeloFinal_AM_Carbono,
                                  ModeloFinal_EM_Carbono,
                                  ModeloFinal_ErM_Carbono,
                                  ModeloFinal_OM_Carbono,
                                  ModeloFinal_WanNM_Carbono,
                                  ModeloFinal_Nfix_Carbono,
                                  auto.label = TRUE)
#Names of variables in models

Variables<-c("Annual Mean Temperature"="BIO1","Annual Precipitation"="BIO12","Soil organic carbon stock"="OCSTHA","Nitrogen"="Nitrogeno","pH=PHIHOX")

#Graphic summary of models 

shannonGLM_plots<-jtools::plot_summs(shannonNitrogenoGLM, shannonCarbonoGLMFINAL, 
                                     scale=TRUE, inner_ci_level = .95, plot.distributions = FALSE, 
                                     coefs = c("Annual Mean Temperature"="BIO1","Annual Precipitation"="BIO12",
                                               "Soil organic carbon stock"="OCSTHA","Nitrogen"="Nitrogeno"), 
                                     model.names = c("GLM with Nitrogen","GLM with Carbono stock"))

Nitrogeno_micorrizas<-jtools::plot_summs(ModeloFinal_AM_Nitrogeno,
                                         ModeloFinal_EM_Nitrogeno,
                                         ModeloFinal_ErM_Nitrogeno,
                                         ModeloFinal_OM_Nitrogeno,
                                         ModeloFinal_WanNM_Nitrogeno,
                                         ModeloFinal_Nfix_Nitrogeno,
                                         scale=TRUE, inner_ci_level = .95, plot.distributions = FALSE, 
                                         coefs = Variables,model.names = c("Am ratio per cell", "Em ratio per cell","Erm ratio per cell","Om ratio per cell","WanNm ratio per cell","Nfixing ratio per cell"))

Carbono_stock_micorrizas<-jtools::plot_summs(ModeloFinal_AM_Carbono,
                                             ModeloFinal_EM_Carbono,
                                             ModeloFinal_ErM_Carbono,
                                             ModeloFinal_OM_Carbono,
                                             ModeloFinal_WanNM_Carbono,
                                             ModeloFinal_Nfix_Carbono,
                                             scale=TRUE, robust = TRUE, inner_ci_level = .95, plot.distributions = FALSE, 
                                             coefs = c("Annual Mean Temperature"="BIO1","Annual Precipitation"="BIO12",
                                                       "Soil organic carbon stock"="OCSTHA","Nitrogen"="Nitrogeno"),
                                             model.names = c("Am ratio per cell", "Em ratio per cell","Erm ratio per cell",
                                                             "Om ratio per cell","WanNm ratio per cell","Nfix ratio per cell"))


apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Helvetica'),
        legend.title=element_blank(), 
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size = 13))

#Figure 5 

NitrogenGLM_final<-Nitrogeno_micorrizas  + apatheme + labs(x = "\n Estimated value of the variable’s coefficient in the model \n ", y = NULL)  


CarbonGLM_final<-Carbono_stock_micorrizas  + apatheme + labs(x = NULL, y = NULL)  


shannonGLM_final<-shannonGLM_plots+ apatheme + labs(x = NULL, y = NULL) 

plot_grid(shannonGLM_final, CarbonGLM_final, NitrogenGLM_final,
          ncol = 1, nrow = 3, align = "v" , labels = c("a", "b", "c"))

Graph5<-plot_grid(shannonGLM_final, CarbonGLM_final, NitrogenGLM_final,
                  ncol = 1, nrow = 3, align = "v" , labels = c("a", "b", "c"))


############
# Figure 7 #
############

library(sp)
library(sf)
library(vegan)

#Load spatial inputs

andes<-st_read("./insumos/Andes.shp")
andes<-as(andes, "Spatial")

Colombia<-st_read("./insumos/Colombia.shp")
Colombia<-as(Colombia, "Spatial")

spat_mico<-st_read("./insumos/Spat_mico.shp")
spat_mico<-as(spat_mico,"Spatial")


# 

bb10 <- bbox(andes) #This command retrieves the bounding box of the spatial object andes and stores it in the variable bb10
cs10 <- c(1/111.325,1/111.325)*10  # This line defines the cell size for a grid.
cc10 <- bb10[, 1] + (cs10/2)  # calculate the cell center offset.
cd10 <- ceiling(diff(t(bb10))/cs10)  # calculate the number of cells per direction
grd10 <- GridTopology(cellcentre.offset=cc10, cellsize=cs10, cells.dim=cd10)#creates a grid topology with the cell center offset, cell size, and number of cells per direction.
grd10
sp_grd10 <- SpatialGridDataFrame(grd10,
                                 data=data.frame(id=1:prod(cd10)),
                                 proj4string=CRS(proj4string(andes))) #This line creates a SpatialGridDataFrame from the grid topology and some associated data.
filtrado_100<-Data10kmFiltrado13_12_20%>%filter(Elev10kmcell>101)


xy10<-as.data.frame(coordinates(sp_grd10)) #This command extracts the coordinates of the SpatialGridDataFrame and converts them into a dataframe.
xy10km<-xy10[c(Data10kmFiltrado13_12_20$id),] #This line selects the rows of the dataframe xy10 that correspond to the ids in Data10kmFiltrado13_12_20
Data10kmFiltrado13_12_20<-cbind(Data10kmFiltrado13_12_20, xy10km) #This command adds the selected coordinates to the dataframe Data10kmFiltrado13_12_20
shannon10kmEspeciesAM_EM_ERM<-c(diversity(Data10kmFiltrado13_12_20[,9:11], index = "shannon")) #This line calculates the Shannon diversity of columns 9 to 11 of Data10kmFiltrado13_12_20
Data10kmFiltrado13_12_20<-cbind(Data10kmFiltrado13_12_20,shannon10kmEspeciesAM_EM_ERM) #This command adds the calculated Shannon diversity to the dataframe Data10kmFiltrado13_12_20.

filtrado_100<-Data10kmFiltrado13_12_20%>%filter(Elev10kmcell>101)

ProporcionMyco10km<-cbind(filtrado_100[,c(45:49)]) #This line is creating a new dataframe ProporcionMyco10km that consists of columns 45 to 49 from the dataframe Data10kmFiltrado13_12_20
ProporcionMyco10km.hel100<-decostand(ProporcionMyco10km, "hellinger") #transforming the data in ProporcionMyco10km using Hellinger transformation
climatico10km100<-cbind(filtrado_100[,c(15,26)])# new dataframe climatico10km that consists of columns 15 and 26 from the dataframe Data10kmFiltrado13_12_20
edafico10km100<-cbind(filtrado_100[,c(35:43,53)])

row.names(xy10km) <- as.numeric(row.names(xy10km))

xy10km_filtrado <- xy10km[row.names(xy10km) %in% filtrado_100$id, ]

espacial10km100<-cbind(xy10km_filtrado) #This line is creating a new dataframe espacial10km that consists of the data in xy10km
ProporcionMyco10km.part100<- varpart(ProporcionMyco10km.hel100,climatico10km100,edafico10km100,espacial10km100) # This line is performing a variation partitioning analysis on the Hellinger-transformed data 
                                                                                                                # with respect to the climatic, edaphic, and spatial dataframes. The results are stored in ProporcionMyco10km.part
ProporcionMyco10km.part100 #printing the results of the variation partitioning analysis
plot(ProporcionMyco10km.part100, digits=2, Xnames= c("Climatic", "Edaphic", "Spatial")
     , bg=c("lightsteelblue3", "orange3", "plum2", "navy"))# This is creating a plot of the variation partitioning analysis results. The digits=2 argument specifies that the numbers on the plot should be rounded to 2 decimal places. 
# The Xnames argument provides the names for the different factors in the analysis. 
# the bg argument specifies the colors to be used in the plot.





