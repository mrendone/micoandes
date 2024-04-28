
###########################################
#                                         #
#           Figures                       #
#                                         #
###########################################

library(dplyr)
library(readxl)


Data10kmFiltrado13_12_20 <- read_excel("insumos/Data10kmFiltrado13_12_20.xlsx")
Data10kmFiltrado13_12_20<- Data10kmFiltrado13_12_20 %>% rename(id="IDceldas10km.IDcelda10km")
mayores_a_100m<-Data10kmFiltrado13_12_20 %>% filter(Elev10kmcell>100)

library(mgcv)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(car)

############
# Figure 1 #
############


library(raster)
library(viridis)
library(ggspatial)
library(sf)


andes<-st_read("insumos/Andes.shp")
andes<-as(andes, "Spatial")

Colombia<-st_read("insumos/Colombia.shp")
Colombia<-as(Colombia, "Spatial")

spat_mico<-st_read("insumos/Spat_mico.shp")
spat_mico<-as(spat_mico,"Spatial")

ColomabiaNew<-raster::bind(andes, Colombia) 


Colombia_sf<-st_as_sf(ColomabiaNew)
cordillera_sf<-st_as_sf(andes)



spat_mico<-st_read("insumos/Spat_mico.shp")
spat_mico<-as(spat_mico,"Spatial")

#This command retrieves the bounding box of the spatial object andes and stores it in the variable bb10
bb10 <- bbox(andes) 
# This line defines the cell size for a grid
cs10 <- c(1/111.325,1/111.325)*10  
# calculate the cell center offset
cc10 <- bb10[, 1] + (cs10/2)  
# calculate the number of cells per direction
cd10 <- ceiling(diff(t(bb10))/cs10)  
#creates a grid topology with the cell center offset, cell size, and number of cells per direction.
grd10 <- GridTopology(cellcentre.offset=cc10, cellsize=cs10, cells.dim=cd10)
grd10
#This line creates a SpatialGridDataFrame from the grid topology and some associated data.
sp_grd10 <- SpatialGridDataFrame(grd10,
                                 data=data.frame(id=1:prod(cd10)),
                                 proj4string=CRS(proj4string(andes))) 



# Convert 'sp_grd10' to a 'SpatialPolygonsDataFrame' object
spPoly <- as(sp_grd10, "SpatialPolygonsDataFrame")

# Join the data from 'Data10kmFiltrado13_12_20' to the data of 'spPoly' using 'id' as the key
spPoly@data <- full_join(spPoly@data, Data10kmFiltrado13_12_20, by="id")

# Convert 'spPoly' to an 'sf' object and remove rows with NA values
andes_sf <- st_as_sf(spPoly)
andes_sf <- na.omit(andes_sf)

# Convert 'ColomabiaNew' and 'andes' to 'sf' objects
Colombia_sf <- st_as_sf(ColomabiaNew)
cordillera_sf <- st_as_sf(andes)

# Rename the columns of 'andes_sf' for better interpretation
andes_sf <- andes_sf %>% rename(Elevation=Elev10kmcell, AM_ratio=AM_ratiocell, EM_ratio=EM_ratiocell,
                                ErM_ratio=ErM_ratiocell, OM_ratio=OM_ratiocell, WanNm_ratio=WanNm_ratiocell,
                                Nfix_ratio=Nfix_ratiocell, Shannon_Index=shannon10kmRegistros)

# Begin the construction of the plot 
Elev_C <- ggplot(data = Colombia_sf) +
  # Add a layer for the data of Colombia
  geom_sf() +
  # Add a layer for the data of the cordillera
  geom_sf(data = cordillera_sf) +
  # Add a layer for the data of the Andes, with colors based on the elevation
  geom_sf(data = andes_sf, aes(fill = Elevation, color = Elevation), lwd = 0) +
  # Add a color scale for the fill based on the "plasma" palette
  scale_fill_viridis_c(option = "plasma") +
  # Add a color scale for the lines based on the "plasma" palette
  scale_color_viridis_c(option = "plasma") +
  # Add a scale annotation at the bottom left of the plot
  annotation_scale(location = "bl", width_hint = 0.4) +
  # Apply a black and white theme to the plot
  theme_bw()

# Display the plot
Elev_C


############
# Figure 2 #
############


#command is reading a CSV file named “tabShaEnvirRaref_50.csv” that is a rarefiated data to 50 samples per cell

tabShaEnvirRaref_50<-read.csv("insumos/tabShaEnvirRaref_50.csv")

#Here it is extracting the row names from the dataframe tabShaEnvirRaref_50 and storing them in a variable called index

index<-row.names(tabShaEnvirRaref_50)
tabShaEnvirRaref_50$id<-as.numeric(index)

#Here it is selecting columns 1 and 53 from the dataframe Data10kmFiltrado13_12_20 and storing them in a variable called nitro
id_nitro<-which(Data10kmFiltrado13_12_20$id==tabShaEnvirRaref_50$X)

nitro<-subset(Data10kmFiltrado13_12_20, id %in% c(tabShaEnvirRaref_50$X))

#This command is performing an “inner join” of the dataframes tabShaEnvirRaref_50 and nitro based on the id column. The result is stored back in tabShaEnvirRaref_50.

tabShaEnvirRaref_50<-data.frame(tabShaEnvirRaref_50,nitro)
#Shannon graph and shannon rarefied graphs  vs elevation and Carbon stock

modelo_gamElev10kmcell <- gam(shannon10kmRegistros ~ s(Elev10kmcell, bs = "cs"), data = mayores_a_100m)

# Calculate confidence intervals for predictions
intervalos_confianzaElev10kmcell <- predict(modelo_gamElev10kmcell, type = "response", interval = "confidence", se.fit = TRUE)

# Extract predictions and margin of error
prediccionesElev10kmcell <- intervalos_confianzaElev10kmcell$fit
margen_errorElev10kmcell <- intervalos_confianzaElev10kmcell$se.fit

# Calculate the lower and upper limits of the confidence intervals
lwrElev10kmcell <- prediccionesElev10kmcell - 1.96 * margen_errorElev10kmcell  
uprElev10kmcell <- prediccionesElev10kmcell + 1.96 * margen_errorElev10kmcell  

# Limit the predicted values so that they are within the expected range
prediccionesElev10kmcell <- pmin(pmax(prediccionesElev10kmcell, 0), 1.277034)
lwrElev10kmcell<- pmin(pmax(lwrElev10kmcell, 0), 1.277034)
uprElev10kmcell <- pmin(pmax(uprElev10kmcell, 0), 1.277034)


# Create a dataframe with the predictions and confidence intervals
datos_filtrados_gamElev10kmcell <- data.frame(Elev10kmcell = mayores_a_100m$Elev10kmcell,
                                              predicciones= prediccionesElev10kmcell,  
                                              lwr=lwrElev10kmcell ,  
                                              upr=uprElev10kmcell )  


shannon10kmRegistrosGraphNElev10kmcell_gam <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Elev10kmcell, y = shannon10kmRegistros)) +
  geom_line(data = datos_filtrados_gamElev10kmcell, aes(x = Elev10kmcell, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +
  geom_ribbon(data = datos_filtrados_gamElev10kmcell, aes(x = Elev10kmcell, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +  # Intervalo de confianza
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Elev10kmcell", y = "Diversity of types of mycorrhizal associations")
  geom_text(data = NULL, aes(label = paste("r =", round(cor(mayores_a_100m$Elev10kmcell, mayores_a_100m$shannon10kmEspecies), 2), ", p < 2.2e-16")), x = min(mayores_a_100m$Elev10kmcell), y = max(mayores_a_100m$shannon10kmRegistros), hjust = 0, vjust = 1, size = 4, color = "black")




RAFshannonRegistrosGraphElev<-ggplot(tabShaEnvirRaref_50, aes(x = Elev10kmcell, y = mean_sha)) +
  geom_point()+
  geom_smooth( color = "#1f77b4")+
  geom_segment(aes(xend = Elev10kmcell, yend = s_p5), linetype = "dotdash", color = "grey") +
  geom_segment(aes(xend = Elev10kmcell, yend = s_p95), linetype = "dotdash", color = "grey") +
  
  geom_point(color = "black", size = 0.75)+ stat_cor(method = "pearson",cor.coef.name=c("r")) + labs(x="Elevation", y= "Diversity of types of mycorrhizal associations")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text= element_text(size=12), axis.title.x= element_text(color="black", size=15), axis.title.y= element_text(color="black", size=15))+
  ggtitle("Rarefied Data")



####################   gam model for Carbon 

modelo_gamOCSTHA <- gam(shannon10kmRegistros ~ s(OCSTHA, bs = "cs"), data = mayores_a_100m)

intervalos_confianzaOCSTHA <- predict(modelo_gamOCSTHA, type = "response", interval = "confidence", se.fit = TRUE)
prediccionesOCSTHA <- intervalos_confianzaOCSTHA$fit
margen_errorOCSTHA  <- intervalos_confianzaOCSTHA$se.fit

lwrOCSTHA <- prediccionesOCSTHA - 1.96 * margen_errorOCSTHA   
uprOCSTHA <- prediccionesOCSTHA + 1.96 * margen_errorOCSTHA   

prediccionesOCSTHA <- pmin(pmax(prediccionesOCSTHA, 0.1988719), 0.7029917)
lwrOCSTHA<- pmin(pmax(lwrOCSTHA, 0.1988719), 0.7029917)
uprOCSTHA <- pmin(pmax(uprOCSTHA, 0.1988719), 0.7029917)


datos_filtrados_gamOCSTHA <- data.frame(OCSTHA = mayores_a_100m$OCSTHA,
                                        predicciones= prediccionesOCSTHA , 
                                        lwr=lwrOCSTHA,  
                                        upr=uprOCSTHA ) 

shannon10kmRegistrosGraphNOCSTHA_gam <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = OCSTHA, y = shannon10kmRegistros)) +
  geom_line(data = datos_filtrados_gamOCSTHA, aes(x = OCSTHA, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +
  geom_ribbon(data = datos_filtrados_gamOCSTHA, aes(x = OCSTHA, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +  # Intervalo de confianza
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Soil organic carbon stock (ton/ha)", y = "Diversity of types of mycorrhizal associations")
shannon10kmRegistrosGraphNOCSTHA_gam <- shannon10kmRegistrosGraphNOCSTHA_gam +
  geom_text(data = NULL, aes(label = paste("r =", round(cor(mayores_a_100m$OCSTHA,  mayores_a_100m$shannon10kmEspecies), 2), ", p < 2.2e-16")), x = min(mayores_a_100m$OCSTHA), y = max(mayores_a_100m$shannon10kmRegistros), hjust = 0, vjust = 1, size = 4, color = "black")


RAFshannonRegistrosGraphOCSTHA<-ggplot(tabShaEnvirRaref_50, aes(x = OCSTHA, y = mean_sha)) +
  geom_point()+
  geom_smooth(color = "#1f77b4")+
  geom_segment(aes(xend = OCSTHA, yend = s_p5), linetype = "dotdash", color = "grey") +
  geom_segment(aes(xend = OCSTHA, yend = s_p95), linetype = "dotdash", color = "grey") +
  
  geom_point(color = "black", size = 0.75)+ stat_cor(method = "pearson",cor.coef.name=c("r")) + labs(x="Soil organic carbon stock (ton/ha)", y= "Diversity of types of mycorrhizal associations")+
  theme_classic()+
  theme(axis.text= element_text(size=12), axis.title.x= element_text(color="black", size=15), axis.title.y= element_text(color="black", size=15))

####################  gam model for Nitrogen

modelo_gamNitrogeno <- gam(shannon10kmRegistros ~ s(Nitrogeno, bs = "cs"), data = mayores_a_100m)


intervalos_confianzaNitrogeno <- predict(modelo_gamNitrogeno, type = "response", interval = "confidence", se.fit = TRUE)
prediccionesNitrogeno <- intervalos_confianzaNitrogeno$fit
margen_errorNitrogeno <- intervalos_confianzaNitrogeno$se.fit

lwrNitrogeno<- prediccionesNitrogeno - 1.96 * margen_errorNitrogeno   
uprNitrogeno <- margen_errorNitrogeno + 1.96 * margen_errorNitrogeno   
lwrNitrogeno<- pmin(pmax(lwrNitrogeno, 0.2384304), 0.7743406)
uprNitrogeno <- pmin(pmax(uprNitrogeno, 0.2384304), 0.7743406)


datos_filtrados_gamNitrogeno <- data.frame(Nitrogeno = mayores_a_100m$Nitrogeno,
                                           predicciones= prediccionesNitrogeno ,  
                                           lwr=lwrNitrogeno,  
                                           upr=uprNitrogeno )  


shannon10kmRegistrosGraphNNitrogeno_gam <- ggplot(data = mayores_a_100m, aes(x = Nitrogeno, y = shannon10kmRegistros)) +
  geom_point() +
  geom_smooth(color = "#1f77b4")+
  theme_classic()
theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x =  "Total Nitrogeno in soil (cg/kg)", y = "Diversity of types of mycorrhizal associations")

# Agregar coeficiente de correlación de Pearson y valor p
shannon10kmRegistrosGraphNNitrogeno_gam <- shannon10kmRegistrosGraphNNitrogeno_gam +
  geom_text(data = NULL, aes(label = paste("r =", round(cor(mayores_a_100m$Nitrogeno,  mayores_a_100m$shannon10kmEspecies), 2), ", p < 2.2e-16")), x = min(mayores_a_100m$Nitrogeno), y = max(mayores_a_100m$shannon10kmRegistros), hjust = 0, vjust = 1, size = 4, color = "black")


RAFshannonRegistrosGraphNitrogeno<-ggplot(tabShaEnvirRaref_50, aes(x = Nitrogeno, y = mean_sha)) +
  geom_point()+
  geom_smooth(method = "gam", color = "#1f77b4")+
  geom_segment(aes(xend = Nitrogeno, yend = s_p5), linetype = "dotdash", color = "grey") +
  geom_segment(aes(xend = Nitrogeno, yend = s_p95), linetype = "dotdash", color = "grey") +
  
  geom_point(color = "black", size = 0.75)+ stat_cor(method = "pearson",cor.coef.name=c("r")) + labs(x="Total Nitrogeno in soil (cg/kg)", y= "Diversity of types of mycorrhizal associations")+
  theme_classic()+
  theme(axis.text= element_text(size=12), axis.title.x= element_text(color="black", size=15), axis.title.y= element_text(color="black", size=15))

fig2<-plot_grid(shannon10kmRegistrosGraphNElev10kmcell_gam,RAFshannonRegistrosGraphElev,shannon10kmRegistrosGraphNOCSTHA_gam,RAFshannonRegistrosGraphOCSTHA
                ,shannon10kmRegistrosGraphNNitrogeno_gam,RAFshannonRegistrosGraphNitrogeno, labels=c("a","b","c","d","e","f"), ncol=2,label_size = 15)



############
# Figure 3 #
############

#  AM_ratiocell10kmGLM gam model
modelo_AM_ratiocell10kmGLM <- gam(AM_ratiocell ~ s(Elev10kmcell), data = mayores_a_100m) 

# Calculate confidence intervals for the model
intervalos_confianza_AM_ratiocell10kmGLM <- predict(modelo_AM_ratiocell10kmGLM, type = "response", interval = "confidence", se.fit = TRUE)
predicciones_AM_ratiocell10kmGLM <- intervalos_confianza_AM_ratiocell10kmGLM$fit
margen_error_AM_ratiocell10kmGLM <- intervalos_confianza_AM_ratiocell10kmGLM$se.fit
lwr_AM_ratiocell10kmGLM <- predicciones_AM_ratiocell10kmGLM - 1.96 * margen_error_AM_ratiocell10kmGLM
upr_AM_ratiocell10kmGLM <- predicciones_AM_ratiocell10kmGLM + 1.96 * margen_error_AM_ratiocell10kmGLM


predicciones_AM_ratiocell10kmGLM <- pmin(pmax(predicciones_AM_ratiocell10kmGLM, 0), 1)
lwr_AM_ratiocell10kmGLM <- pmin(pmax(lwr_AM_ratiocell10kmGLM, 0), 1)
upr_AM_ratiocell10kmGLM <- pmin(pmax(upr_AM_ratiocell10kmGLM, 0), 1)

# Create a dataframe with the predictions and confidence intervals
datos_filtrados_AM_ratiocell10kmGLM <- data.frame(Elev10kmcell = mayores_a_100m$Elev10kmcell,
                                                  predicciones = predicciones_AM_ratiocell10kmGLM,
                                                  lwr = lwr_AM_ratiocell10kmGLM,
                                                  upr = upr_AM_ratiocell10kmGLM)



grafico_AM_ratiocell10kmGLM <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Elev10kmcell, y = AM_ratiocell)) +
  geom_line(data = datos_filtrados_AM_ratiocell10kmGLM, aes(x = Elev10kmcell, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +
  geom_ribbon(data = datos_filtrados_AM_ratiocell10kmGLM, aes(x = Elev10kmcell, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Elevation", y = "Proportion of Am")+ stat_cor(method = "pearson",cor.coef.name=c("r"),label.x =1000 , label.y = 0.25)


# Model for EM_ratiocell10kmGLM
modelo_EM_ratiocell10kmGLM <- gam(EM_ratiocell ~ s(Elev10kmcell), data = mayores_a_100m) 

intervalos_confianza_EM_ratiocell10kmGLM <- predict(modelo_EM_ratiocell10kmGLM, type = "response", interval = "confidence", se.fit = TRUE)
predicciones_EM_ratiocell10kmGLM <- intervalos_confianza_EM_ratiocell10kmGLM$fit
margen_error_EM_ratiocell10kmGLM <- intervalos_confianza_EM_ratiocell10kmGLM$se.fit
lwr_EM_ratiocell10kmGLM <- predicciones_EM_ratiocell10kmGLM - 1.96 * margen_error_EM_ratiocell10kmGLM
upr_EM_ratiocell10kmGLM <- predicciones_EM_ratiocell10kmGLM + 1.96 * margen_error_EM_ratiocell10kmGLM

predicciones_EM_ratiocell10kmGLM <- pmin(pmax(predicciones_EM_ratiocell10kmGLM, 0), 1)
lwr_EM_ratiocell10kmGLM <- pmin(pmax(lwr_EM_ratiocell10kmGLM, 0), 1)
upr_EM_ratiocell10kmGLM <- pmin(pmax(upr_EM_ratiocell10kmGLM, 0), 1)

datos_filtrados_EM_ratiocell10kmGLM <- data.frame(Elev10kmcell = mayores_a_100m$Elev10kmcell,
                                                  predicciones = predicciones_EM_ratiocell10kmGLM,
                                                  lwr = lwr_EM_ratiocell10kmGLM,
                                                  upr = upr_EM_ratiocell10kmGLM)



grafico_EM_ratiocell10kmGLM <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Elev10kmcell, y = EM_ratiocell)) +
  geom_line(data = datos_filtrados_EM_ratiocell10kmGLM, aes(x = Elev10kmcell, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +
  geom_ribbon(data = datos_filtrados_EM_ratiocell10kmGLM, aes(x = Elev10kmcell, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Elevation", y = "Proportion of Em")+ stat_cor(method = "pearson",cor.coef.name=c("r"),label.x =1000 , label.y = 0.25)


# Model for ErM_ratiocell10kmGLM
modelo_ErM_ratiocell10kmGLM <- gam(ErM_ratiocell ~ s(Elev10kmcell), data = mayores_a_100m) 

intervalos_confianza_ErM_ratiocell10kmGLM <- predict(modelo_ErM_ratiocell10kmGLM, type = "response", interval = "confidence", se.fit = TRUE)
predicciones_ErM_ratiocell10kmGLM <- intervalos_confianza_ErM_ratiocell10kmGLM$fit
margen_error_ErM_ratiocell10kmGLM <- intervalos_confianza_ErM_ratiocell10kmGLM$se.fit
lwr_ErM_ratiocell10kmGLM <- predicciones_ErM_ratiocell10kmGLM - 1.96 * margen_error_ErM_ratiocell10kmGLM
upr_ErM_ratiocell10kmGLM <- predicciones_ErM_ratiocell10kmGLM + 1.96 * margen_error_ErM_ratiocell10kmGLM

predicciones_ErM_ratiocell10kmGLM <- pmin(pmax(predicciones_ErM_ratiocell10kmGLM, 0), 1)
lwr_ErM_ratiocell10kmGLM <- pmin(pmax(lwr_ErM_ratiocell10kmGLM, 0), 1)
upr_ErM_ratiocell10kmGLM <- pmin(pmax(upr_ErM_ratiocell10kmGLM, 0), 1)

datos_filtrados_ErM_ratiocell10kmGLM <- data.frame(Elev10kmcell = mayores_a_100m$Elev10kmcell,
                                                   predicciones = predicciones_ErM_ratiocell10kmGLM,
                                                   lwr = lwr_ErM_ratiocell10kmGLM,
                                                   upr = upr_ErM_ratiocell10kmGLM)



grafico_ErM_ratiocell10kmGLM <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Elev10kmcell, y = ErM_ratiocell)) +
  geom_line(data = datos_filtrados_ErM_ratiocell10kmGLM, aes(x = Elev10kmcell, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +
  geom_ribbon(data = datos_filtrados_ErM_ratiocell10kmGLM, aes(x = Elev10kmcell, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Elevation", y = "Proportion of Erm")+ stat_cor(method = "pearson",cor.coef.name=c("r"),label.x =1000 , label.y = 0.25)


# Model for OM_ratiocell10kmGLM
modelo_OM_ratiocell10kmGLM <- gam(OM_ratiocell ~ s(Elev10kmcell), data = mayores_a_100m) 

intervalos_confianza_OM_ratiocell10kmGLM <- predict(modelo_OM_ratiocell10kmGLM, type = "response", interval = "confidence", se.fit = TRUE)
predicciones_OM_ratiocell10kmGLM <- intervalos_confianza_OM_ratiocell10kmGLM$fit
margen_error_OM_ratiocell10kmGLM <- intervalos_confianza_OM_ratiocell10kmGLM$se.fit
lwr_OM_ratiocell10kmGLM <- predicciones_OM_ratiocell10kmGLM - 1.96 * margen_error_OM_ratiocell10kmGLM
upr_OM_ratiocell10kmGLM <- predicciones_OM_ratiocell10kmGLM + 1.96 * margen_error_OM_ratiocell10kmGLM

predicciones_OM_ratiocell10kmGLM <- pmin(pmax(predicciones_OM_ratiocell10kmGLM, 0), 1)
lwr_OM_ratiocell10kmGLM <- pmin(pmax(lwr_OM_ratiocell10kmGLM, 0), 1)
upr_OM_ratiocell10kmGLM <- pmin(pmax(upr_OM_ratiocell10kmGLM, 0), 1)

datos_filtrados_OM_ratiocell10kmGLM <- data.frame(Elev10kmcell = mayores_a_100m$Elev10kmcell,
                                                  predicciones = predicciones_OM_ratiocell10kmGLM,
                                                  lwr = lwr_OM_ratiocell10kmGLM,
                                                  upr = upr_OM_ratiocell10kmGLM)



grafico_OM_ratiocell10kmGLM <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Elev10kmcell, y = OM_ratiocell)) +
  geom_line(data = datos_filtrados_OM_ratiocell10kmGLM, aes(x = Elev10kmcell, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +
  geom_ribbon(data = datos_filtrados_OM_ratiocell10kmGLM, aes(x = Elev10kmcell, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Elevation", y = "Proportion of Om")+ stat_cor(method = "pearson",cor.coef.name=c("r"),label.x =1000 , label.y = 0.25)


# Modelo for  WanNm_ratiocell10kmGLM
modelo_WanNm_ratiocell10kmGLM <- gam(WanNm_ratiocell ~ s(Elev10kmcell), data = mayores_a_100m) 

intervalos_confianza_WanNm_ratiocell10kmGLM <- predict(modelo_WanNm_ratiocell10kmGLM, type = "response", interval = "confidence", se.fit = TRUE)
predicciones_WanNm_ratiocell10kmGLM <- intervalos_confianza_WanNm_ratiocell10kmGLM$fit
margen_error_WanNm_ratiocell10kmGLM <- intervalos_confianza_WanNm_ratiocell10kmGLM$se.fit
lwr_WanNm_ratiocell10kmGLM <- predicciones_WanNm_ratiocell10kmGLM - 1.96 * margen_error_WanNm_ratiocell10kmGLM
upr_WanNm_ratiocell10kmGLM <- predicciones_WanNm_ratiocell10kmGLM + 1.96 * margen_error_WanNm_ratiocell10kmGLM

predicciones_WanNm_ratiocell10kmGLM <- pmin(pmax(predicciones_WanNm_ratiocell10kmGLM, 0), 1)
lwr_WanNm_ratiocell10kmGLM <- pmin(pmax(lwr_WanNm_ratiocell10kmGLM, 0), 1)
upr_WanNm_ratiocell10kmGLM <- pmin(pmax(upr_WanNm_ratiocell10kmGLM, 0), 1)

datos_filtrados_WanNm_ratiocell10kmGLM <- data.frame(Elev10kmcell = mayores_a_100m$Elev10kmcell,
                                                     predicciones = predicciones_WanNm_ratiocell10kmGLM,
                                                     lwr = lwr_WanNm_ratiocell10kmGLM,
                                                     upr = upr_WanNm_ratiocell10kmGLM)



grafico_WanNm_ratiocell10kmGLM <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Elev10kmcell, y =  WanNm_ratiocell)) +
  geom_line(data = datos_filtrados_WanNm_ratiocell10kmGLM, aes(x = Elev10kmcell, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +
  geom_ribbon(data = datos_filtrados_WanNm_ratiocell10kmGLM, aes(x = Elev10kmcell, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Elevation", y = "Proportion of WanNm")+ stat_cor(method = "pearson",cor.coef.name=c("r"),label.x =1000 , label.y = 0.25)



# Model for  Nfix_ratiocell10kmGLM
modelo_Nfix_ratiocell10kmGLM <- gam( Nfix_ratiocell ~ s(Elev10kmcell), data = mayores_a_100m) 

intervalos_confianza_Nfix_ratiocell10kmGLM <- predict(modelo_Nfix_ratiocell10kmGLM, type = "response", interval = "confidence", se.fit = TRUE)
predicciones_Nfix_ratiocell10kmGLM <- intervalos_confianza_Nfix_ratiocell10kmGLM$fit
margen_error_Nfix_ratiocell10kmGLM <- intervalos_confianza_Nfix_ratiocell10kmGLM$se.fit
lwr_Nfix_ratiocell10kmGLM <- predicciones_Nfix_ratiocell10kmGLM - 1.96 * margen_error_Nfix_ratiocell10kmGLM
upr_Nfix_ratiocell10kmGLM <- predicciones_Nfix_ratiocell10kmGLM + 1.96 * margen_error_Nfix_ratiocell10kmGLM

predicciones_Nfix_ratiocell10kmGLM <- pmin(pmax(predicciones_Nfix_ratiocell10kmGLM, 0), 1)
lwr_Nfix_ratiocell10kmGLM <- pmin(pmax(lwr_Nfix_ratiocell10kmGLM, 0), 1)
upr_Nfix_ratiocell10kmGLM <- pmin(pmax(upr_Nfix_ratiocell10kmGLM, 0), 1)

datos_filtrados_Nfix_ratiocell10kmGLM <- data.frame(Elev10kmcell = mayores_a_100m$Elev10kmcell,
                                                    predicciones = predicciones_Nfix_ratiocell10kmGLM,
                                                    lwr = lwr_Nfix_ratiocell10kmGLM,
                                                    upr = upr_Nfix_ratiocell10kmGLM)



grafico_Nfix_ratiocell10kmGLM <- ggplot() +
  geom_point(data = mayores_a_100m, aes(x = Elev10kmcell, y =  Nfix_ratiocell)) +
  geom_line(data = datos_filtrados_Nfix_ratiocell10kmGLM, aes(x = Elev10kmcell, y = predicciones), color = "#1f77b4", linetype = "solid", size = 1.3) +
  geom_ribbon(data = datos_filtrados_Nfix_ratiocell10kmGLM, aes(x = Elev10kmcell, ymin = lwr, ymax = upr), fill = "lightgrey", alpha = 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title.x = element_text(color = "black", size = 13), axis.title.y = element_text(color = "black", size = 13)) +
  labs(x = "Elevation", y = "Proportion of Nfix")+ stat_cor(method = "pearson",cor.coef.name=c("r"),label.x =1000 , label.y = 0.25)

grap3<-plot_grid(grafico_AM_ratiocell10kmGLM, grafico_EM_ratiocell10kmGLM, grafico_ErM_ratiocell10kmGLM,  grafico_OM_ratiocell10kmGLM, 
                 grafico_WanNm_ratiocell10kmGLM, grafico_Nfix_ratiocell10kmGLM, labels=c("a","b","c","d","e","f"))


#############
#  Figure 4 #
#############

# The process is the same as for Figure 1, we only change the object to fill,
# in this case we are going to use the Shannon index.

Shannon_colombia<-ggplot(data=Colombia_sf)+ geom_sf()+
  geom_sf(data = cordillera_sf)+
  geom_sf(data = andes_sf,aes(fill = Shannon_Index, color=Shannon_Index), lwd=0) +
  scale_fill_viridis_c(option = "plasma")+
  scale_color_viridis_c(option = "plasma") +
  annotation_scale(location = "bl", width_hint = 0.4)+
  theme_bw()



#######################################
#  Figure 5 and supplementary tables  #
#######################################

library(MASS)
library(car)
library(sjPlot)

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


shannonNitrogenoGLM <- glm(shannon10kmRegistros ~ BIO1 + BIO12 + Nitrogeno,
                           data = mayores_a_100m,
                           gaussian(link = "identity"))

shannonCarbonoGLMFINAL <- glm(shannon10kmRegistros ~ BIO1 + BIO12 + OCSTHA,
                              data = mayores_a_100m, family = gaussian(link = "identity"))

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


NitrogenGLM_final<-Nitrogeno_micorrizas  + apatheme + labs(x = "\n Estimated value of the variable’s coefficient in the model \n ", y = NULL)  


CarbonGLM_final<-Carbono_stock_micorrizas  + apatheme + labs(x = NULL, y = NULL)  


shannonGLM_final<-shannonGLM_plots+ apatheme + labs(x = NULL, y = NULL) 

cowplot::plot_grid(shannonGLM_final, CarbonGLM_final, NitrogenGLM_final,
                   ncol = 1, nrow = 3, align = "v" , labels = c("a", "b", "c"))


Graph5<-cowplot::plot_grid(shannonGLM_final, CarbonGLM_final, NitrogenGLM_final,
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
#This command retrieves the bounding box of the spatial object andes and stores it in the variable bb10
bb10 <- bbox(andes) 
# This line defines the cell size for a grid.
cs10 <- c(1/111.325,1/111.325)*10  
# calculate the cell center offset.
cc10 <- bb10[, 1] + (cs10/2) 
# calculate the number of cells per direction
cd10 <- ceiling(diff(t(bb10))/cs10) 
#creates a grid topology with the cell center offset, cell size, and number of cells per direction.
grd10 <- GridTopology(cellcentre.offset=cc10, cellsize=cs10, cells.dim=cd10)
grd10
#This line creates a SpatialGridDataFrame from the grid topology and some associated data.
sp_grd10 <- SpatialGridDataFrame(grd10,
                                 data=data.frame(id=1:prod(cd10)),
                                 proj4string=CRS(proj4string(andes))) 
filtrado_100<-Data10kmFiltrado13_12_20%>%filter(Elev10kmcell>101)

#This command extracts the coordinates of the SpatialGridDataFrame and converts them into a dataframe.
xy10<-as.data.frame(coordinates(sp_grd10))
#This line selects the rows of the dataframe xy10 that correspond to the ids in Data10kmFiltrado13_12_20
xy10km<-xy10[c(Data10kmFiltrado13_12_20$id),] 
#This command adds the selected coordinates to the dataframe Data10kmFiltrado13_12_20
Data10kmFiltrado13_12_20<-cbind(Data10kmFiltrado13_12_20, xy10km) 
#This line calculates the Shannon diversity of columns 9 to 11 of Data10kmFiltrado13_12_20
shannon10kmEspeciesAM_EM_ERM<-c(diversity(Data10kmFiltrado13_12_20[,9:11], index = "shannon")) 
#This command adds the calculated Shannon diversity to the dataframe Data10kmFiltrado13_12_20.
Data10kmFiltrado13_12_20<-cbind(Data10kmFiltrado13_12_20,shannon10kmEspeciesAM_EM_ERM) 

filtrado_100<-Data10kmFiltrado13_12_20%>%filter(Elev10kmcell>101)
#This line is creating a new dataframe ProporcionMyco10km that consists of columns 45 to 49 from the dataframe Data10kmFiltrado13_12_20
ProporcionMyco10km<-cbind(filtrado_100[,c(45:49)]) 
#transforming the data in ProporcionMyco10km using Hellinger transformation
ProporcionMyco10km.hel100<-decostand(ProporcionMyco10km, "hellinger") 
# new dataframe climatico10km that consists of columns 15 and 26 from the dataframe Data10kmFiltrado13_12_20
climatico10km100<-cbind(filtrado_100[,c(15,26)])
edafico10km100<-cbind(filtrado_100[,c(35:43,53)])

row.names(xy10km) <- as.numeric(row.names(xy10km))

xy10km_filtrado <- xy10km[row.names(xy10km) %in% filtrado_100$id, ]

#This line is creating a new dataframe espacial10km that consists of the data in xy10km
espacial10km100<-cbind(xy10km_filtrado) 
# This line is performing a variation partitioning analysis on the Hellinger-transformed data 
# with respect to the climatic, edaphic, and spatial dataframes. The results are stored in ProporcionMyco10km.part
ProporcionMyco10km.part100<- varpart(ProporcionMyco10km.hel100,climatico10km100,edafico10km100,espacial10km100) 

ProporcionMyco10km.part100 

# This is creating a plot of the variation partitioning analysis results. The digits=2 argument specifies that the numbers on the plot should be rounded to 2 decimal places. 
# The Xnames argument provides the names for the different factors in the analysis. 
# the bg argument specifies the colors to be used in the plot.
plot(ProporcionMyco10km.part100, digits=2, Xnames= c("Climatic", "Edaphic", "Spatial")
     , bg=c("lightsteelblue3", "orange3", "plum2", "navy"))




