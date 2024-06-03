#################################
#    Models                     #
#################################

library(Metrics)
library(readxl)
Data10kmFiltrado13_12_20 <-read_excel("../insumos/Data10kmFiltrado13_12_20.xlsx")

library(dplyr)
mayores_a_100m<-Data10kmFiltrado13_12_20 %>% filter(Elev10kmcell>100)


mayores_a_100m<-Data10kmFiltrado13_12_20 %>% filter(Elev10kmcell>100)

Data10kmFiltrado13_12_20<- Data10kmFiltrado13_12_20 %>% rename(id="IDceldas10km.IDcelda10km")
Data10kmFiltrado_101 <-Data10kmFiltrado13_12_20%>%filter(Elev10kmcell>101)

mayores_a_100m<-mutate(mayores_a_100m, Cuadrado_BIO1=(BIO1)**2)
mayores_a_100m<-mutate(mayores_a_100m, Cuadrado_BIO12=(BIO12)**2)
mayores_a_100m<-mutate(mayores_a_100m, Cuadrado_AWCh1=(AWCh1)**2)
mayores_a_100m<-mutate(mayores_a_100m, Cuadrado_PHIHOX=(PHIHOX)**2)
mayores_a_100m<-mutate(mayores_a_100m, Cuadrado_Nitrogeno=(Nitrogeno)**2)
mayores_a_100m<-mutate(mayores_a_100m, Cuadrado_OCSTHA=(OCSTHA)**2)



Data10kmFiltrado_101$BIO12_log <- log(Data10kmFiltrado_101$BIO12)
Data10kmFiltrado_101$PHIHOX_log <- log(Data10kmFiltrado_101$PHIHOX)
Data10kmFiltrado_101$Nitrogeno_log <- log(Data10kmFiltrado_101$Nitrogeno)

# Load necessary libraries
library(gtsummary)
library(car)
library(MASS)

# Fit the initial model with quadratic terms
shannonCarbonoGLM_cuadrado <- glm(shannon10kmRegistros ~ BIO1 + BIO12 + AWCh1 +
                                    PHIHOX + OCSTHA + Cuadrado_BIO1 + Cuadrado_BIO12 +
                                    Cuadrado_AWCh1 + Cuadrado_PHIHOX + Cuadrado_OCSTHA,
                                  data = mayores_a_100m, family = gaussian(link = "identity"))

# Perform stepwise AIC selection for the initial model
shannonCarbonoGLM_cuadradoAIC <- stepAIC(shannonCarbonoGLM_cuadrado, trace = TRUE, direction = "backward")

# Fit a simplified model with only significant quadratic terms
shannonCarbonoGLM_cuadrado2 <- glm(shannon10kmRegistros ~ Cuadrado_BIO1 + Cuadrado_BIO12 +
                                     Cuadrado_AWCh1 + Cuadrado_PHIHOX + Cuadrado_OCSTHA,
                                   data = mayores_a_100m, family = gaussian(link = "identity"))

# Perform stepwise AIC selection for the simplified model
shannonCarbonoGLM_cuadrado2AIC <- stepAIC(shannonCarbonoGLM_cuadrado2, trace = TRUE, direction = "backward")

# Fit a final model without quadratic terms
shannonCarbonoGLM <- glm(shannon10kmRegistros ~ BIO1 + BIO12 + AWCh1 + PHIHOX + OCSTHA,
                         data = mayores_a_100m, family = gaussian(link = "identity"))

# Perform stepwise AIC selection for the final model
shannonCarbonoGLM_AIC <- stepAIC(shannonCarbonoGLM, trace = TRUE, direction = "backward")

# Fit a simplified final model with only significant predictors
shannonCarbonoGLMFINAL <- glm(shannon10kmRegistros ~ BIO1 + BIO12 + OCSTHA,
                              data = mayores_a_100m, family = gaussian(link = "identity"))

# Extracting coefficient names as strings
N_glmC <- paste(names(coef(shannonCarbonoGLM)), collapse = ", ")
N_glmcuadraC <- paste(names(coef(shannonCarbonoGLM_cuadradoAIC)), collapse = ", ")
N_glmcuadra2C <- paste(names(coef(shannonCarbonoGLM_cuadrado2AIC)), collapse = ", ")

# Defining model names
Modelo_C <- c("Model 1: shannonCarbonoGLM",
              "Model 2: shannonCarbonoGLM_Mixto_cuadrado",
              "Model 3: shannonCarbonoGLM_cuadrado")

# Calculating root mean square error (RMSE)
require(Metrics)
Rmse_ShannonCarbono <- c(rmse(Data10kmFiltrado_101$shannon10kmRegistros, predict(shannonCarbonoGLMFINAL)),
                         rmse(Data10kmFiltrado_101$shannon10kmRegistros, predict(shannonCarbonoGLM_cuadradoAIC)),
                         rmse(Data10kmFiltrado_101$shannon10kmRegistros, predict(shannonCarbonoGLM_cuadrado2AIC)))

# Collecting variable names for each model
Variables_C_AIC <- c(N_glmC, N_glmcuadraC, N_glmcuadra2C)

# Calculating AIC values for each model
AIC_C <- c(shannonCarbonoGLM$aic, 
           shannonCarbonoGLM_cuadradoAIC$aic, 
           shannonCarbonoGLM_cuadrado2AIC$aic)

# Combining all information into a data frame for viewing
result_df <- data.frame(Modelo_C, Variables_C_AIC, Rmse_ShannonCarbono, AIC_C)
View(result_df)


####################################
# coefficients of Shannon's models #
####################################


library(QuantPsyc)
library(glmmTMB)
library(MASS)
library(car)
library(gtsummary)

# Adjusting AM_ratiocell values for perfect overlap
Data10kmFiltrado_101$AM_ratiocell <- replace(Data10kmFiltrado_101$AM_ratiocell, 
                                             Data10kmFiltrado_101$AM_ratiocell == 1, 
                                             Data10kmFiltrado_101$AM_ratiocell[Data10kmFiltrado_101$AM_ratiocell == 1] - 0.00001)

# GLM Model for AM_ratiocell with Nitrogen
Modelo_final_AM_ratiocellNitrogenoGLM_AIC <- glm(AM_ratiocell ~ BIO1 + BIO12 + Nitrogeno,
                                                 data = Data10kmFiltrado_101, 
                                                 family = gaussian(link = "identity")) 


Modelo_final_AM_ratiocellCarbonoGLM<-glm(AM_ratiocell~ BIO1 + BIO12 + OCSTHA,
                                         data = Data10kmFiltrado_101, family=gaussian(link = "identity"))


Modelo_final_AM_ratiocellCarbonoZIB<-glmmTMB(AM_ratiocell ~ BIO1 + BIO12_log  + OCSTHA,
                                             ziformula=~1, 
                                             family=beta_family(), 
                                             data = Data10kmFiltrado_101)

# Boxplots to visualize outliers
boxplot(Data10kmFiltrado_101$BIO1, main = "BIO1")
boxplot(Data10kmFiltrado_101$BIO12, main = "BIO12")
boxplot(Data10kmFiltrado_101$AWCh1, main = "AWCh1")
boxplot(Data10kmFiltrado_101$PHIHOX, main = "PHIHOX")
boxplot(Data10kmFiltrado_101$Nitrogeno, main = "Nitrogeno")

# Detecting outliers
BIO12_outliers <- boxplot.stats(Data10kmFiltrado_101$BIO12)$out
PHIHOX_outliers <- boxplot.stats(Data10kmFiltrado_101$PHIHOX)$out
Nitrogeno_outliers <- boxplot.stats(Data10kmFiltrado_101$Nitrogeno)$out

# Log-transforming variables
Data10kmFiltrado_101$BIO12_log <- log(Data10kmFiltrado_101$BIO12)
Data10kmFiltrado_101$PHIHOX_log <- log(Data10kmFiltrado_101$PHIHOX)
Data10kmFiltrado_101$Nitrogeno_log <- log(Data10kmFiltrado_101$Nitrogeno)

# ZIB Model for AM_ratiocell with Nitrogen (including zero-inflation)
Modelo_final_AM_ratiocellNitrogenoZIB_AIC <- glmmTMB(AM_ratiocell ~ BIO1 + BIO12_log + Nitrogeno_log,
                                                     ziformula = ~1, 
                                                     family = beta_family(), 
                                                     data = Data10kmFiltrado_101)

# ZIB Model for AM_ratiocell with Nitrogen (without zero-inflation)
Modelo_final_AM_ratiocellNitrogenoB_AIC <- glmmTMB(AM_ratiocell ~ BIO1 + BIO12_log + Nitrogeno_log,
                                                   family = beta_family(), 
                                                   data = Data10kmFiltrado_101)

# Defining and fitting the GLM model for nitrogen
shannonNitrogenoGLM <- glm(shannon10kmRegistros ~ BIO1 + BIO12 + Nitrogeno,
                           data = mayores_a_100m,
                           gaussian(link = "identity"))

# Calculating standardized coefficients for nitrogen
standardized.coefficientsNitrogeno <- lm.beta(shannonNitrogenoGLM)

# Printing standardized coefficients for nitrogen
print(standardized.coefficientsNitrogeno)

# Carbon
# Defining and fitting the final GLM model for carbon
shannonCarbonoGLMFINAL <- glm(shannon10kmRegistros ~ BIO1 + BIO12 + OCSTHA,
                              data = mayores_a_100m,
                              gaussian(link = "identity"))

# Calculating standardized coefficients for carbon
standardized.coefficientsCarbono <- lm.beta(shannonCarbonoGLMFINAL)

# Printing standardized coefficients for carbon
print(standardized.coefficientsCarbono)


# Nitrogen Proportion of EcM

## GLM Model for Nitrogen Proportion of EcM
Modelo_Final_EM_ratiocell_NitrogenoGLM <- glm(EM_ratiocell ~ BIO1 + BIO12 + Nitrogeno,
                                              data = Data10kmFiltrado_101,
                                              family = gaussian(link = "identity"))

## ZIB Model for Nitrogen Proportion of EcM
Modelo_Final_EM_ratiocell_NitrogenoZIB <- glmmTMB(EM_ratiocell ~ BIO1 + BIO12_log + Nitrogeno_log,
                                                  ziformula = ~1,
                                                  family = beta_family(),
                                                  data = Data10kmFiltrado_101)

### Model Comparison Using AIC
aic_EM_ratiocell_NitrogenoGLM <- AIC(Modelo_Final_EM_ratiocell_NitrogenoGLM)
aic_EM_ratiocell_NitrogenoZIB <- AIC(Modelo_Final_EM_ratiocell_NitrogenoZIB)

print(paste("AIC EM_ratiocell Nitrogeno GLM: ", aic_EM_ratiocell_NitrogenoGLM))
print(paste("AIC EM_ratiocell Nitrogeno ZIB: ", aic_EM_ratiocell_NitrogenoZIB))


# Carbon Proportion of EcM

## GLM Model for Carbon Proportion of EcM
Modelo_final_EM_ratiocell_CarbonoGLM <- glm(EM_ratiocell ~ BIO1 + BIO12 + OCSTHA,
                                            data = Data10kmFiltrado_101,
                                            family = gaussian(link = "identity"))

## ZIB Model for Carbon Proportion of EcM
Modelo_final_EM_ratiocell_CarbonoZIB <- glmmTMB(EM_ratiocell ~ BIO1 + BIO12_log + OCSTHA,
                                                ziformula = ~1,
                                                family = beta_family(),
                                                data = Data10kmFiltrado_101)

### Model Comparison Using AIC
aic_EM_ratiocell_CarbonoGLM <- AIC(Modelo_final_EM_ratiocell_CarbonoGLM)
aic_EM_ratiocell_CarbonoGLMZIB <- AIC(Modelo_final_EM_ratiocell_CarbonoZIB)

print(paste("AIC EM_ratiocell Carbono GLM: ", aic_EM_ratiocell_CarbonoGLM))
print(paste("AIC EM_ratiocell Carbono ZIB: ", aic_EM_ratiocell_CarbonoGLMZIB))


# ErM Proportion

## GLM Model for Nitrogen Proportion of ErM
Modelo_Final_ErM_ratiocell_NitrogenoGLM <- glm(ErM_ratiocell ~ BIO1 + BIO12 + Nitrogeno,
                                               data = Data10kmFiltrado_101,
                                               family = gaussian(link = "identity"))

## ZIB Model for Nitrogen Proportion of ErM
Modelo_Final_ErM_ratiocell_NitrogenoZIB <- glmmTMB(ErM_ratiocell ~ BIO1 + BIO12_log + Nitrogeno_log,
                                                   ziformula = ~1,
                                                   family = beta_family(),
                                                   data = Data10kmFiltrado_101)

### Model Comparison Using AIC
aic_ErM_ratiocell_NitrogenoGLM <- AIC(Modelo_Final_ErM_ratiocell_NitrogenoGLM)
aic_ErM_ratiocell_NitrogenoZIB <- AIC(Modelo_Final_ErM_ratiocell_NitrogenoZIB)

print(paste("AIC ErM_ratiocell Nitrogeno GLM: ", aic_ErM_ratiocell_NitrogenoGLM))
print(paste("AIC ErM_ratiocell Nitrogeno ZIB: ", aic_ErM_ratiocell_NitrogenoZIB))


## GLM Model for Carbon Proportion of ErM
Modelo_Final_ErM_ratiocell_CarbonoGLM <- glm(ErM_ratiocell ~ BIO1 + BIO12 + OCSTHA,
                                             data = Data10kmFiltrado_101,
                                             family = gaussian(link = "identity"))

## ZIB Model for Carbon Proportion of ErM
Modelo_Final_ErM_ratiocell_CarbonoZIB <- glmmTMB(ErM_ratiocell ~ BIO1 + BIO12_log + OCSTHA,
                                                 ziformula = ~1,
                                                 family = beta_family(),
                                                 data = Data10kmFiltrado_101)

### Model Comparison Using AIC
aic_ErM_ratiocell_CarbonoGLM <- AIC(Modelo_Final_ErM_ratiocell_CarbonoGLM)
aic_ErM_ratiocell_CarbonoZIB <- AIC(Modelo_Final_ErM_ratiocell_CarbonoZIB)

print(paste("AIC ErM_ratiocell Carbono GLM: ", aic_ErM_ratiocell_CarbonoGLM))
print(paste("AIC ErM_ratiocell Carbono ZIB: ", aic_ErM_ratiocell_CarbonoZIB))


# OM Proportion

## GLM Model for Nitrogen Proportion of OM
Modelo_Final_OM_ratiocell_NitrogenoGLM <- glm(OM_ratiocell ~ BIO1 + BIO12 + Nitrogeno,
                                              data = Data10kmFiltrado_101,
                                              family = gaussian(link = "identity"))

## ZIB Model for Nitrogen Proportion of OM
Modelo_Final_OM_ratiocell_NitrogenoZIB <- glmmTMB(OM_ratiocell ~ BIO1 + BIO12_log + Nitrogeno_log,
                                                  ziformula = ~1,
                                                  family = beta_family(),
                                                  data = Data10kmFiltrado_101)






Modelo_Final_OM_ratiocell_CarbonoGLM<-glm(OM_ratiocell~ BIO1 + BIO12 + OCSTHA,
                                          data = Data10kmFiltrado_101, family=gaussian(link = "identity")) 


Modelo_Final_OM_ratiocell_CarbonoZIB<- glmmTMB(OM_ratiocell ~ BIO1 + BIO12_log + OCSTHA,
                                               ziformula=~1, 
                                               family=beta_family(), 
                                               data = Data10kmFiltrado_101)

aic_OM_ratiocell_CarbonoGLM<- AIC(Modelo_Final_OM_ratiocell_CarbonoGLM)
aic_OM_ratiocell_CarbonoZIB <- AIC(Modelo_Final_OM_ratiocell_CarbonoZIB)

print(paste("AIC OM_ratiocell Carbono GLM: ", aic_OM_ratiocell_CarbonoGLM))
print(paste("AIC OM_ratiocell Carbono ZIB: ", aic_OM_ratiocell_CarbonoZIB))


Modelo_Final_WanNm_ratiocell_NitrogenoGLM<-glm(WanNm_ratiocell~ BIO1 + BIO12 + Nitrogeno,
                                               data = Data10kmFiltrado_101, family=gaussian(link = "identity")) 



Modelo_Final_WanNm_ratiocell_NitrogenoZIB<- glmmTMB(WanNm_ratiocell ~ BIO1 + BIO12_log + Nitrogeno_log,
                                                    ziformula=~1, 
                                                    family=beta_family(), 
                                                    data = Data10kmFiltrado_101)

aic_WanNm_ratiocell_NitrogenoGLM<- AIC(Modelo_Final_WanNm_ratiocell_NitrogenoGLM)
aic_WanNm_ratiocelll_NitrogenoZIB <- AIC(Modelo_Final_WanNm_ratiocell_NitrogenoZIB)

print(paste("AIC ErM_ratiocell Nitrogeno GLM: ", aic_WanNm_ratiocell_NitrogenoGLM))
print(paste("AIC ErM_ratiocell Nitrogeno ZIB: ", aic_WanNm_ratiocelll_NitrogenoZIB))

Modelo_Final_WanNm_ratiocell_CarbonoGLM<-glm(WanNm_ratiocell~ BIO1 + BIO12 + OCSTHA,
                                             data = Data10kmFiltrado_101, family=gaussian(link = "identity")) 



Modelo_Final_WanNm_ratiocell_CarbonoZIB<- glmmTMB(WanNm_ratiocell ~ BIO1 + BIO12_log + OCSTHA,
                                                  ziformula=~1, 
                                                  family=beta_family(), 
                                                  data = Data10kmFiltrado_101)



aic_WanNm_ratiocell_CarbonoGLM<- AIC(Modelo_Final_WanNm_ratiocell_CarbonoGLM)
aic_WanNm_ratiocell_CarbonoZIB <- AIC(Modelo_Final_WanNm_ratiocell_CarbonoZIB)

print(paste("AIC WanNm_ratiocell Carbono GLM: ", aic_WanNm_ratiocell_CarbonoGLM))
print(paste("AIC WanNm_ratiocell Carbono ZIB: ", aic_WanNm_ratiocell_CarbonoZIB))


Modelo_Final_Nfix_ratiocell_NitrogenoGLM<-glm(Nfix_ratiocell~ BIO1 + BIO12 + Nitrogeno,
                                              data = Data10kmFiltrado_101, family=gaussian(link = "identity")) 



Modelo_Final_Nfix_ratiocell_NitrogenoZIB<- glmmTMB(Nfix_ratiocell ~ BIO1 + BIO12_log + Nitrogeno_log,
                                                   ziformula=~1, 
                                                   family=beta_family(), 
                                                   data = Data10kmFiltrado_101)

aic_Nfix_ratiocell_NitrogenooGLM<- AIC(Modelo_Final_Nfix_ratiocell_NitrogenoGLM)
aic_Nfix_ratiocell_NitrogenoZIB <- AIC(Modelo_Final_Nfix_ratiocell_NitrogenoZIB)

print(paste("AIC Nfix_ratiocell Nitrogeno GLM: ", aic_Nfix_ratiocell_NitrogenooGLM))
print(paste("AIC Nfix_ratiocell Nitrogeno ZIB: ", aic_Nfix_ratiocell_NitrogenoZIB))

Modelo_Final_Nfix_ratiocell_OCSTHAGLM<-glm(Nfix_ratiocell~ BIO1 + BIO12 + OCSTHA,
                                           data = Data10kmFiltrado_101, family=gaussian(link = "identity")) 



Modelo_Final_Nfix_ratiocell_OCSTHAZIB<- glmmTMB(Nfix_ratiocell ~ BIO1 + BIO12_log + OCSTHA,
                                                ziformula=~1, 
                                                family=beta_family(), 
                                                data = Data10kmFiltrado_101)

aic_Nfix_ratiocell_CarbonoGLM<- AIC(Modelo_Final_Nfix_ratiocell_OCSTHAGLM)
aic_Nfix_ratiocell_CarbonoZIB <- AIC(Modelo_Final_Nfix_ratiocell_OCSTHAZIB)

print(paste("AIC Nfix_ratiocell Carbono GLM: ", aic_Nfix_ratiocell_CarbonoGLM))
print(paste("AIC Nfix_ratiocell Carbono ZIB: ", aic_Nfix_ratiocell_CarbonoZIB))

#Effect size proporción micorrizas

library(QuantPsyc)

#Nfix
coefficients_NfixCarbono<-lm.beta(Modelo_Final_Nfix_ratiocell_OCSTHAGLM)
coefficients_NfixNitro<-lm.beta(Modelo_Final_Nfix_ratiocell_NitrogenoGLM)

#AM
coefficients_AMCarbono<-lm.beta(Modelo_final_AM_ratiocellCarbonoGLM)
coefficients_AMNitrogeno<-lm.beta(Modelo_final_AM_ratiocellNitrogenoGLM_AIC)

#EM
coefficients_EMNitrogeno<-lm.beta(Modelo_Final_EM_ratiocell_NitrogenoGLM)
coefficients_EMCarbono<-lm.beta(Modelo_final_EM_ratiocell_CarbonoGLM)


#ErM
coefficients_ErMNitrogeno<-lm.beta(Modelo_Final_ErM_ratiocell_NitrogenoGLM)
coefficients_ErMCarbono<-lm.beta(Modelo_Final_ErM_ratiocell_CarbonoGLM)

#OM
coefficients_OMNitrogeno<-lm.beta(Modelo_Final_OM_ratiocell_NitrogenoGLM)
coefficients_OMCarbono<-lm.beta(Modelo_Final_OM_ratiocell_CarbonoGLM)

#WamNm

coefficients_WanNmNitrogeno<-lm.beta(Modelo_Final_WanNm_ratiocell_NitrogenoGLM)
coefficients_WanNmCarbono<-lm.beta(Modelo_Final_WanNm_ratiocell_CarbonoGLM)


print(coefficients_WanNmNitrogeno)
  
library(MuMIn)
library(Metrics)

Modelos_Nitrogeno<-list(Modelo_final_AM_ratiocellNitrogenoGLM_AIC,Modelo_final_AM_ratiocellNitrogenoZIB_AIC,
                        Modelo_Final_EM_ratiocell_NitrogenoGLM,Modelo_Final_EM_ratiocell_NitrogenoZIB,
                        Modelo_Final_ErM_ratiocell_NitrogenoGLM, Modelo_Final_ErM_ratiocell_NitrogenoZIB,
                        Modelo_Final_Nfix_ratiocell_NitrogenoGLM,Modelo_Final_Nfix_ratiocell_NitrogenoZIB,
                        Modelo_Final_OM_ratiocell_NitrogenoGLM,Modelo_Final_OM_ratiocell_NitrogenoZIB,
                        Modelo_Final_WanNm_ratiocell_NitrogenoGLM, Modelo_Final_WanNm_ratiocell_NitrogenoZIB)

Modelos_Nitrogeno_<-c("Modelo_final_AM_ratiocellNitrogenoLM_AIC",
                      "Modelo_final_AM_ratiocellNitrogenoZIB_AIC",
                      "Modelo_Final_EM_ratiocell_NitrogenoLM","Modelo_Final_EM_ratiocell_NitrogenoZIB",
                      "Modelo_Final_ErM_ratiocell_NitrogenoLM", "Modelo_Final_ErM_ratiocell_NitrogenoZIB",
                      "Modelo_Final_Nfix_ratiocell_NitrogenoLM","Modelo_Final_Nfix_ratiocell_NitrogenoZIB",
                      "Modelo_Final_OM_ratiocell_NitrogenoLM","Modelo_Final_OM_ratiocell_NitrogenoZIB",
                      "Modelo_Final_WanNm_ratiocell_NitrogenoLM","Modelo_Final_WanNm_ratiocell_NitrogenoZIB")

Modelos_OCSTHA<-list(Modelo_final_AM_ratiocellCarbonoGLM,Modelo_final_AM_ratiocellCarbonoZIB,
                     Modelo_final_EM_ratiocell_CarbonoGLM,
                     Modelo_final_EM_ratiocell_CarbonoZIB,
                     Modelo_Final_ErM_ratiocell_CarbonoGLM, Modelo_Final_ErM_ratiocell_CarbonoZIB,
                     Modelo_Final_Nfix_ratiocell_OCSTHAGLM,Modelo_Final_Nfix_ratiocell_OCSTHAZIB,
                     Modelo_Final_OM_ratiocell_CarbonoGLM,Modelo_Final_OM_ratiocell_CarbonoZIB,
                     Modelo_Final_WanNm_ratiocell_CarbonoGLM, Modelo_Final_WanNm_ratiocell_CarbonoZIB)


Modelos_Carbono<-c("Modelo_final_AM_ratiocellCarbonoLM","Modelo_final_AM_ratiocellCarbonoZIB",
                   "Modelo_final_EM_ratiocell_CarbonoLM",
                   "Modelo_final_EM_ratiocell_CarbonoZIB",
                   "Modelo_Final_ErM_ratiocell_CarbonoLM", "Modelo_Final_ErM_ratiocell_CarbonoZIB",
                   "Modelo_Final_Nfix_ratiocell_OCSTHALM","Modelo_Final_Nfix_ratiocell_OCSTHAZIB",
                   "Modelo_Final_OM_ratiocell_CarbonoLM","Modelo_Final_OM_ratiocell_CarbonoZIB",
                   "Modelo_Final_WanNm_ratiocell_CarbonoLM", "Modelo_Final_WanNm_ratiocell_CarbonoZIB")

AIC_Nitrogeno<-lapply(Modelos_Nitrogeno, FUN = AIC)

#Revisar 


Rmse_Nitrogeno<-c(rmse(Data10kmFiltrado_101$AM_ratiocell,predict(Modelo_final_AM_ratiocellNitrogenoGLM_AIC)),
                  rmse(Data10kmFiltrado_101$AM_ratiocell,predict(Modelo_final_AM_ratiocellNitrogenoZIB_AIC)),
                  rmse(Data10kmFiltrado_101$EM_ratiocell,predict(Modelo_Final_EM_ratiocell_NitrogenoGLM)),
                  rmse(Data10kmFiltrado_101$EM_ratiocell,predict(Modelo_Final_EM_ratiocell_NitrogenoZIB)),
                  rmse(Data10kmFiltrado_101$ErM_ratiocell,predict(Modelo_Final_ErM_ratiocell_NitrogenoGLM)),
                  rmse(Data10kmFiltrado_101$ErM_ratiocell,predict(Modelo_Final_ErM_ratiocell_NitrogenoZIB)),
                  rmse(Data10kmFiltrado_101$Nfix_ratiocell,predict(Modelo_Final_Nfix_ratiocell_NitrogenoGLM)),
                  rmse(Data10kmFiltrado_101$Nfix_ratiocell,predict(Modelo_Final_Nfix_ratiocell_NitrogenoZIB)),
                  rmse(Data10kmFiltrado_101$OM_ratiocell,predict(Modelo_Final_OM_ratiocell_NitrogenoGLM)),
                  rmse(Data10kmFiltrado_101$OM_ratiocell,predict(Modelo_Final_OM_ratiocell_NitrogenoZIB)),
                  rmse(Data10kmFiltrado_101$WanNm_ratiocell,predict(Modelo_Final_WanNm_ratiocell_NitrogenoGLM)),
                  rmse(Data10kmFiltrado_101$WanNm_ratiocell,predict(Modelo_Final_WanNm_ratiocell_NitrogenoZIB)))

Rmse_Carbono<-c(rmse(Data10kmFiltrado_101$AM_ratiocell,
                     predict(Modelo_final_AM_ratiocellCarbonoGLM)),
                rmse(Data10kmFiltrado_101$AM_ratiocell,predict(Modelo_final_AM_ratiocellCarbonoZIB)),
                rmse(Data10kmFiltrado_101$EM_ratiocell,predict(Modelo_final_EM_ratiocell_CarbonoGLM)),
                rmse(Data10kmFiltrado_101$EM_ratiocell,predict(Modelo_final_EM_ratiocell_CarbonoZIB)),
                rmse(Data10kmFiltrado_101$ErM_ratiocell,predict(Modelo_Final_ErM_ratiocell_CarbonoGLM)),
                rmse(Data10kmFiltrado_101$ErM_ratiocell,predict(Modelo_Final_ErM_ratiocell_CarbonoZIB)),
                rmse(Data10kmFiltrado_101$Nfix_ratiocell,predict(Modelo_Final_Nfix_ratiocell_OCSTHAGLM)),
                rmse(Data10kmFiltrado_101$Nfix_ratiocell,predict(Modelo_Final_Nfix_ratiocell_OCSTHAZIB)),
                rmse(Data10kmFiltrado_101$OM_ratiocell,predict(Modelo_Final_OM_ratiocell_CarbonoGLM)),
                rmse(Data10kmFiltrado_101$OM_ratiocell,predict(Modelo_Final_OM_ratiocell_CarbonoZIB)),
                rmse(Data10kmFiltrado_101$WanNm_ratiocell,predict(Modelo_Final_WanNm_ratiocell_CarbonoGLM)),
                rmse(Data10kmFiltrado_101$WanNm_ratiocell,predict(Modelo_Final_WanNm_ratiocell_CarbonoZIB)))








resultados_Nitrogeno<-cbind(Modelos_Nitrogeno_,Rmse_Nitrogeno)

resultados_Carbono<-cbind(Modelos_Carbono, Rmse_Carbono)

View(resultados_Nitrogeno)
View(resultados_Carbono)



### Small test (to remove later)


predShannon<-predict(shannonCarbonoGLMFINAL)
elev<-mayores_a_100m$Elev10kmcell

cut_elev<-cut(elev,breaks = seq(100,4600,by=200))
boxplot(predShannon~cut_elev)
predShannon<-predict(shannonNitrogenoGLM)
boxplot(predShannon~cut_elev)



plot(mayores_a_100m$Nitrogeno~mayores_a_100m$Elev10kmcell)
plot(mayores_a_100m$OCSTHA~mayores_a_100m$Elev10kmcell)
plot(mayores_a_100m$shannon10kmEspecies~mayores_a_100m$Elev10kmcell)

plot(predShannon~mayores_a_100m$Elev10kmcell)
