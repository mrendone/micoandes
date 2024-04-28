# micoandes

Here we share the code used for Rendon et al. “Diversity of mycorrhizal types along
altitudinal gradients in the tropical andes”.

# The databases

The database used for this study compiles 266,543 plant records from herbarium collections
across the Colombian Andes and includes 3,094 plant genera. This database includes all the
records included in previously published by Bottin et al. (2020) including herbarium data from
the global aggregator GBIF (http://doi.org/10.15468/dl.xqndaq), the database “Tropicos” of
Missouri Botanical Garden (https://www.tropicos.org), the database from the National
Colombian Herbarium of the Instituto de Ciencias Naturales (COL, Raz and Agudelo 2021).
Details about the data and data standardization protocols are described in Bottin et al. (2020).
We assigned type of mycorrhizal association to plant genera based on the FungalRoot database
1.0 (Soudzilovskaia et al., 2020) and Wang and Qiu, 2006; Steidinger et al., 2019; Vetrovsky et
al., 2020. From all the types of associations, we chose to work with the most abundant types of
mycorrhizal symbiosis, following Brundrett and Tedersoo (2018).
Soil variables were downloaded in raster format from the International Union of Soil Sciences
website (ISRIC, https://www.isric.org/explore/soilgrids/faq-soilgrids ) with a resolution 250 m at
four depths: 0, 0.05, 0.15 and 0.3 m. The soil variable values for all depths were averaged for
each cell using Orfeo Toolbox (Grizonnet et al., 2017). Climatic variables were downloaded in
raster format with a resolution of 30 arc seconds from the Worldclim website
( https://www.worldclim.org/data/worldclim21.html ). Additionally, a digital elevation model
raster was downloaded from the same website to extract the elevation data from the logs with the
same resolution.

# Code organization

## Model.R

The <Models.R> file includes all the code to repeat the linear models as in the article; The script includes linear models that evaluate which are the best environmental and soil variables to explain the diversity of mycorrhizas in the Colombian Andes. Additionally, identify which variables best explain the proportion of mycorrhizae within the same region. The script also incorporates the use of root mean square error (RMSE) to select the most suitable models.

## Null models and rarefactions  

To evaluate the spatial distribution patterns of the mycorrhizal association types, we generated a
spatial grid data frame with a resolution of 10 km 2 for the entire Colombian Andes using the sp
package v1.4 (Bivand et al., 2013). To extract soil, climatic, and elevation data, a resampling was
performed to adjust the resolution to that of the other grid types using the raster package v3.0
(Hijmans, 2019). Once the resampling for the variables was carried out, the values were
extracted per cell. Geographic coordinates of the cell (longitude and latitude) were extracted
using the coordinates function of the raster package. We generated the mean of 10,000 random
estimates of the number of mycorrhizal types per cell at each level of rarefaction.
This index used the relative abundance of each type of mycorrhizal association calculated as the
number of records for each mycorrhizal association type divided by the number of plants
occurrence records in the cell.

GLMs with a Gaussian distribution and identity link were used to explore the influence of our
spatially explicit predictor variables on several response variables: 1) the proportions of the
different mycorrhizal association types per grid cell; 2) the diversity index of the different types
of mycorrhizal associations per grid cell.

The <Nullmodels_rarefactions_and_elevationgraphics.R> provides the code necessary to perform the data rarification process and the null model used in the article. First, it filters the data by removing observations with no known association with microorganisms and cells with fewer than 10 plant records. It then generates a matrix of unique fungal data and calculates rarefaction for different numbers of individuals and replicates, saving the results to files. Next, calculate the Shannon diversity index for each rarefied data set, plotting the relationship between this index and cell elevation. Additionally, it calculates the Pearson correlation between these two values, determining its statistical significance. Finally, create a list of graphs showing the relationship between the Shannon index and elevation for different rarefaction values.

## Figures and supplementary tables 


