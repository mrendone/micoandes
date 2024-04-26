
##################################
#                                #
#  Null models and rarefactions  #
#                                #
##################################

file_pl_fg <- "../insumos/Idcell_plantas_hongos.csv"
pl_fg <- read.csv(file_pl_fg)

#Here we apply a filter:
#1. we remove the data with no microorganism known association
#2. we remove the cells with less than 10 plant registers

pl_fg<-pl_fg[apply(pl_fg[,c("wam_nm", "nfix", "am", "om", "em", "erm")],1,sum,na.rm=T)>0,]
moreThan10<-table(pl_fg$grid10km_cell)>=10
pl_fg<-pl_fg[pl_fg$grid10km_cell%in%names(moreThan10)[moreThan10],]
fungi <- unique(pl_fg[,c("name_gn", "wam_nm", "nfix", "am", "om", "em", "erm")])
any(duplicated(fungi$name_gn))
fungi <- as.matrix(fungi[,-which(colnames(fungi)=="name_gn")])
fungi[which(is.na(fungi), arr.ind = T)] <- FALSE
list_gn <- as.list(tapply(pl_fg$name_gn,pl_fg$grid10km_cell,function(x)x))
hist(sapply(list_gn, length), nclass=100)

# Assign the row names of the 'fungi' dataframe to the values in the 'name_gn' column
rownames(fungi) <- fungi$name_gn

# Convert the 'fungi' dataframe to a matrix, excluding the 'name_gn' column
fungi <- as.matrix(fungi[,-which(colnames(fungi)=="name_gn")])

# Replace all NA values in the 'fungi' matrix with FALSE
fungi[which(is.na(fungi), arr.ind = T)] <- FALSE

# Create a list 'list_gn' where each element is a vector of fungus names ('name_gn') grouped by 'grid10km_cell'
list_gn <- as.list(tapply(pl_fg$name_gn,pl_fg$grid10km_cell,function(x)x))
hist(sapply(list_gn, length), nclass=100)

##################################################################################################################
#Function for the rarefaction                                                                                    #
#The function for the rarefaction get various arguments:                                                         #
#  • list_tax is a simple list of the genera which are present in each sampling units (in our case the cells of  #
#                                                                                      the grid)                 #
#• nb_tax is the number of individuals that we want to rarefy to                                                 #
#• nb_rep is the number of repetitions of the process                                                            #
#• min is the minimal number of individuals that is required for a sampling unit to enter the process of         #
#rarefaction                                                                                                     #
#• replace is whether the rarefaction is made with reputting the whole individual list at each draw, in          #
#case we want to use some sort of bootstrap process                                                              #
##################################################################################################################


raref <- function(list_tax, nb_tax=10, nb_rep=10000, min=nb_tax*1.5, replace=F){
  list_tax_filtered <- list_tax[sapply(list_tax,length) >= min]
  lapply(1:nb_rep,function(x,lt,s,r)
    lapply(lt,sample,size=s, replace=replace),
    lt=list_tax_filtered,s=nb_tax,r=replace)
}

require(parallel)

numbers_raref<-c(5,10,20,30,40,50,75,100,200)
mclapply(numbers_raref,function(n,lt){
  assign(paste("raref",n,sep="_"),raref(lt,nb_tax=n))
  save(list=paste("raref",n,sep="_"),file=paste0("raref_",n,".RData"))
},lt=list_gn
)



# Combine two lists from the 'matMycoRaref_5' list into a new list 'all_matrices_lists'
all_matrices_lists<-c(list(matMycoRaref_5[[1]]), list(matMycoRaref_5[[2]]))

# Read an Excel file named 'Data10kmFiltrado13_12_20.xlsx' from the current directory and assign it to 'Data10kmFiltrado13_12_20'
library(readxl)
Data10kmFiltrado13_12_20 <- read_excel("../insumos/Data10kmFiltrado13_12_20.xlsx")

# Define a function 'export_results_plots_function' with two parameters 'y' and 'h'
export_results_plots_function <- function(y, h) {
  # Define a function 'MycoShannon' that calculates the Shannon diversity index of a given vector 'x'
  MycoShannon <- function(x) {
    shannon <- as.vector(vegan::diversity(x, index = "shannon"))
    return(shannon)
  }
  
  # Define a function 'unir_shannon' that combines a given vector 'x' with the first, 43rd, and 44th columns of 'Data10kmFiltrado13_12_20'
  unir_shannon <- function(x) {
    Altura_Carbono <- as.matrix(Data10kmFiltrado13_12_20[, c(1, 43, 44)])
    datos <- as.data.frame(cbind(x, Altura_Carbono))
    return(datos)
  }
  
  
  # Initialize several lists to store the results
  results_list <- list()
  results_list_final <- list()
  combined_results_list <- list()
  summary_list <- list()
  plot_list <- list()
  
  # Apply the 'MycoShannon' function to each list of matrices in 'y' and store the results
  for (i in seq_along(y)) {
    result <- lapply(y[[i]], MycoShannon)
    results_list[[i]] <- result
    
    result_final <- lapply(result, unir_shannon)
    results_list_final[[i]] <- result_final
    
    df_shannon <- unlist(result_final)
    combined_Shannon <- bind_rows(unlist(df_shannon), .id = "IDceldas10km.IDcelda10km")
    combined_results_list[[i]] <- combined_Shannon
    
    # Calculate the mean, 5th percentile, and 95th percentile of the Shannon diversity index for each group in 'combined_Shannon'
    summary_Shannon <- combined_Shannon %>%
      group_by(IDceldas10km.IDcelda10km) %>%
      summarize(s_mean = mean(x),
                s_p5 = quantile(x, probs = 0.05),
                s_p95 = quantile(x, probs = 0.95),
                Elev = first(Elev10kmcell))
    summary_list[[i]] <- summary_Shannon
    
    # Create a plot of the mean Shannon diversity index against 'Elev10kmcell'
    plot <- ggplot(summary_Shannon, aes(x = Elev, y = s_mean)) +
      geom_line(color = "blue", size = 1) +
      geom_point(color = "black", size = 3) +
      geom_segment(aes(xend = Elev, yend = s_p5), linetype = "solid", color = "black") +
      geom_segment(aes(xend = Elev, yend = s_p95), linetype = "solid", color = "black") +
      labs(title = "shannon mean vs Elev10kmcell",
           x = "Elev10kmcell",
           y = "s_mean") +
      theme_minimal()
    plot_list[[i]] <- plot
    
    # Save the summary and plot to a .RData file in the directory specified by 'h'
    filename <- paste0(h, "/summary_plot_", i, ".RData")
    data <- list(
      summary = summary_list[[i]],
      plot = plot_list[[i]]
    )
    save(data, file = filename)
  }
}

# Call the 'export_results_plots_function' function with 'all_matrices_lists' and 'save_directory' as arguments
export_results_plots_function(all_matrices_lists, save_directory)

(load("../insumos/matMycoRaref_20.RData"))
# Apply the 'unique' function to the row names of each matrix in 'matMycoRaref_20'
apply(sapply(matMycoRaref_20,rownames),1,unique)

# Check if the 'IDceldas10km.IDcelda10km' column in 'Data10kmFiltrado13_12_20' matches the row names of the first matrix in 'matMycoRaref_20'
Data10kmFiltrado13_12_20$IDceldas10km.IDcelda10km==rownames(matMycoRaref_20[[1]])

# Define a function 'MycoShannon' that calculates the Shannon diversity index of a given matrix 'matMyco'
MycoShannon <- function(matMyco) {
  shannon <- vegan::diversity(matMyco, index = "shannon")
  return(shannon)
}

# Apply the 'MycoShannon' function to each matrix in 'matMycoRaref_20' and calculate the mean and 5th and 95th percentiles of the results
sha_raref<-sapply(matMycoRaref_20,MycoShannon)
mean_sha_raref<-apply(sha_raref,1,mean)
quant_sha_raref<-Reduce(rbind,apply(sha_raref,1,quantile,probs=c(.05,.95),simplify=F))
colnames(quant_sha_raref)<-c("s_p5","s_p95")
rownames(quant_sha_raref)<-rownames(matMycoRaref_20[[1]])
tabShaRaref<-data.frame(mean_sha=mean_sha_raref,
                        quant_sha_raref,
                        Data10kmFiltrado13_12_20[match(names(mean_sha_raref),Data10kmFiltrado13_12_20$IDceldas10km.IDcelda10km),c("Elev10kmcell","OCSTHA")])

# Load each 'matMycoRaref_' file, apply the 'MycoShannon' function to it, calculate the mean and 5th and 95th percentiles of the results, and save the results to a .RData file
for(i in c(5,10,20,30,40,50,75,100,200))
{
  load(paste0("matMycoRaref_",i,".RData"))
  varName <- paste0("matMycoRaref_",i)
  sha_raref <- sapply(get(varName),MycoShannon)
  mean_sha_raref <- apply(sha_raref,1,mean)
  quant_sha_raref<-Reduce(rbind,apply(sha_raref,1,quantile,probs=c(.05,.95),simplify=F))
  colnames(quant_sha_raref)<-c("s_p5","s_p95")
  rownames(quant_sha_raref)<-rownames(get(varName)[[1]])
  tabShaRaref<-data.frame(mean_sha=mean_sha_raref,
                          quant_sha_raref,
                          Data10kmFiltrado13_12_20[match(names(mean_sha_raref),Data10kmFiltrado13_12_20$IDceldas10km.IDcelda10km),c("Elev10kmcell","OCSTHA")])
  resName <- paste0("tabShaEnvirRaref_",i)
  assign(resName,tabShaRaref,envir=.GlobalEnv)
  save(list=resName,file=paste0(resName,".RData"))
  rm(list=varName)
}

}


(load("../insumos/matMycoRaref_100.RData"))

# Define a function 'ratiocell' that calculates the ratio of each column to the row sum of selected columns in a given matrix 'x'
ratiocell<-function(x){
  res<-cbind(
    AM_SP_ratiocell = x[,"am"] / rowSums(x[,c("am","em","erm","om","wam_nm")]),
    EM_SP_ratiocell = x[,"em"] / rowSums(x[,c("am","em","erm","om","wam_nm")]),
    ErM_SP_ratiocell = x[,"erm"] / rowSums(x[,c("am","em","erm","om","wam_nm")]),
    OM_SP_ratiocell = x[,"om"] / rowSums(x[,c("am","em","erm","om","wam_nm")]),
    WaNM_SP_ratiocell = x[,"wam_nm"] / rowSums(x[,c("am","em","erm","om","wam_nm")]),
    Nfix_SP_ratiocell = x[,"nfix"] / rowSums(x[,c("am","em","erm","om","wam_nm","nfix")])
  )
  rownames(res)<-rownames(x)
  return(res)
}

# Apply the 'ratiocell' function to each matrix in 'matMycoRaref_100' and assign the results to 'ratiocell_raref_100'
ratiocell_raref_100 <- lapply((matMycoRaref_100),ratiocell) 

# Calculate the mean of 'ratiocell_raref_100' and assign it to 'matMeanRatioCell_raref100'
matMeanRatioCell_raref100 <- Reduce("+", ratiocell_raref_100) / length(ratiocell_raref_100)

# Save 'ratiocell_raref_100' and 'matMeanRatioCell_raref100' to a .RData file named 'ratiocell_raref_100.RData'
save(list=c("ratiocell_raref_100","matMeanRatioCell_raref100"), file = "ratiocell_raref_100.RData")


# Define a vector 'archivos_a_cargar' containing the names of the RData files to be loaded
archivos_a_cargar <- c(
  "tabShaEnvirRaref_5", "tabShaEnvirRaref_10", "tabShaEnvirRaref_20",
  "tabShaEnvirRaref_30", "tabShaEnvirRaref_40", "tabShaEnvirRaref_50",
  "tabShaEnvirRaref_75", "tabShaEnvirRaref_100", "tabShaEnvirRaref_200"
)

# Define a function 'cargar_y_asignar_contenido' that loads an RData file and assigns its content to the global environment
cargar_y_asignar_contenido <- function(archivo) {
  load(paste0(archivo, ".RData")) # Load the RData file
  assign(archivo, get(archivo), envir = .GlobalEnv) # Assign the content of the RData file to the global environment
}

# Use 'lapply' to apply the 'cargar_y_asignar_contenido' function to each file name in 'archivos_a_cargar'
lapply(archivos_a_cargar, cargar_y_asignar_contenido)

library(ggplot2)
library(dplyr)


#############
#  Graphics #
#############

# Create a list 'plot_list_Elev' to store the plots
plot_list_Elev <- list()

# Create a list 'summary_list' containing the data frames to be plotted
summary_list <- list(
  tabShaEnvirRaref_5, tabShaEnvirRaref_10, tabShaEnvirRaref_20, tabShaEnvirRaref_30, tabShaEnvirRaref_40, tabShaEnvirRaref_50,
  tabShaEnvirRaref_75, tabShaEnvirRaref_100, tabShaEnvirRaref_200
)

# Define a vector 'n_values' containing the corresponding n values
n_values <- c(5, 10, 20, 30, 40, 50, 75, 100, 200)

# Apply the plot code to each element of 'summary_list' in a loop
for (i in seq_along(summary_list)) {
  # Calculate the Pearson correlation coefficient and its significance
  cor_test_result <- cor.test(summary_list[[i]]$Elev10kmcell, summary_list[[i]]$mean_sha)
  
  # Determine whether the result is significant
  is_significant <- cor_test_result$p.value < 0.05
  
  # Create the title with the n value, the correlation coefficient, and its significance
  if (is_significant) {
    significance_label <- "Significant"
  } else {
    significance_label <- "Not Significant"
  }
  
  title <- paste(
    "shannon mean vs Elev10kmcell (n =", n_values[i], ")",
    "\nPearson Correlation =", round(cor_test_result$estimate, 2),
    "\np-value =", format(cor_test_result$p.value, scientific = TRUE, digits = 2),
    "\nResult:", significance_label
  )
  
  # Create the plot
  plot <- ggplot(summary_list[[i]], aes(x = Elev10kmcell, y = mean_sha)) +
    geom_segment(aes(xend = Elev10kmcell, yend = s_p5), linetype = "dotdash", color = "grey") +
    geom_segment(aes(xend = Elev10kmcell, yend = s_p95), linetype = "dotdash", color = "grey") +
    geom_point(color = "black", size = 0.75) +
    labs(
      title = title,
      x = "Elev10kmcell",
      y = "s_mean"
    ) +
    theme_minimal()
  
  # Add the plot to 'plot_list_Elev'
  plot_list_Elev[[i]] <- plot
}

#plot list 

for (i in seq_along(plot_list_Elev)) {
  print(plot_list_Elev[[i]])
}


