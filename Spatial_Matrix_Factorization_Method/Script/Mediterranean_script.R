# Set the working directory
setwd("~/Desktop/Spatial_Matrix_Factorization_Method")

# Load necessary libraries
library(maps)
library(tess3r)
library(mapplots)

# Read data
metadata_motu_global <- read_csv("Data/metadata_motu_global_Leti.csv")
motu_global <- read.csv("Data/motu_global.csv")
meta_data_med <- read_csv("Data/protection_med.csv")

# Constants
num_PCR <- 12

# Preprocess metadata
metadata_motu_global <- metadata_motu_global[order(metadata_motu_global$spy_code),]

# Prepare Data
Data <- motu_global[,7:ncol(motu_global)]
rownames(Data) <- motu_global$sequence
Data[is.na(Data)] <- 0
Data <- Data[, metadata_motu_global$region == "Mediterranean_Sea"]
Data <- subset(Data, rowSums(Data) > 0) 
Data <- t(Data)

# Get dimensions
n_samples <- dim(Data)[1]
n_MOTU <- dim(Data)[2]

# Extract taxonomic information
taxo_identification <- motu_global[motu_global$sequence %in% colnames(Data),1:6]

# Filter metadata
metadata_motu_global <- metadata_motu_global[metadata_motu_global$region == "Mediterranean_Sea",]

# Order Data by spygen code
Data <- Data[order(row.names(Data)),]

# Divide Data per number of PCR
Data <- Data / num_PCR

# Construct Data_double for tess3 function
Data_double <- matrix(0, n_samples, 2 * n_MOTU)   
for (i in 1:ncol(Data)) {
  j <- 2 * i
  Data_double[,j-1] <- Data[,i] 
  Data_double[,j] <- 1 - Data[,i]
}

# Run tess3 function to perform spatial clustering
coord <- cbind(long = as.numeric(metadata_motu_global$longitud), lat = as.numeric(metadata_motu_global$latitud))
tess3.obj <- tess3(X = NULL, XProba = Data_double, coord = coord, K = (2:10), ploidy = 1, lambda =1, rep = 10, W = NULL,
                   method = "projected.ls", max.iteration = 100, tolerance = 1e-06, 
                   openMP.core.num = 1, Q.init = NULL, mask = 0, algo.copy = TRUE,
                   keep = "best", verbose = FALSE)

# Plot mean cross-validation score
fun <- function(x) mean(x$crossentropy)
mean_cross_score <- sapply(tess3.obj, fun)
plot(mean_cross_score, pch = 19, col = "blue", xlab="Number of K clusters", ylab="Cross-validation score")

# Retrieve tess3 Q matrix for K = 5 clusters 
s.matrix <- qmatrix(tess3.obj, K = 5)

# STRUCTURE-like barplot for the Q-matrix
col_5 <- c("palegreen3", "sienna1", '#FF3E96', "yellow", "royalblue3")


m_s_matrix <- as.matrix(t(s.matrix[order(metadata_motu_global$latitud),]))
barplot(m_s_matrix, border = NA, space = 0, xlab = "", ylab = "proportions of samples", xaxt='n', col = col_5)

# Pie-plot World
plot(coord[,], xlab = "Longitude", ylab = "Latitude", type = "n")
map(add = T, col = "grey90", fill = TRUE)
for (i in 1:nrow(coord)){
  if (meta_data_med$protection_med[i] == "reserve"){
    add.pie(z = s.matrix[i,], x = coord[i,1] , y = coord[i,2] - 0.4, labels = "", radius = 0.1, col = col_5)}}
for (i in 1:nrow(coord)){
  if (meta_data_med$protection_med[i] == "outside"){
    add.pie(z = s.matrix[i,], x = coord[i,1], y = coord[i,2], labels ='', radius = 0.07, col = col_5)}}
points(coord[meta_data_med$protection_med == "reserve",], cex = .5, pch= 8, col = "red")

# Model fitting
protection_med_binary <- ifelse(meta_data_med$protection_med == "reserve", 1, 0)
s.matrix_comp <- ilr(s.matrix)
my_glm <- glm(protection_med_binary ~ s.matrix_comp, family = binomial)
Data_presence <- as.matrix((Data > 0) + 0)
motu_richness <- rowSums(Data_presence)
lm <- lm(motu_richness ~ s.matrix_comp)
pseudo_r_squared <- pR2(my_glm)

# Principle Component Analysis
res.pca <- prcomp(Data, scale = TRUE)
eigs <- res.pca$sdev^2
pca_variance <- eigs[1] / sum(eigs)
