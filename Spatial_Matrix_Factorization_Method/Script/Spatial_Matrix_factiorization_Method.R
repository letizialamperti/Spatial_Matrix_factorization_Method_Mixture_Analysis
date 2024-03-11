# Set the working directory
setwd("path_to_your_directory")

# Load necessary libraries
library(maps)
library(tess3r)
library(mapplots)


# Read in the data:your community matrix could be abundance or occurrence matrix of species or MOTU..
# Please replace "path_to_your_directory" with the appropriate path to your data files

# Community matrix
motu_global<- read.csv("")

# Metadata
metadata_motu_global <- read.csv("")
taxo_identification <- motu_global[ ,1:6]

# Define the number of PCR replicas
num_PCR <- 12

# Transpose the data
Data <- t(motu_global[, 7:ncol(motu_global)])
dim(Data)
colnames(Data) <- motu_global$sequence

# Get the dimensions
n_samples <- dim(Data)[1]
n_MOTU <- dim(Data)[2]

# Replace NA values with 0
Data[is.na(Data)] <- 0


# Order Data by spygen code
Data <- Data[order(row.names(Data)), ]
metadata_motu_global <- metadata_motu_global[order(metadata_motu_global$spy_code), ]

# Divide Data by number of PCR
Data <- Data / num_PCR

# Check uniqueness
unique(metadata_motu_global$spy_code == row.names(Data))
unique(colnames(Data) == taxo_identification$sequence)


###################################################################################################################

# Construct Data_double for tess3 function
#See Appendix 1 from paper for more information

Data_double <- matrix(0, n_samples, 2 * n_MOTU)

for (i in 1:ncol(Data)) {
  j = 2 * i
  Data_double[, j - 1] <- Data[, i]
  Data_double[, j] <- 1 - Data[, i]
}

###################################################################################################################

# tess3

###################################################################################################################


# Run tess3 function to perform spatial clustering

# The parameter 'K' defines the range of sources to test
# 'lambda' controls clustering based on spatial information

coord <- cbind(long = as.numeric(metadata_motu_global$longitud), lat = as.numeric(metadata_motu_global$latitud))

# tess3 funtion 
tess3.obj <- tess3(X = NULL, XProba = Data_double, coord = coord, K = (3:9), ploidy = 1, lambda = 10, rep = 10, W = NULL,  
                   method = "projected.ls", max.iteration = 200, tolerance = 1e-05, 
                   openMP.core.num = 1, Q.init = NULL, mask = 0, algo.copy = TRUE,
                   keep = "best", verbose = FALSE)

# Define a function to calculate mean cross-entropy score
fun = function(x) mean(x$crossentropy)

# Calculate the mean cross-entropy score for each K
mean_cross_score <- sapply(tess3.obj, fun)

# Plot the mean cross-entropy score against the number of sources (K)
plot(mean_cross_score[1:9], pch = 19, col = "blue", xlab = "Number of K sources", 
     ylab = "Cross-validation score") # to choose K

# retrieve tess3 S matrix through the qmatrix function of tess3r package for K = 6 sources 

s.matrix <- qmatrix(tess3.obj, K = 6)


###################################################################################################################

# STRUCTURE-like barplot for the S-matrix 

#colors
col_6 = c("yellow2","purple","sienna1",'#FF3E96',"royalblue3","palegreen3" )

m_s_matrix <- as.matrix(t(s.matrix[order(metadata_motu_global$region),]))
barplot(m_s_matrix, border = NA, space = 0, 
        xlab = "", ylab = "samples proportions", xaxt='n',
        col = col_6) -> bp
bp
axis(1, at = 1:nrow(s.matrix), labels = metadata_motu_global$region[order(metadata_motu_global$region)], las = 2,  cex.axis = .4)  #axes


###################################################################################################################

#Pie-plot World

# Plot the coordinates on a map
library(mapplots)
plot(coord[,], xlab = "Longitude", ylab = "Latitude", type = "n")
map(add = T, col = "grey90", fill = TRUE)
points(coord[,], cex = 1)

# Add pie charts representing source proportions
for (i in 1:nrow(coord)){
  add.pie(z = s.matrix[i,], x = coord[i,1], y = coord[i,2], labels ='', radius = 4,
          col = col_6)}


###################################################################################################################

#P-values

###################################################################################################################

#extracting p-values through the pvalue funciton of the tess3r package
p.values <- pvalue(tess3.obj, K = 6)
hist(p.values, col = "lightblue") 

p_values <- data.frame (id = c(1:n_MOTU),pvalues = -log10(as.numeric(p.values)), taxo_identification = taxo_identification_motu_global$scientific_name_ncbi_corrected)

# Benjamini-Hochberg algorithm
L = length(p.values)
fdr.level = 1e-30
w = which(sort(p.values) < fdr.level * (1:L)/L)
candidates = order(p.values)[w]
length(candidates)

# manhattan plot 
plot(p.values, main = "Manhattan plot", 
     xlab = "MOTUs", 
     ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates, -log10(p.values)[candidates], 
       pch = 19, cex = .2, col = "blue")

###################################################################################################################

# G matrix

###################################################################################################################

# Set the number of sources
K = 6

# Extract M matrix from tess3 object
m.matrix <- tess3.obj[[K]][["tess3.run"]][[1]][["G"]]

# Get taxa indices from taxo_identification
taxa_indecies_ncbi <- taxo_identification$scientific_name_ncbi_corrected


# Create a clean M matrix
clean_m_matrix <- matrix (0, n_MOTU, K)
rownames(clean_m_matrix) <- taxo_identification$species_name_corrected
colnames(clean_m_matrix) <- col_6

# Fill in the clean M matrix with values from the M matrix
# We select the frequencies of the MOTUs that correspond to the p detection probabilities coming from D_double

for (i in 1:n_MOTU) {
  clean_m_matrix[i,] = m.matrix[(2*i -1),]
}

# Find the maximum value for each row in the clean G matrix

maxiu <- apply(clean_m_matrix, 1, which.max)
clean_m_matrix_check <- NULL
clean_m_matrix_check <- cbind(clean_m_matrix, maxiu, 1:n_MOTU)

# Sort the maxiu values

maxiu<- sort(maxiu)
family_id <- taxo_identification$family_name_corrected



# Convert clean_m_matrix to data frame

clean_m_matrix <- as.data.frame(clean_m_matrix)

# Barplot for MOTU proportions
barplot(t(clean_m_matrix[,]), border = NA, space = 0, 
        xlab = "", ylab = "MOTU proportions", xaxt='n',
        col = col_6)
axis(1, at = 1:nrow(clean_m_matrix), labels = 1:n_MOTU, las = 2,  cex.axis = .2)  #axes

# Create clean_m_matrix_1 to visualize the most aboundant MOTUs

clean_m_matrix_1 <- matrix(0, n_MOTU, K+2)

# Loop through rows to filter clean_m_matrix

for (i in 1:n_MOTU) {
  if (sum(clean_m_matrix[i,]) > 0.3 ){
    clean_m_matrix_1[i,]  <- clean_m_matrix[i,]  }
}

# Set rownames for clean_m_matrix_1
rownames(clean_m_matrix_1) <- rownames(clean_m_matrix)


# Barplot for MOTU proportions
barplot(t(clean_m_matrix_1[,1:K]), border = NA, space = 0, 
        xlab = "", ylab = "MOTU proportions", xaxt='n',
        col = col_6)

axis(1, at = 1:nrow(clean_m_matrix_1), outer = FALSE, labels = 1:nrow(clean_m_matrix_1), las = 2,  cex.axis = 1)  #axes


###############################################################################################################################
###############################################################################################################################
