library(dplyr)

#ANLISYS with AGNES DATA (Duhamet et al. 2023)
load("Data/Data_Agnes.RData")

dim(Mat_Pa_actinopterygii)
dim(Mat_Pa_Chondrichthyes)

colnames(Mat_Pa_Chondrichthyes)[1] <- "Longitude"
colnames(Mat_Pa_Chondrichthyes)[2] <- "Latitude"


prova <- inner_join(Mat_Pa_actinopterygii, Mat_Pa_Chondrichthyes, by = c("Longitude" , "Latitude"))


metadata_motu_global <- read.csv("Data/4A_metadata_motu_global_Leti.csv")
MOTU_F_matrix_Species <- read.csv("~/Desktop/K ancestral/4A_MOTU_F_matrix_Species.csv", row.names=1)
MOTU_F_matrix_Species <- subset(MOTU_F_matrix_Species[,], !is.na(MOTU_F_matrix_Species$species))

length(unique(MOTU_F_matrix_Species$species)) #number of species = 791

# Select columns by name using dplyr
selected_columns_prova <- prova %>% dplyr::select(intersect(MOTU_F_matrix_Species$species, colnames((prova))))
species_inter <- colnames(selected_columns_prova)

MOTU_F_matrix_Species <- MOTU_F_matrix_Species[MOTU_F_matrix_Species$species %in% species_inter,]

selected_columns_prova <- selected_columns_prova [,order(species_inter)]
selected_columns_prova <- data.frame(Longitude = prova$Longitude, Latitude = prova$Latitude, selected_columns_prova)



longitude_patch_per_region <- metadata_motu_global %>%
  group_by(region) %>%
  summarize(min_longitude = min(longitud), max_longitude = max(longitud))


latitude_patch_per_region <- metadata_motu_global %>%
  group_by(region) %>%
  summarize(min_latitude = min(latitud), max_latitude = max(latitud))


patch_per_region <- data.frame ( longitude_patch_per_region, min_latitude =  latitude_patch_per_region$min_latitude, max_latitude = latitude_patch_per_region$max_latitude)
coord <- prova[,1:2]

#plotting patch
library(mapplots)
library(maps)
plot(coord[,], xlab = "Longitude", ylab = "Latitude", type = "n")
map(add = T, col = "grey90", fill = TRUE)
col_region = c("royalblue3","palegreen3" , "palegreen3" ,"royalblue3","purple","sienna1","sienna1","yellow2",'#FF3E96',"purple" )
rect(patch_per_region$min_longitude, patch_per_region$min_latitude,patch_per_region$max_longitude,patch_per_region$max_latitude, col = col_region, border = "black")

selected_rows_and_cols_patch <- NULL

#selecting patches all together

longitude_spezzata <- c(-146:-146, -82:-54, -6:-4, 3:9, 11:31, 40:55, 130:135, 163:167)
latitude_spezzata <- c(-65:-58, -23:-11, -5:0, 3:5, 11:17, 39:44, 47:49, 73:82 )


for (i in 1:dim(selected_columns_prova)[1]) {
  if (selected_columns_prova$Longitude[i] %in% longitude_spezzata && selected_columns_prova$Latitude[i] %in% latitude_spezzata) {
    selected_rows_and_cols_patch <- rbind(selected_rows_and_cols_patch, selected_columns_prova[i,])
  }
}

#grouping together MOTUs of same species
MOTU_F_matrix_Species <- MOTU_F_matrix_Species %>%
  group_by(species) %>%
  summarize(across(everything(), mean))

MOTU_F_matrix_Species <- MOTU_F_matrix_Species %>% arrange(species)


#defining PATCHES of different six regions
#MEDITERRANEAN SEA
long_med <- c(-5:-4, 3:9)
lat_med  <- c(39:44, 47:49)
Med_patch <- data.frame()

for (i in 1:dim(selected_rows_and_cols_patch)[1]) {
  if (selected_rows_and_cols_patch$Longitude[i] %in% long_med && selected_rows_and_cols_patch$Latitude[i] %in% lat_med) {
    Med_patch <- rbind(Med_patch, selected_rows_and_cols_patch[i,])
  }
}

#POLES
long_poles <- c(-64:-54, 11:31)
lat_poles  <- c( -65:-58, 73:82)
Poles_patch <- data.frame()

for (i in 1:dim(selected_rows_and_cols_patch)[1]) {
  if (selected_rows_and_cols_patch$Longitude[i] %in% long_poles && selected_rows_and_cols_patch$Latitude[i] %in% lat_poles) {
    Poles_patch <- rbind(Poles_patch, selected_rows_and_cols_patch[i,])
  }
}


#TROPIC ATLANTIC
long_atl <- c(-82:-61, 11:31)
lat_atl  <- c( 3:5, 11:17)
Altlantic_patch <- data.frame()

for (i in 1:dim(selected_rows_and_cols_patch)[1]) {
  if (selected_rows_and_cols_patch$Longitude[i] %in% long_atl && selected_rows_and_cols_patch$Latitude[i] %in% lat_atl) {
    Altlantic_patch <- rbind(Altlantic_patch, selected_rows_and_cols_patch[i,])
  }
}

#WESTERN CORAL TRIANGLE
long_coral_tr <- c(130:135)
lat_coral_tr  <- c(-5:0)
Coral_triangle_patch <- data.frame()

for (i in 1:dim(selected_rows_and_cols_patch)[1]) {
  if (selected_rows_and_cols_patch$Longitude[i] %in% long_coral_tr && selected_rows_and_cols_patch$Latitude[i] %in% lat_coral_tr) {
    Coral_triangle_patch <- rbind(Coral_triangle_patch, selected_rows_and_cols_patch[i,])
  }
}

#NEW CALEDONIA
long_new_cal<- c(163:167)
lat_new_cal  <- c( -23:-19)
new_cal_patch <- data.frame()

for (i in 1:dim(selected_rows_and_cols_patch)[1]) {
  if (selected_rows_and_cols_patch$Longitude[i] %in% long_new_cal  && selected_rows_and_cols_patch$Latitude[i] %in% lat_new_cal ) {
    new_cal_patch <- rbind(new_cal_patch, selected_rows_and_cols_patch[i,])
  }
}


#INDIAN OCEAN
long_indian <- c(-146:-145,40:55)
lat_indian  <- c( -23:-11 )
indian_patch <- data.frame()

for (i in 1:dim(selected_rows_and_cols_patch)[1]) {
  if (selected_rows_and_cols_patch$Longitude[i] %in% long_indian && selected_rows_and_cols_patch$Latitude[i] %in% lat_indian) {
    indian_patch <- rbind(indian_patch, selected_rows_and_cols_patch[i,])
  }
}



final_patch_per_region <- data.frame(species = MOTU_F_matrix_Species$species, new_cal =c(colSums(new_cal_patch[,3:732])),
                                     indian = c(colSums(indian_patch[,3:732])), atlantic = c(colSums(Altlantic_patch[,3:732])),
                                     western_coral_triangle = c(colSums(Coral_triangle_patch[,3:732])), poles = c(colSums(Poles_patch[,3:732])),
                                     med = c(colSums(Med_patch[,3:732])) )


col_6 = c("yellow2","purple","sienna1",'#FF3E96',"royalblue3","palegreen3" )
barplot(t(final_patch_per_region[,2:7]), border = NA, space = 0, 
        xlab = "", ylab = "MOTU proportions", xaxt='n',
        col = col_6)


mean(as.matrix(MOTU_F_matrix_Species[,2:7]))
median(as.matrix(MOTU_F_matrix_Species[,2:7]))
quantile(as.matrix(MOTU_F_matrix_Species[,2:7]))


MOTU_F_matrix_Species_pres_abs <- binary_matrix <- as.numeric(MOTU_F_matrix_Species[,2:7] > 9.999895e-06 )
final_patch_per_region_binary <- as.numeric(final_patch_per_region[,2:7] > 0 )


cor(rowSums(as.matrix(MOTU_F_matrix_Species_pres_abs )), rowSums(as.matrix(final_patch_per_region_binary)))


