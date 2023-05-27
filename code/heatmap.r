###### R code for Complex Heatmap used in Responder-Nonresponder Analysis ######

library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(methods)

###### Reading input files ####
df1 = fread("atleast_2.csv", header = TRUE)
df2 = fread("BL1.csv", header = TRUE)
df3 = fread("BL2.csv", header = TRUE)
df4 = fread("LAR.csv", header = TRUE)
df5 = fread("Mes.csv", header = TRUE)

mat1 <- tibble::column_to_rownames(df1,"Subtype") %>% as.matrix
mat2 <- tibble::column_to_rownames(df2,"Subtype") %>% as.matrix
mat3 <- tibble::column_to_rownames(df3,"Subtype") %>% as.matrix
mat4 <- tibble::column_to_rownames(df4,"Subtype") %>% as.matrix
mat5 <- tibble::column_to_rownames(df5,"Subtype") %>% as.matrix

#### Remove "cluster_columns=FALSE" if want clustering
h1 <- Heatmap(mat1,  name = "Common", column_title = "Common", width = 2, height = 3, column_names_gp = grid::gpar(fontsize = 5), row_names_gp = grid::gpar(fontsize = 8), cluster_columns=FALSE)
h2 <- Heatmap(mat2,  name = "BL1", column_title = "BL1 Unique", width = 5, height = 3, column_names_gp = grid::gpar(fontsize = 5), row_names_gp = grid::gpar(fontsize = 8), cluster_columns=FALSE)
h3 <- Heatmap(mat3,  name = "BL2", column_title = "BL2 Unique", width = 5, height = 3, column_names_gp = grid::gpar(fontsize = 5), row_names_gp = grid::gpar(fontsize = 8), cluster_columns=FALSE)
h4 <- Heatmap(mat4,  name = "LAR", column_title = "LAR Unique", width = 5, height = 3, column_names_gp = grid::gpar(fontsize = 5), row_names_gp = grid::gpar(fontsize = 8), cluster_columns=FALSE)
h5 <- Heatmap(mat5,  name = "M", column_title = "M Unique", width = 5, height = 3, column_names_gp = grid::gpar(fontsize = 5), row_names_gp = grid::gpar(fontsize = 8), cluster_columns=FALSE)

#### Storing multiple heatmaps in one

ht_list =  h1 + h2 + h3 + h4 + h5

#### Drawing the Combined heatmap

draw(ht_list, column_title = "Frequency Comparison", column_title_gp = gpar(fontsize = 15))

