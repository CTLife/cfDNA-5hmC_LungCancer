##################################################################################################################
## Suffixes of all self-defined global variables must be "_g".
## Example:  
## Rscript  clutering_heatmap.R     3UTR    


library(factoextra) # clustering algorithms & visualization
library(flexclust)
suppressPackageStartupMessages( library(ggplot2)  )
suppressPackageStartupMessages( library(randomcoloR) )
suppressPackageStartupMessages( library(RColorBrewer) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(dendextend) )
suppressPackageStartupMessages( library(factoextra) ) 
suppressPackageStartupMessages( library(MASS) )
suppressPackageStartupMessages( library(ggplot2)  )
suppressPackageStartupMessages( library(scales) )
suppressPackageStartupMessages( library(corrplot) )
suppressPackageStartupMessages( library(gplots) )
suppressPackageStartupMessages( library(Hmisc) )


args_g <- commandArgs(TRUE)
print("##########################")
print("args: ")
print(args_g[1])   
print("##########################")

input_matrix_g= args_g[1];     ## Input matrix file

# input_matrix_g = "3UTR"

outDir_g = paste(input_matrix_g , ".Results",  sep="")
if( ! file.exists(outDir_g)          ) { dir.create(outDir_g, recursive = TRUE) }
##################################################################################################################





###################
DF1 <- read.table( paste(input_matrix_g,   ".txt",  sep="") , header=F, sep="\t", quote = "", comment.char = "") 
dim(DF1)
DF1[1:5,]

matrix   = DF1[, -1]
COL_NAME = DF1[, 1 ]
matrix[1:5,]
COL_NAME[1:5]
max(matrix)
min(matrix)


for(i in c(1:length(COL_NAME)) ){
  COL_NAME[i] = paste(COL_NAME[i], i, sep="...")
}
COL_NAME[1:5]



reset_outliers2 <- function(x, na.rm = TRUE ) {
  qnt <- quantile(x, probs=c(0.01, 0.99) , type=1,  na.rm = na.rm )  
  y <- x
  y[x < qnt[1] ] <- qnt[1]
  y[x > qnt[2] ] <- qnt[2]    
  y
}

myScaleMatrix2 <- function( matrix_temp8, upper_temp8 = 1, lower_temp8 = -1 ) {
  rawMatrix_2 = matrix_temp8  ## reset_outliers2(matrix_temp8)  
  rawMatrix_2 = lower_temp8 + (upper_temp8 - lower_temp8) * ( rawMatrix_2 - min(rawMatrix_2) )/( max(rawMatrix_2)- min(rawMatrix_2) )
  return(rawMatrix_2)
}


  

################################################################################
matrixA = reset_outliers2(matrix) 
max(matrixA)
min(matrixA)

matrixA[1:5, ]
matrix2 <- matrix(as.numeric(unlist(matrixA)), ncol = ncol(matrixA))   # Convert to numeric matrix 
dim( matrix2 )
matrix2 = log2(matrix2+1)
rownames(matrix2) = COL_NAME
## matrix2[1:5, ]


library("pheatmap")
pdf( file = paste(outDir_g, "1.Hierarchical.pdf", sep="/"),  width=3, height=10  )
     presutls = pheatmap(matrix2, cutree_rows = 29 , cluster_cols = F, clustering_distance_rows = "euclidean", show_rownames = F, show_colnames = F  )
dev.off() 

names(presutls)
index1 = presutls$tree_row$order

 
matrix3 = matrix2[index1 , ]
pdf( file = paste(outDir_g, "2.heatmap.pdf", sep="/"),  width=3, height=10  )
   presutls = pheatmap(matrix3,   cluster_rows = F,  cluster_cols = F , show_rownames = F, show_colnames = F  )
dev.off() 


write.table(matrix3 ,  file = paste(outDir_g,   "2.matrix.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".",  row.names = TRUE,  col.names = TRUE )



 
DF2 <- read.table( paste(input_matrix_g, ".SortedRegions",  sep="") , header=F, sep="\t", quote = "", comment.char = "") 
#dim(DF2)
#DF2[1:5,]


DF3 = DF2[index1 , ]
matrix4 = cbind(DF3[,1:6] , matrix3)
#matrix4[1:5,]


write.table(matrix4 ,  file = paste(outDir_g,   "3.matrix.bed",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".",  row.names = F,  col.names = F )

