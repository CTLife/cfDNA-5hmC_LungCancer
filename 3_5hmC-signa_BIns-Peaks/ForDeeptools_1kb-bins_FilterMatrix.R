



###########################################################################################################################################################################
rawMatrix_1 <- read.table("B.1-RawCounts.1000bp-Bin.txt", header=TRUE,   sep="\t", comment.char = "" )  
dim(rawMatrix_1 )
rawMatrix_1[1:10,1:10]

myBED = rawMatrix_1[, c(1,2,3)]
dim(myBED)
myBED[1:10,]

rawMatrix_2 = rawMatrix_1[, -c(1,2,3)]
dim(rawMatrix_2 )
rawMatrix_2[1:10,1:10]

library(stringr)
rownames1 = apply( myBED, 1, paste, sep="" , collapse = "..." )
length(rownames1)
head( rownames1 )
rownames1[1:100]

myBED[,4] = rownames1
rownames(rawMatrix_2) = rownames1
dim(myBED)
myBED[1:10,]
dim(rawMatrix_2 )
rawMatrix_2[1:10,1:10]


myLen = myBED[,3] - myBED[,2]
myBED[,5] = myLen
dim(myBED)
myBED[1:10,]

myBED[,6] = rep(".",   length(myLen) )
dim(myBED)
myBED[1:10,]




AllResults_g <- "ForDeeptools_FilterMatrix"
if( ! file.exists(AllResults_g) ) { dir.create(path=AllResults_g, recursive = TRUE) }

write.table( myBED ,   file = paste(AllResults_g, "All.Autosome.peaks.bed", sep="/"), 
             append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".",  row.names = F,  col.names = F )





##########################
info_1 <- read.table("Classified_291_samples.txt", header=TRUE,   sep="\t", comment.char = "" )  
dim(info_1)
info_1[1:10,1:10]
head(info_1)

colnames(info_1)
names_1 = info_1$Sample_name
names_1_order = order(names_1)
names_1[names_1_order]

info_2 = info_1[names_1_order,]
dim(info_2)
info_2[1:10,1:10]
head(info_2)

info_3 = t(info_2)
dim(info_3)
info_3[,1:5]
head(info_3)

##########################





colnames_1 = colnames( rawMatrix_2 )
colnames_1
colnames_1_order = order(colnames_1)
colnames_1_order
colnames_1[colnames_1_order]

rawMatrix_3 = rawMatrix_2[ , colnames_1_order]
rawMatrix_3[1:10,1:10]
dim(rawMatrix_3)

dim(rawMatrix_3)
dim(info_3)
rawMatrix_3[1:5,1:5]
info_3[1:10,1:5]
mybool1 = (info_3[6,] == colnames(rawMatrix_3))
length(mybool1[mybool1])

colnames(info_3) = colnames(rawMatrix_3)
rawMatrix_4 = rbind( info_3, rawMatrix_3 )
dim(rawMatrix_4)
rawMatrix_4[1:20,1:10]

colnames( rawMatrix_4 )
rownames( rawMatrix_4 )[1:20]
myBool = ( colnames( rawMatrix_4 ) == rawMatrix_4[6, ] )
length( myBool[myBool]  )
length( myBool[!myBool]  )



write.table( rawMatrix_4 ,   file = paste(AllResults_g, "All.info.txt", sep="/"), 
             append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".",  row.names = T,  col.names = T )






######################################################################################################
## Number of zeros of each row
numOfZero <- function(x) {
  return(length(which(x == 0)))
}
numOfZero1 = apply(rawMatrix_4,  1, numOfZero )
length(numOfZero1)
length(numOfZero1[numOfZero1>290])
length(numOfZero1[numOfZero1>250]) ##used
length(numOfZero1[numOfZero1>200])
length(numOfZero1[numOfZero1>150])
length(numOfZero1[numOfZero1>100])



mybool = (  numOfZero1<250 )
length( mybool )
length( mybool[mybool] )


rawMatrix_5 = rawMatrix_4[mybool,]
dim(rawMatrix_5)
rawMatrix_5[1:20,1:10]

mybool2 = mybool[-c(1:16)]
myBED2 = myBED[mybool2, ]    
dim(myBED2)
myBED2[1:10,]
rawMatrix_5[15:20,1:3]




write.table( myBED2 ,   file = paste(AllResults_g, "Kept.Autosome.peaks.bed", sep="/"), 
             append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".",  row.names = F,  col.names = F )



write.table( rawMatrix_5 ,   file = paste(AllResults_g, "Kept.info.txt", sep="/"), 
             append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".",  row.names = T,  col.names = T )












