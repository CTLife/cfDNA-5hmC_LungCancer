





library(ROTS)
library(GenomicRanges)
library(rtracklayer)
library(Rsubread)
library(Rsamtools)
library(ade4)
library(made4)
library(edgeR)
library(EDASeq)
library(RUVSeq)
library(RColorBrewer)
library(ggplot2)




outDir_g      <- "DPs"
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g,   recursive = TRUE) }



##############################################################################################################################################################################################
DF_A <- read.table("../3c.normalizedCounts.withHeader.txt", header=TRUE,   sep="\t", comment.char = "" )  
dim( DF_A )
colnames(DF_A)
## DF_A[1:23,1:5]


## Participant_ID               P-00054784          P-00054071          P-00065404          P-00054495           P-00054314
## Sample_ID                    3481647801          S041700106          S021705891          3481638627           S011702512
## Batch                            Batch1              Batch1              Batch1              Batch1               Batch1
## Sequencing_ID                       JK1                JK10                JK12                JK13                JK14 
## Disease_state                CONTROLLED          CONTROLLED          CONTROLLED          CONTROLLED           CONTROLLED
## Sample_name          Batch1_Ctrl.JK1_S1 Batch1_Ctrl.JK10_S9 Batch1_Ctrl.JK12_S3 Batch1_Ctrl.JK13_S4 Batch1_Ctrl.JK14_S11
## Gender                           Female              Female              Female              Female               Female
## Race                              White               Black               Black               Black                White
## Race_3groups                      White               Black               Black               Black                White
## Age                                  66                  80                  62                  67                   81
## Age_2groups                         Old                 Old               Adult                 Old                  Old
## Smoking_status                   former              former               never              former                  n/a
## Pack_years                            7                  15                   0                  30                  n/a
## Smoking_3groups                   Light               Light               Never               Heavy                  n/a
## EGFRmutation_5groups         T790M-like      Classical-like      Classical-like              Others       Classical-like
## EGFR_mutation               L858R T790M              Del 19       E746_A750del                P741A               Del 19 


DF_A1 = DF_A[-c(1:16),]
myHeader = DF_A[c(1:16),]
dim(DF_A1)
dim(myHeader)
## DF_A1[1:23,1:5]
## myHeader[,1:5]


myGroups = myHeader[7, ]
length( myGroups )
myGroups



myGroups = as.factor( as.vector(unlist(myGroups) ) )
##  myGroups
length( myGroups )
table(myGroups)

myColors <- brewer.pal(6, "Set2")
myColors[myGroups]


## DF_A1[1:5,1:3]
dim( DF_A1 )
myRPKM2 <- matrix(as.numeric(unlist(DF_A1)), ncol = ncol(DF_A1))   # Convert to numeric matrix 
dim( myRPKM2 )
## myRPKM2[1:3,1:3]
colnames(myRPKM2) = colnames(DF_A1)
rownames(myRPKM2) = rownames(DF_A1)
## hist( as.vector(myRPKM2) , breaks=1000)
## hist( as.vector(myRPKM2) , breaks=1000, xlim=c(0, 40))
## hist(  myRPKM2[1,]  , breaks=1000)







 




############################################################################################################################################
AllResults_g <- outDir_g
if( ! file.exists(AllResults_g) ) { dir.create(path=AllResults_g, recursive = TRUE) }


bool3 = ( myGroups == names(table(myGroups))[1] )
bool4 = ( myGroups == names(table(myGroups))[2] )
length( bool3 )
length( bool4 )
length( bool3[bool3] )
length( bool4[bool4] )

myGroup1 = myRPKM2[ , bool3]
myGroup2 = myRPKM2[ , bool4]
dim(myGroup1)
dim(myGroup2)

myHeader1 = myHeader[ , bool3]
myHeader2 = myHeader[ , bool4]
dim(myHeader1)
dim(myHeader2)

write.table( rbind(myHeader1 , myGroup1),    file = paste(AllResults_g, "1a.Group1.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE )

write.table(rbind(myHeader2 , myGroup2),    file = paste(AllResults_g, "1b.Group2.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE )


pvalues = vector( length = nrow(myRPKM2) )
for(i in c(1:nrow(myRPKM2)) ) {  
  x1  =  as.numeric( as.vector( myGroup1[i,] ) ) 
  y1  =  as.numeric( as.vector( myGroup2[i,] ) ) 
  re1 = wilcox.test(x=x1, y = y1,  alternative = "two.sided" )
  pvalues[i] = re1$p.value
}

length( pvalues[ pvalues < 0.05] )
length( pvalues[ pvalues < 0.001] )
pvalues[1]
pvalues[nrow(myRPKM2)]
FDR = p.adjust( pvalues , method = "fdr", n = length(pvalues) )  ## # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",  "fdr", "none")
length( FDR[ FDR < 0.2] )

matrix6 <- cbind(pvalues, FDR)
rownames(matrix6) = rownames( myRPKM2 )

write.table(matrix6,    file = paste(AllResults_g, "2.wilcox.test.FDR.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE )






############################################### DPs
print("ROTS call......")

myGroups2 = as.vector(myGroups) 
myGroups2[bool3] = 1
myGroups2[bool4] = 2
myGroups2
rots_out <- ROTS(data= as.matrix(myRPKM2) , groups=myGroups2, B = 1000,   paired = FALSE,   seed = 14,   log = FALSE, progress = TRUE, verbose = TRUE)
names(rots_out)
save(rots_out, file = paste(AllResults_g, "6a.rots_out.RData", sep="/")   )
 
write.table(rots_out$data,    file = paste(AllResults_g, "6b.rots_out.data.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE )

length(rots_out$pvalue[rots_out$pvalue<0.05])
length(rots_out$FDR[rots_out$FDR<0.2])
rots_re <- cbind( rots_out$d, rots_out$logfc, rots_out$pvalue, rots_out$FDR)
dim(rots_re) 
rots_re[1:5,]
colnames(rots_re) <- c( " ROTS-statistic",  "ROTS.logFC", "ROTS.p-value",  "ROTS.FDR"  )

write.table(rots_re,    file = paste(AllResults_g, "6c.rots_results.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE )

 
   
bool8  <-  (rownames(matrix6) ==  rownames(rots_re )  )
length(bool8)
length(bool8[bool8])

rots_re2 = cbind(rots_re, matrix6)   
rots_re3 = cbind(myRPKM2, rots_re2 )   
dim( rots_re2 )
dim( rots_re3 )

write.table(rots_re2,    file = paste(AllResults_g, "6d.rots_results.wilcox.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE )


write.table(rots_re3,    file = paste(AllResults_g, "6e.rots_results.withMatrix.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE )


 

dim( myHeader )
temp1 = colnames( rots_re2 )
temp2 = rbind( temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1, temp1  )
dim(temp2)
colnames(temp2) = temp1
rownames(temp2) = rownames(myHeader)
temp2

myHeader2 = cbind(myHeader, temp2)
dim( myHeader2 )
dim( rots_re3 )


rots_re4 = rbind(myHeader2, rots_re3 )
dim(rots_re4)   


write.table(rots_re4,    file = paste(AllResults_g, "7.all.Matrix.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE )











   
