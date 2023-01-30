## example:  Rscript   1_filterChromXYM.R     0_rawPeaks     1_filterChromXYM   



args <- commandArgs(TRUE)
print("args: ")
print(args[1])   
print(args[2])   
print("#############")

inputDir= args[1];
outDir  = args[2];     
#  inputDir= "2-noBlacklist";
#  outDir  = "3-filtered";     


outPath1   = paste(outDir, "1-noXYM", sep="/" )
outPath2   = paste(outDir, "2-X", sep="/" )
outPath3   = paste(outDir, "3-Y", sep="/" )
outPath4   = paste(outDir, "4-M", sep="/" )
outPath5   = paste(outDir, "5-runLog", sep="/" )

if( ! file.exists(outPath1) ) { dir.create(outPath1, recursive = TRUE)  }
if( ! file.exists(outPath2) ) { dir.create(outPath2, recursive = TRUE)  }       
if( ! file.exists(outPath3) ) { dir.create(outPath3, recursive = TRUE)  }
if( ! file.exists(outPath4) ) { dir.create(outPath4, recursive = TRUE)  }
if( ! file.exists(outPath5) ) { dir.create(outPath5, recursive = TRUE)  }

peakFiles          <- list.files(path = inputDir, pattern = ".bed$",  full.names = TRUE  )
peakFiles_onlyName <- list.files(path = inputDir, pattern = ".bed$",  full.names = FALSE )


write.table(x=cbind(peakFiles, peakFiles_onlyName), file = "1_filterChromXYM.runLog.txt", 
                append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = FALSE )



###########################################################################################
for(i in c(1:length(peakFiles)) ) { 


matrix_1 <- read.table( peakFiles[i]  , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
myChrom = matrix_1[,1]
matrix_1A <- matrix_1[(myChrom!="chrX" & myChrom!="chrY"  & myChrom!="chrM" ), ] 
matrix_1B <- matrix_1[myChrom=="chrX", ]
matrix_1C <- matrix_1[myChrom=="chrY", ] 
matrix_1D <- matrix_1[myChrom=="chrM", ] 


write.table(matrix_1A, file = paste(outPath1, peakFiles_onlyName[i], sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =  FALSE,   col.names =  FALSE)
write.table(matrix_1B, file = paste(outPath2, peakFiles_onlyName[i], sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =  FALSE,   col.names =  FALSE)
write.table(matrix_1C, file = paste(outPath3, peakFiles_onlyName[i], sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =  FALSE,   col.names =  FALSE)
write.table(matrix_1D, file = paste(outPath4, peakFiles_onlyName[i], sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =  FALSE,   col.names =  FALSE)



sink( file = paste(outPath5, "/",  "1-dimensions-", peakFiles_onlyName[i], sep="")  )
print("raw:")
print( dim(matrix_1) )
print("rmXY:")
print( dim(matrix_1A) )
print("onlyXY:")
print( dim(matrix_1B) )
print("onlyX:")
print( dim(matrix_1C) )
print("onlyY:")
print( dim(matrix_1D) )
sink() 



}



