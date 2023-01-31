## Bar Plots with gapped y-axis for the number of reads
## suppressPackageStartupMessages( library(plotrix)  )



## Bar Plots for the number of reads
df_1 <- read.table("narrow_run1.number.txt",  header=F,   sep="\t" )  
dim( df_1 )
df_1[1:5 , ]

num_peaks <- as.numeric( as.vector( df_1[ , 1] ) )
length(num_peaks)
num_peaks[1:20]

num_peaks2 = sort( num_peaks, decreasing = TRUE)  
num_peaks3 = sort( num_peaks[1:30], decreasing = TRUE)  

pdf( "bar_plot.all.pdf"  , width= 15, height=3)
    barplot(  num_peaks2, col="red" )
dev.off()

pdf( "bar_plot.bottom30.pdf"  , width= 7, height=3)
    p1 <- barplot(  num_peaks3, col="red" )
    text(x = p1, y = num_peaks3/1.1 ,   labels = num_peaks3, srt=90, col="blue" )
dev.off()





