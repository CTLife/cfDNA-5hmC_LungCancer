## plot mapping ratio and error rate

suppressPackageStartupMessages( library(ggplot2)  )

df_1 <- read.table("multiqc_general_stats.txt",  header=T,   sep="\t" )  
dim( df_1 )
df_1[1:5 , ]

samples       <- as.vector( df_1[ , 1] )
mapping_ratio <- as.numeric( as.vector( df_1[ , 2] ) )
error_rate    <- as.numeric( as.vector( df_1[ , 3] ) )

length(samples)
length(mapping_ratio)
length(error_rate)

samples[1:5]
mapping_ratio[1:5]
error_rate[1:5]





#df_2 <- data.frame(x=c(samples,samples),  y=c(mapping_ratio, error_rate), 
#                   type=c( rep("mapping_ratio", length=length(mapping_ratio)), rep("error_rate", length=length(error_rate)) ) )
#dim( df_2 )
#df_2[1:5 , ]


df_2A <- data.frame( x=samples,  y=mapping_ratio  )
dim( df_2A )
df_2A[1:5 , ]

df_2B <- data.frame( x=samples,  y=error_rate  )
dim( df_2B )
df_2B[1:5 , ]


pdf( "1.pdf" , width= 10, height=3 )
ggplot(df_2A, aes(x=x, y=y)) + geom_point(size=1, shape=19, color="red")
ggplot(df_2B, aes(x=x, y=y)) + geom_point(size=1, shape=19, color="cyan")
dev.off()



