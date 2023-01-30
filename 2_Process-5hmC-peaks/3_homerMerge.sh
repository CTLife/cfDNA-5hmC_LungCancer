out="3_Merged"
mkdir  -p  $out

mergePeaks  -strand  -d given   2_pooled/pooled.bed  -prefix  $out       > $out/allPeaks.bed  


