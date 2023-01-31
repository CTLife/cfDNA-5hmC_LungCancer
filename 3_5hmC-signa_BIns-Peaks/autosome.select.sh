grep -P "(^chr\d)|(^\#)"     D.1-RawCounts.2000bp-Bin.txt   > Autosome.FPKM.txt

cut -f 1  Autosome.FPKM.txt  | sort | uniq




