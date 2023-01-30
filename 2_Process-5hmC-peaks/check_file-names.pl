#!/usr/bin/env perl
use strict;
use warnings;





########## Keys and Values ##########
my %args = @ARGV;



########## Set Defaults ##########
my $file1_g      = 'All.txt';
my $file2_g      = 'merged.bed';
my $file3_g      = 'T.genotype.txt';

########## Get Arguments ##########
if ( exists $args{'-file1'}   )    { $file1_g  = $args{'-file1'};   }
if ( exists $args{'-file2'}   )    { $file2_g  = $args{'-file2'};   }
if ( exists $args{'-file3'}   )    { $file3_g  = $args{'-file3'};   }

########### Conditions #############
$file1_g      =~ m/^\S+$/ or die;
$file2_g      =~ m/^\S+$/ or die;
$file3_g      =~ m/^\S+$/ or die;


######### Example ###########
# perl   check_file-names.pl   -file1  1_297BW_from_deepTools/0_Bigwig      -file2   2_Peaks_297/0_rawPeaks      -file3 3_297BAM_Batch_Effects/0_BAM






opendir(INPUT1,     $file1_g)     or die "$!"; 
opendir(INPUT2,     $file2_g)     or die "$!"; 
opendir(INPUT3,     $file3_g)     or die "$!"; 
my @files1 = readdir(INPUT1);
my @files2 = readdir(INPUT2);
my @files3 = readdir(INPUT3);

for (my $i=0; $i<=$#files1; $i++) {
    next unless $files1[$i] =~ s/\.bw//;
    print("...... $files1[$i]\n");
    my $peakFile = "$files1[$i].bed";
    my $bamFile  = "$files1[$i].bam";
    my $bool_peak = 0; 
    my $bool_bam  = 0; 
    
    for (my $j=0; $j<=$#files2; $j++) {
        if($files2[$j] eq $peakFile){$bool_peak = 1; }
    }
    for (my $k=0; $k<=$#files3; $k++) {
        if($files3[$k] eq $bamFile){$bool_bam = 1; }
    }

    if($bool_peak == 0) {print "peak: $files1[$i]\n"; }
    if($bool_bam  == 0) {print "bam:  $files1[$i]\n"; }
}




