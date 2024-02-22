#!/usr/bin/env perl
use warnings;
use strict;

die "\n\tperl parse_script.pl web_summary.html out\n\n" if @ARGV != 2;

my $in = shift @ARGV;
my $out = shift @ARGV;

##提取饱和度数据
`cp $in xx`;
`sed -i 's#data#\\n#g' xx`;
`grep 'This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth' xx >yy`;

open IN,"yy" or die $!;
open OUT,">$out" or die $!;
while(<IN>){
	chomp;
	if(/\"x\":\[(\S+)\],\"y\":\[(\S+)\],/){
		my $xaxis = $1;
		my $yaxis = $2;
		print OUT "Saturation_x\t$xaxis\n";
		print OUT "Saturation_y\t$yaxis\n"
	}
}
close IN;

##提取中位基因数据
`grep 'This plot shows the Median Genes per Cell as a function of downsampled sequencing depth in mean reads per cell' xx >yy`;
open IN,"yy" or die $!;
while(<IN>){
	chomp;
	if(/\"x\":\[(\S+)\],\"y\":\[(\S+)\],/){
		my $xaxis = $1;
		my $yaxis = $2;
		print OUT "Median_Genes_per_Cell_x\t$xaxis\n";
		print OUT "Median_Genes_per_Cell_y\t$yaxis\n"
	}
}
close IN;
close OUT;

##清理中间数据
`rm xx yy`;
