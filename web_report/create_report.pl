#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my %opts;
GetOptions(\%opts,"i=s","c=s","h");
if(!defined($opts{i}) || !defined($opts{c}) ||defined($opts{h})){
	print <<"	Usage End.";
Description:

Usage:
    -i        input dir         <dir>     must be given
    -c        config            <file>    must be given

	-h        Help document

	Usage End.
	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
&show_log("Start Time :[$Time_Start]");
################
$| = 1 ;
## get parameters
my $cfgfile = abs_path($opts{c});
my $outdir = abs_path($opts{i});
## reading cfg file
my %hcfg = ();
&reading_cfg_file($cfgfile, \%hcfg);


my $result_dir = $hcfg{OUTDIR} ;
my $prefix = $hcfg{PREFIX} ;
##10x web report
mkdir "$outdir/10x" unless -d "$outdir/10x";
my $dir = "$result_dir/02.cellranger/$prefix/outs";
`cp -r $dir/filtered_feature_bc_matrix $dir/analysis $outdir/10x`;
`cp $dir/metrics_summary.csv $dir/web_summary.html $outdir/10x`;

##bmk web report
mkdir "$outdir/bmk" unless -d "$outdir/bmk";
#`cp -r $dir/filtered_feature_bc_matrix $dir/raw_feature_bc_matrix $dir/metrics_summary.csv $outdir/bmk`;
`cp $dir/metrics_summary.csv $outdir/bmk`;
`cp -r $dir/filtered_feature_bc_matrix $outdir/bmk/$prefix.filtered_feature_bc_matrix`;
`cp -r $dir/raw_feature_bc_matrix $outdir/bmk/$prefix.raw_feature_bc_matrix`;
`cp $result_dir/03.cluster/umap_df.xls $result_dir/03.cluster/tsne_df.xls $outdir/bmk`;
#my $cmd = "cd $outdir/bmk && $python_env python $Bin/../cellcalling/rank_plot.py -m $outdir/bmk/$prefix.raw_feature_bc_matrix -f $outdir/bmk/$prefix.filtered_feature_bc_matrix\n";
#&run_or_die($cmd) ;

mkdir "$outdir/bmk/tmp" unless -d "$outdir/bmk/tmp";
`cp $Bin/* $outdir/bmk/tmp`;
#`cp $outdir/bmk/metrics_summary.csv $outdir/bmk/Qual.stat $outdir/bmk/rank_barcode_umi.html $outdir/bmk/tmp`;
`cp $outdir/bmk/metrics_summary.csv $outdir/bmk/tmp`;
`cp $result_dir/01.fastq2BcUmiSC/$prefix.qual.stat $outdir/bmk/tmp/Qual.stat`;
`cp $result_dir/03.cluster/umap_df.xls $result_dir/03.cluster/tsne_df.xls $result_dir/03.cluster/marker_gene.csv $outdir/bmk/tmp`;
`cp -r $dir/filtered_feature_bc_matrix $dir/raw_feature_bc_matrix $outdir/bmk/tmp`;
my $cmd = "perl $Bin/parse_script.pl $dir/web_summary.html $outdir/bmk/tmp/Saturation.xls\n";
&run_or_die($cmd) ;

my $R_env = "";
$R_env = "export PATH=".$hcfg{Rscript}.":\$PATH && " if defined $hcfg{Rscript};
$cmd = "cd $outdir/bmk/tmp && $R_env Rscript -e \"rmarkdown::render(input = 'report_v3.Rmd',output_file = 'report.html')\" && mv $outdir/bmk/tmp/report.html $outdir/bmk/$prefix.report.html\n";
&run_or_die($cmd) ;

#####修改报告中report
`sed -i 's/样品信息:/样品信息: $prefix/' $outdir/bmk/$prefix.report.html`;
`sed -i 's/质控报告/质控报告-$prefix/' $outdir/bmk/$prefix.report.html`;
##########结果汇总提交产品
&stat();
`R --slave < $outdir/bmk/$prefix.R`;

#####清理R脚本和临时数据
`rm $outdir/bmk/$prefix.R`;
`rm -r $outdir/bmk/tmp`;
#########for产品用于数据汇总
sub stat{
	open (OUT,">$outdir/bmk/$prefix.R") or die "$!";
	print OUT << "	R.end";
	setwd("$outdir/bmk/")
	bc = read.table("tmp/Qual.stat",sep ="\t")
	metrics <- read.table("metrics_summary.csv",sep = ",",header = T,check.names = F)
	bc\$V4 = paste0(bc\$V2,"(",bc\$V3,")")
	bc_sub = bc[,c("V4")]
	names(bc_sub) <- bc\$V1
	bc_sub <- t(as.data.frame(bc_sub))
	df = cbind(bc_sub,metrics)
	write.table(x = df,file = "$prefix.stat.xls",row.names= F,quote=F,sep = "\t",col.names=T)
	q("no")
	R.end
	close(OUT);
	return("$outdir/bmk/$prefix.R");
}


###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&show_log("End Time :[$Time_End]");


############################################
sub reading_cfg_file()
{
    my ($cfgfile, $ahcfg) = @_ ;

    open (IN, $cfgfile) || die "$cfgfile, $!\n" ;
    while(<IN>){
        chomp ;
        next if (m/^\#|^\s*$/);
        my ($key, $value) = split /\s+/ ;
        $ahcfg->{$key} = $value ;
    }
    close(IN);
    return ;
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
        $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub show_log()
{
        my ($txt) = @_ ;
        my $time = time();
        my $Time = &sub_format_datetime(localtime($time));
        print "$Time:\t$txt\n" ;
        return ($time) ;
}

sub run_or_die()
{
        my ($cmd) = @_ ;
        my $start_time = &show_log($cmd);
        my $flag = system($cmd) ;
        if ($flag != 0){
                my $end_time = &show_log("Error: command fail: $cmd");
                &show_log("Elaseped time: ".($end_time - $start_time)."s\n");
                exit(1);
        }
        my $end_time = &show_log("done.");
        &show_log("Elaseped time: ".($end_time - $start_time)."s\n");
        return ;
}
