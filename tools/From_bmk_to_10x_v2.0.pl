#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(shuffle);
use threads;
use threads::shared;

my %opts;
GetOptions(\%opts,"b=s","s=s","c=s","h");
if(!defined($opts{b}) || !defined($opts{s}) || !defined($opts{c}) ||defined($opts{h})){
	print <<"	Usage End.";
	
	-b        barcode info        [bmk_10x_barcode.xls]
	-s        bc_umi_read info    [*.umi]
	-c        config

	-h        Help document
	
	
	Usage End.
	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
&show_log("Start Time :[$Time_Start]");

$| = 1;

my $barcode = $opts{b};
my $bc_umi_r = $opts{s};
my $config = $opts{c};
###########################读取配置文件
open IN,$config or die $!;
my ($Fq1,$Fq2,$prefix,$fa,$gtf,$ram,$EC,$threads,$Nobam);
while(<IN>){
	chomp;
	next if /^\s*$/;
	next if /^\#/;
	my ($key,$value) = split;
	$Fq1 = $value if /^FQ1/;
	$Fq2 = $value if /^FQ2/;
	$threads = $value if /^Threads/;
	$prefix = $value if /^PREFIX/;
	$fa = $value if /^FA/;
	$gtf = $value if /^GTF/;
	$ram = $value if /^RAM/;
	$EC = $value if /^EC/;
	$Nobam = $value if /^Nobam/;
}
close IN;

############识别测序仪器信息(兼容压缩和非压缩数据),兼容BGI数据
my ($info,$bgi);
if($Fq1 =~ /.gz$/){
	my $line = `zcat $Fq1 | head -n 1`;
	chomp $line;
	if($line =~ /\/1/){
		$bgi = "yes";
	}else{
		$info = (split /\s+/,$line)[1];
		$bgi = "no";
	}
}else{
	my $line = `head -n 1 $Fq1`;
	chomp $line;
	if($line =~ /\/1/){
		$bgi = "yes";
	}else{
		$info = (split /\s+/,$line)[1];
		$bgi = "no";
	}
}

########################统计 10x barcode vs bmk的barcode对照表
my %table;
open IN,"$barcode" or die $!;
while(<IN>){
	chomp;
	next if /^\s*$/;
	my ($bmk,$bc10x) = split;
	$table{$bmk} = $bc10x;
}
close IN;

##########################################输出R1的信息(多线程并行化处理)
# 计算文件总行数
my $total_num = `awk 'END{print NR}' $bc_umi_r`;
# 每个子文件的行数
my $size = (int($total_num / ($threads * 4)) + 1)*4;
##########快速拆分*.umi数据
`awk -v sz=$size  'BEGIN{i=1}{ print \$0 > "$prefix." i ".umi"; if (NR>=i*sz){close("$prefix." i ".umi");i++}}' $bc_umi_r`;

#######################并行运行,确保捕获的umi数据是按照顺序来的
my @files;
foreach(1..$threads){
	push @files,"$prefix.$_.umi";
}
# 开始多线程
#my @atcg = qw /A T C G/ x 23;
#my @fred1 = qw /F F F F F F F F F F F F F F F F F F F F F F F F F F F F ; :/;
#my @fred2 =  ('F','F','F',':',',',';') x 20;
foreach (0..$threads-1){
	my $thr = threads -> create(\&dealumi, $files[$_], $_);
}

# join各个线程的结果
my @num_array2;
while(threads -> list()){
    foreach my $t(threads -> list(threads::joinable)){
        my $tmp = $t -> join();
        push @num_array2,$tmp;
    }
}
print "threads:@num_array2\n";
unlink @files;

##############并行结果整理合并,并清理中间结果
`rm R1.fq.gz` if -e "R1.fq.gz";
`rm ReadID.xls` if -e "ReadID.xls";
foreach (0..$threads-1){
	`cat R1.$_.fq.gz >>R1.fq.gz`;
	`cat ReadID.$_.xls >>ReadID.xls`;
	`rm R1.$_.fq.gz ReadID.$_.xls`;
}

##################################################R1 和 R2 readID 匹配
my $cmd = "seqkit grep -f ReadID.xls $Fq2 -o R2.fq.gz -j $threads\n";
&run_or_die($cmd);

####################################整理数据,用于cellranger分析
if(-d "data"){
	`rm -rf data`;
	mkdir("data",0755) unless -d "data";
	`ln R1.fq.gz data/$prefix\_S1_L001_R1_001.fastq.gz`;
	`ln R2.fq.gz data/$prefix\_S1_L001_R2_001.fastq.gz`;
}else{
	mkdir("data",0755) unless -d "data";
	`ln R1.fq.gz data/$prefix\_S1_L001_R1_001.fastq.gz`;
	`ln R2.fq.gz data/$prefix\_S1_L001_R2_001.fastq.gz`;
}

#######################################执行cellranger分析
### create genome index and cellranger count
if(-e "INDEX/reference.json"){
	if($Nobam == 0){
		`rm -rf $prefix` if -d $prefix;
		my $cmd = "cellranger count --id $prefix --transcriptome INDEX --fastqs data  --localcores=$threads --localmem=$ram --include-introns true --expect-cells $EC --no-bam\n";
		&run_or_die($cmd);
	}else{
		`rm -rf $prefix` if -d $prefix;
		my $cmd = "cellranger count --id $prefix --transcriptome INDEX --fastqs data  --localcores=$threads --localmem=$ram --include-introns true --expect-cells $EC\n";
		&run_or_die($cmd);
	}
}else{
	`rm -rf INDEX` if -d "INDEX";
	my $cmd = "cellranger mkref --genome=INDEX --fasta=$fa --genes=$gtf --nthreads=$threads --memgb=$ram\n";
	&run_or_die($cmd);
	if($Nobam == 0){
		`rm -rf $prefix` if -d $prefix;
		$cmd = "cellranger count --id $prefix --transcriptome INDEX --fastqs data  --localcores=$threads --localmem=$ram --include-introns true --expect-cells $EC --no-bam\n";
		&run_or_die($cmd);
	}else{
		`rm -rf $prefix` if -d $prefix;
		$cmd = "cellranger count --id $prefix --transcriptome INDEX --fastqs data  --localcores=$threads --localmem=$ram --include-introns true --expect-cells $EC\n";
		&run_or_die($cmd);
	}
}



###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&show_log("End Time :[$Time_End]");
################################sub
sub rand_seq(){
	my $i = 0;
	my @dna;
	while($i < 92){
		my @seq = qw/A C G T/;
		my $base = $seq[rand(@seq)];
		push @dna,$base;
		$i++;
	}
	my $seq = join "",@dna;
	return($seq);
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

sub dealumi()
{
	 my ($umi, $num) = @_ ;
	 open IN,$umi or die $!;
	 open OUT,"| gzip >R1.$num.fq.gz";
	 open R1ID,">ReadID.$num.xls" or die $!;
	 while(<IN>){
		chomp;
		next if /^\s*$/;
		next if /^\#/;
		my ($read,$bc,$cigar,$bcpos,$umi,$umipos) = split;
		my $bc10x;
		if($table{$bc}){
			$bc10x = $table{$bc};
		}else{
			die "please check $bc is not in $barcode!!!\n";
		}
		#my $polyT = "T"x 30;
		#my $atcg_shuf = join "",shuffle @atcg;
		#my $fred1_shuf = join "",shuffle @fred1;
		#my $fred2_shuf = join "",shuffle @fred2;
		#my $seq = "$bc10x"."$umi"."$polyT"."$atcg_shuf";
		my $seq = "$bc10x"."$umi";
		#my $qual = "$fred1_shuf"."$fred2_shuf";
		my $qual = "FFFFFFFFFFFFFFFFFFFFFFFFFFFF";
		if($bgi eq "yes"){
			print OUT "\@${read}\/1\n$seq\n+\n$qual\n";
			print R1ID "$read\/2\n";
		}else{
			print OUT "\@$read $info\n$seq\n+\n$qual\n";
			print R1ID "$read\n";
		}
	 }
	close IN;
	close OUT;
	close R1ID;
	return($num);
}


