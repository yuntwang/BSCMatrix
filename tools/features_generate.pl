#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2021
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2021
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2021
my $ver="1.0";

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#####################

my %opts;
GetOptions(\%opts,"i=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-i           gtf file                      <infile>     must be given
		-o           features file                 <outfile>    must be given
		-h           Help document

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
my $gtffile = $opts{i} ;
my $outfile = $opts{o} ;

## process
&reading_and_generate_features_file($gtffile, $outfile);

###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
&show_log("End Time :[$Time_End]");

###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir [$in]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
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

## qsub
sub qsub()
{
	my ($shfile, $maxproc) = @_ ;
	my $dir = dirname($shfile) ;
	my $cmd = "cd $dir && sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --reqsub --independent $shfile" ;
	&run_or_die($cmd);

	return ;
}

## qsub_mult($shfile, $max_proc, $job_num)
sub qsub_mult()
{
	my ($shfile, $max_proc, $job_num) = @_ ;
	if ($job_num > 500){
		my @shfiles = &cut_shfile($shfile);
		for my $file (@shfiles){
			&qsub($file, $max_proc);
		}
	}
	else{
		&qsub($shfile, $max_proc) ;
	}
}

#my @shfiles = &cut_shfile($shfile);
sub cut_shfile()
{
	my ($file) = @_ ;
	my @files = ();
	my $num = 0 ;
	open (IN, $file) || die "Can't open $file, $!\n" ;
	(my $outprefix = $file) =~ s/.sh$// ;
	while (<IN>){
		chomp ;
		if ($num % 500 == 0){
			close(OUT);
			my $outfile = "$outprefix.sub_".(int($num/500)+1).".sh" ;
			push @files, $outfile ;
			open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
		}
		print OUT $_, "\n" ;
		$num ++ ;
	}
	close(IN);
	close(OUT);

	return @files ;
}

sub sub_normal_dis_ran(){#$average,$standard_deviation
        my ($aver,$stand_dev) = @_ ;
        my $ran1 = rand() ;
        my $ran2 = rand() ;
        my $y = ((-2*log($ran1))**(0.5))*sin(2*3.1415926535*$ran2) ;
        return ($aver + $stand_dev*$y) ;
}

#&reading_and_generate_features_file($gtffile, $outfile);
sub reading_and_generate_features_file()
{
    my ($gtffile, $outfile) = @_ ;

    open (IN, $gtffile) || die "$gtffile, $!\n" ;
    open (OUT, ">$outfile") || die "$outfile, $!\n" ;
    my %ht = ();
    while(<IN>){
        chomp ;
        next if (m/^\#|^\s*$/);
        my ($chr, undef, $type, $start, $end, undef, $ori, $phase, $info) = split /\t+/, $_, 9;
        if ($type eq 'gene' || $type eq 'exon'){
            #my %hinfo = @infos ;
            my @infos = split /\;\s*/, $info ;
            my %hinfo = ();
            for (my $i=0; $i<@infos; $i++){
                my ($key, $value) = split /\s+/, $infos[$i], 2 ;
                $hinfo{$key} = $value ;
            }
            next if (!defined $hinfo{'gene_id'});
            my $gene_id = $hinfo{'gene_id'} ;
            $gene_id =~ s/\"|\;//g ;
            $gene_id =~ s/^gene[\:\-\_]// ;
            next if ($gene_id =~ m/^\s*$/) ;
            if (defined $ht{$gene_id}){
                next ;
            }
            else{
                $ht{$gene_id} = 1;
            }
            if (defined $hinfo{'gene_name'}){
                my $gene_name = $hinfo{'gene_name'} ;
                $gene_name =~ s/\"|\;//g ;
                print OUT "$gene_id\t$gene_name\tGene Expression\n"
            }
            else{
                print OUT "$gene_id\t$gene_id\tGene Expression\n" ;
            }
        }
    }
    close(IN);
    close(OUT);

    return ;
}

