#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2021
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2021
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2021
my $ver="2.1";

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

#####################

my %opts;
GetOptions(\%opts,"i=s","k=s","b=s","d=s","p=s","s=s","l=s","u=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{b}) || !defined($opts{k}) || !defined($opts{d}) || !defined($opts{p}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver
		v1.2:   out put full barcode and umi fq reads for following analysis.
		v1.3:   only out put valid read id instead of fq file.
		v2.1:   change the way of barcode search.
		v2.2:   adapt to multi-threaded version


	Usage:

		-i           reads fq1 file                          <infile>     must be given
		-k           kmer file                               <infile>     must be given
		-b           bc file                                 <infile>     must be given
		-d           result dir                              <outdir>     must be given
		-p           outfile prefix                          <string>     must be given
		-s           max offset number                       [int]        optional[5]  
		-l           kmer length                             [int]        optional[6]  
		-u           umi length                              [int]        optional[12]  
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
my $read1file = $opts{i} ;
my $kmerfile = $opts{k} ;
my $bcfile = $opts{b} ;
my $outdir = $opts{d} ;
mkdir $outdir if (!-d $outdir) ;
$outdir = &ABSOLUTE_DIR($outdir) ;
my $prefix = $opts{p} ;
my $max_offset = defined $opts{s} ? $opts{s} : 5 ;
my $kmerlen = defined $opts{l} ? $opts{l} : 6 ;
my $umi_len = defined $opts{u} ? $opts{u} : 12 ;

my $link1 = "ACTGGTGA" ;
my $link2 = "GGTAGTGACA" ;
my %humi_pos = (
    'bc1' => 57,
    'bc2' => 34,
    'bc3' => 12,
);

my %href_pos = (
	'bc1'   => 0 ,
	'bc2'   => 23 ,
	'bc3'   => 45 ,
);
my %hbc_type_len = (
	'bc1'   => 16 ,
	'bc2'   => 12 ,
	'bc3'   => 12 ,
);
my %hbc_range = (
    'bc1' => [0, 21],
    'bc2' => [18, 39],
    'bc3' => [40, 63],
);

## get kmer
my %hkmer = ();
my %hkmer_pos = ();
&get_kmer($kmerfile, \%hkmer, \%hkmer_pos);

## get bc seq
my %hbc_seq = ();
my %hbc_len = ();
my %hbcs = ();
&get_bc_seq($bcfile, \%hbc_seq, \%hbc_len, \%hbcs);

## process
my %hbc_umi = ();
&read_fq_and_check_barcode($read1file, $outdir, $prefix, $max_offset, \%hbc_seq, \%hbc_len, \%hbcs, \%hkmer, \%hkmer_pos, $kmerlen, \%humi_pos, \%href_pos, $umi_len, \%hbc_umi);

## out put bc umi
#&output_bc_umi(\%hbc_umi, $outdir, $prefix);

##############Time
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

#&get_kmer($kmerfile, \%hkmer, \%hpos);
sub get_kmer()
{
	my ($kmerfile, $ahkmer, $ahpos) = @_ ;

	open (IN, $kmerfile) || die "$kmerfile, $!\n" ;
	while(<IN>){
		chomp ;
		next if (m/^\s*$|^\#/);
		my ($kmer, $num, @ids) = split ;
		for (my $i=0; $i<@ids; $i+=2){
			push @{$ahkmer->{$kmer}}, $ids[$i] ;
			push @{$ahpos->{$kmer}}, $ids[$i+1] ;
		}
	}
	close(IN);

	return ;
}

#&get_bc_seq($bcfile, \%hbc_seq, \%hbc_len, \%hbcs);
sub get_bc_seq()
{
    my ($bcfile, $ahbc_seq, $ahbc_len, $ahbcs) = @_ ;

    open (IN, $bcfile) || die "$bcfile, $!\n" ;
    $/ = ">" ;
    while(<IN>){
        chomp ;
        next if (m/^\s*$/);
        my ($id, $seq) = split(/\n/, $_, 2) ;
        $id = (split /\s+/, $id)[0] ;
        $seq =~ s/\n//g ;
        my $len = length($seq) ;
        $ahbc_seq->{$id} = $seq ;
        $ahbc_len->{$id} = $len ;
        $ahbcs->{$seq} = $id ;
    }
    $/ = "\n" ;
    close(IN);

    return ;
}

#&read_fq_and_check_barcode($read1file, $outdir, $prefix, $max_offset, \%hbc_seq, \%hbc_len, \%hbcs, \%hkmer, \%hkmer_pos, $kmerlen, \%humi_pos, \%href_pos, $umi_len, \%hbc_umi);
sub read_fq_and_check_barcode()
{
	my ($fqfile, $outdir, $prefix, $max_offset, $ahbc_seq, $ahbc_len, $ahbcs, $ahkmer, $ahkmer_pos, $kmerlen, $ahumi_pos, $ahref_pos, $umi_len, $ahbc_umi) = @_ ;

    my $str = $fqfile ;
    if ($fqfile =~ /.gz$/){
        $str = "gunzip -c $fqfile |" ;
    }
	open (IN, $str) || die "$fqfile, $!\n" ;
	open (OUT, ">$outdir/$prefix.umi") || die "$outdir/$prefix.umi, $!\n" ;
	open (FIT, ">$outdir/$prefix.filter") || die "$outdir/$prefix.filter, $!\n" ;
    #open (ID, ">$outdir/$prefix.select_id") || die "$outdir/$prefix.select_id, $!\n" ;
    #open (MAP, ">$outdir/$prefix.id_map") || die "$outdir/$prefix.id_map, $!\n" ;
    my %hbc_stat = ();
    my %hsingle_stat = ();
    my $total = 0 ;
    my $full = 0 ;
    my $part = 0 ;
    my $null = 0 ;
    my $umi_num = 0 ;
    my $final_num = 0 ;
    my $read_index = 0 ;
	while(<IN>){
		my $id = $_ ; my $seq = <IN> ; my $id2 = <IN> ; my $qual = <IN> ;
		chomp($id, $seq, $id2, $qual);
		$id = (split /\s+/, $id)[0] ;
		$id =~ s/^\@// ;
        $id =~ s/\/1$// ;
        #print MAP "$read_index\t$id\n" ;
        my $r_id = $read_index ;
        $read_index++ ;
        my @types = ();
        my @bcs = ();
        my @cigars = ();
        my @end_poses = ();
        my $bc_flag = 0 ;
        my $final_umi_seq = "" ;
        my $final_umi_pos = "" ;
        my ($bc1, $bc2, $bc3, $umi) = ();
        if ($seq =~ m/([ACGTN]{12,17})$link1([ACGTN]{11,13})$link2([ACGTN]{12})([ACGTN]{$umi_len})/g){
            my $pos = pos($seq) ;
            if (abs($pos - (13+8+12+10+12+$umi_len+2)) < $max_offset){
                my ($bc1_seq, $bc2_seq, $bc3_seq, $umi) = ($1, $2, $3, $4) ;
                $final_umi_pos = $pos - $umi_len + 1;
                $final_umi_seq = $umi ;
                my ($bc1, $bc1_cigar) = &judge_bc1('bc1', $bc1_seq, $ahbcs, $ahbc_seq, $ahkmer, $kmerlen, $ahbc_len);
                my $bc1_pos = "" ;
                if ($bc1 eq 'null'){
                    ($bc1, $bc1_cigar, $bc1_pos, my $umi_pos, my $umi_seq) = &judge_bc_kmer('bc1', $seq, $ahbcs, $ahbc_seq, $ahbc_len, $ahref_pos, $ahkmer_pos, $ahkmer, $kmerlen, $ahumi_pos);
                }
                else{
                    $bc1_pos = ($href_pos{bc1}+$ahbc_len->{$bc1}-1).",".($pos - $umi_len - 12 - length($link2) - length($bc2_seq) - length($link1)) ;
                }
                my ($bc2, $bc2_cigar) = &judge_bc('bc2', $bc2_seq, $ahbcs, $ahbc_seq, $ahkmer, $kmerlen);
                my $bc2_pos = ($href_pos{bc2}+$hbc_type_len{bc2}-1).",".($pos - $umi_len - 12 - length($link2)) ;
                if ($bc2 eq 'null'){
                    ($bc2, $bc2_cigar, $bc2_pos, my $umi_pos, my $umi_seq) = &judge_bc_kmer('bc2', $seq, $ahbcs, $ahbc_seq, $ahbc_len, $ahref_pos, $ahkmer_pos, $ahkmer, $kmerlen, $ahumi_pos);
                }
                my ($bc3, $bc3_cigar) = &judge_bc('bc3', $bc3_seq, $ahbcs, $ahbc_seq, $ahkmer, $kmerlen);
                my $bc3_pos = ($href_pos{bc3}+$hbc_type_len{bc3}-1).",".($pos - $umi_len) ;
                if ($bc3 eq 'null'){
                    ($bc3, $bc3_cigar, $bc3_pos, my $umi_pos, my $umi_seq) = &judge_bc_kmer('bc3', $seq, $ahbcs, $ahbc_seq, $ahbc_len, $ahref_pos, $ahkmer_pos, $ahkmer, $kmerlen, $ahumi_pos);
                    if ($bc3 ne 'null'){
                        $final_umi_pos = $umi_pos ;
                        $final_umi_seq = $umi_seq ;
                    }
                }
                if ($bc1 ne 'null'){
                    $bc_flag = 1 ;
                    push @end_poses, $bc1_pos ;
                    push @types, 'bc1';
                    push @bcs, $bc1 ;
                    push @cigars, $bc1_cigar ;
                }
                if ($bc2 ne 'null'){
                    $bc_flag = 1 ;
                    push @end_poses, $bc2_pos ;
                    push @types, 'bc2';
                    push @bcs, $bc2 ;
                    push @cigars, $bc2_cigar ;
                }
                if ($bc3 ne 'null'){
                    $bc_flag = 1 ;
                    push @end_poses, $bc3_pos ;
                    push @types, 'bc3';
                    push @bcs, $bc3 ;
                    push @cigars, $bc3_cigar ;
                }
            }
        }
        if ($bc_flag == 0){
            my $len = length($seq) ;
            my $offlen = $len > 75 ? 75 : $len ;
            my %htypes = ();
            my %hkmer_ref = ();
            for (my $i=0; $i<$offlen-$kmerlen+1; $i++){
                my $kmer = substr($seq, $i, $kmerlen);
                if (defined $ahkmer->{$kmer}){
                    for (my $j=0; $j<@{$ahkmer->{$kmer}}; $j++){
                        my $ref_pos = $ahkmer_pos->{$kmer}[$j] ;
                        my $map_pos = $i + 1;
                        my $offset = abs($ref_pos - $map_pos) ;
                        if ($offset > $max_offset){
                            next ;
                        }
                        if ($ahkmer->{$kmer}[$j] =~ m/(bc\d)\_(\d+)/){
                            $htypes{$1}{$ahkmer->{$kmer}[$j]}++ ;
                        }
                        push @{$hkmer_ref{$ahkmer->{$kmer}[$j]}}, [$ahkmer_pos->{$kmer}[$j], $map_pos] ;
                    }
                }
            }
            for my $type (sort keys %htypes){
                for my $bc (sort {$htypes{$type}{$b} <=> $htypes{$type}{$a} || $a cmp $b} keys %{$htypes{$type}}){
                    my $bc_seq = $ahbc_seq->{$bc} ;
                    my $bc_len = $ahbc_len->{$bc} ;
                    my @poses = @{$hkmer_ref{$bc}} ;
                    my @bc_reg = ([$poses[0][0], $poses[0][0]]);
                    my @read_reg = ([$poses[0][1], $poses[0][1]]);
                    for (my $i=1; $i<@poses; $i++){
                        my ($ref_pos, $read_pos) = @{$poses[$i]} ;
                        if ($ref_pos - $bc_reg[-1][1] == 1 && $read_pos - $read_reg[-1][1] == 1){
                            $bc_reg[-1][1] = $ref_pos ;
                            $read_reg[-1][1] = $read_pos ;
                        }
                        else{
                            push @bc_reg, [$ref_pos, $ref_pos] ;
                            push @read_reg, [$read_pos, $read_pos] ;
                        }
                    }
                    my ($flag, $match_num, $mismatch_num, $cigar, $bc_end_pos, $read_end_pos) = &judge_match(\@bc_reg, \@read_reg, $bc_seq, $bc_len, $seq, $kmerlen, $ahref_pos->{$type});
                    if ($flag == 1){
                        my $umi_pos = $ahumi_pos->{$type} + ($read_end_pos - $bc_end_pos) + $ahref_pos->{$type} ;
                        my $umi_seq = substr($seq, $umi_pos-1, $umi_len);
                        $final_umi_pos = $umi_pos ;
                        $final_umi_seq = $umi_seq ;
                        push @end_poses, "$bc_end_pos,$read_end_pos" ;
                        push @types, $type ;
                        push @bcs, $bc ;
                        push @cigars, $cigar ;
                        $bc_flag = 1 ;
                        last ;
                    }
                }
            }
        }
        $total++ ;
        if ($bc_flag == 1){
            my $flag = &check_umi_seq($final_umi_seq, $umi_len);
            for my $bc (@types){
                $hsingle_stat{$bc}++ ;
            }
            my $bc_type = join("-", @types);
            $hbc_stat{$bc_type}++ ;
            if ($flag == 0){
                print FIT "$id\t", join("-", @bcs),"\t", join(":",@cigars) ;
                print FIT "\t", join("-", @end_poses) ;
                print FIT "\t$final_umi_seq\t$final_umi_pos\n" ;
                if (scalar @types == 3){
                    $full++ ;
                }
                else{
                    $part++ ;
                }
            }
            elsif (scalar @types == 3 && $bcs[0] =~ /bc1/ && $bcs[1] =~ /bc2/ && $bcs[2] =~ /bc3/){
                my $type = join("-", @bcs) ;
                print OUT "$id\t$type\t", join(":",@cigars) ;
                print OUT "\t", join("-", @end_poses) ;
                print OUT "\t$final_umi_seq\t$final_umi_pos\n" ;
                #push @{$ahbc_umi->{$type}{$final_umi_seq}}, $id ;
                $full++ ;
                $final_num++ ;
                $umi_num++ ;
                #print ID "$r_id\n" ;
            }
            else{
                my $type = join("-", @bcs) ;
                print FIT "$id\t$type\t", join(":",@cigars) ;
                print FIT "\t", join("-", @end_poses) ;
                print FIT "\t$final_umi_seq\t$final_umi_pos\n" ;
                $part++ ;
                $umi_num++ ;
            }
        }
        else{
            print FIT "$id\tNull\n" ;
            $null++ ;
        }
	}
	close(IN);
	close(OUT);
    close(FIT);
    #close(ID);
    #close(MAP);
    open(OUT, ">$outdir/$prefix.stat") || die "$outdir/$prefix.stat, $!\n" ;
    print OUT "Total\t$total\t100\n" ;
    print OUT "Full\t$full\t", (int($full/$total*10000+0.5))/100, "\n" ;
    print OUT "Part\t$part\t", (int($part/$total*10000+0.5))/100, "\n" ;
    print OUT "Null\t$null\t", (int($null/$total*10000+0.5))/100, "\n" ;
    print OUT "============\n" ;
    print OUT "UMI\t$umi_num\t", (int($umi_num/$total*10000+0.5))/100, "\n" ;
    print OUT "Final\t$final_num\t", (int($final_num/$total*10000+0.5))/100, "\n" ;
    close(OUT);
    open (OUT, ">$outdir/$prefix.qual.stat") || die "$outdir/$prefix.qual.stat, $!\n" ;
    print OUT "Number of Reads\t$total\t100\n" ;
    print OUT "Valid Barcodes\t$full\t", (int($full/$total*10000+0.5))/100, "\n" ;
    print OUT "Valid UMIs\t$umi_num\t", (int($umi_num/$total*10000+0.5))/100, "\n" ;
    print OUT "Final Valid Reads\t$final_num\t", (int($final_num/$total*10000+0.5))/100, "\n" ;
    close(OUT);
    open(OUT, ">$outdir/$prefix.bc_stat") || die "$outdir/$prefix.bc_stat, $!\n" ;
    print OUT "#Type\tReadsNum\tPercent(%)\n" ;
    for my $type (sort keys %hbc_stat){
        my $num = $hbc_stat{$type} ;
        my $per = (int($num/$total*10000+0.5))/100 ;
        print OUT "$type\t$num\t$per\n" ;
    }
    print OUT "======== single stat ========\n" ;
    for my $bc (sort keys %hsingle_stat){
        my $num = $hsingle_stat{$bc} ;
        my $per = (int($num/$total*10000+0.5))/100 ;
        print OUT "$bc\t$num\t$per\n" ;
    }
    close(OUT);

	return ;
}

#($bc1, $bc1_cigar, $bc1_pos, my $umi_pos, my $umi_seq) = &judge_bc_kmer('bc1', $seq, $ahbcs, $ahbc_seq, $ahbc_len, $ahref_pos, $ahkmer_pos, $ahkmer, $kmerlen, $ahumi_pos);
sub judge_bc_kmer()
{
    my ($type, $seq, $ahbcs, $ahbc_seq, $ahbc_len, $ahref_pos, $ahkmer_pos, $ahkmer, $kmerlen, $ahumi_pos) = @_ ;

    my ($start, $end) = @{$hbc_range{$type}} ;
    my $len = length($seq) ;
    my $offlen = $len > $end ? $end : $len ;
    my %htypes = ();
    my %hkmer_ref = ();
    for (my $i=$start; $i<$offlen-$kmerlen+1; $i++){
        my $kmer = substr($seq, $i, $kmerlen);
        if (defined $ahkmer->{$kmer}){
            for (my $j=0; $j<@{$ahkmer->{$kmer}}; $j++){
                next if ($ahkmer->{$kmer}[$j] !~ m/^$type/);
                my $ref_pos = $ahkmer_pos->{$kmer}[$j] ;
                my $map_pos = $i + 1;
                my $offset = abs($ref_pos - $map_pos) ;
                if ($offset > $max_offset){
                    next ;
                }
                if ($ahkmer->{$kmer}[$j] =~ m/(bc\d)\_(\d+)/){
                    $htypes{$1}{$ahkmer->{$kmer}[$j]}++ ;
                }
                push @{$hkmer_ref{$ahkmer->{$kmer}[$j]}}, [$ahkmer_pos->{$kmer}[$j], $map_pos] ;
            }
        }
    }
    for my $bc (sort {$htypes{$type}{$b} <=> $htypes{$type}{$a} || $a cmp $b} keys %{$htypes{$type}}){
        my $bc_seq = $ahbc_seq->{$bc} ;
        my $bc_len = $ahbc_len->{$bc} ;
        my @poses = @{$hkmer_ref{$bc}} ;
        my @bc_reg = ([$poses[0][0], $poses[0][0]]);
        my @read_reg = ([$poses[0][1], $poses[0][1]]);
        for (my $i=1; $i<@poses; $i++){
            my ($ref_pos, $read_pos) = @{$poses[$i]} ;
            if ($ref_pos - $bc_reg[-1][1] == 1 && $read_pos - $read_reg[-1][1] == 1){
                $bc_reg[-1][1] = $ref_pos ;
                $read_reg[-1][1] = $read_pos ;
            }
            else{
                push @bc_reg, [$ref_pos, $ref_pos] ;
                push @read_reg, [$read_pos, $read_pos] ;
            }
        }
        my ($flag, $match_num, $mismatch_num, $cigar, $bc_end_pos, $read_end_pos) = &judge_match(\@bc_reg, \@read_reg, $bc_seq, $bc_len, $seq, $kmerlen, $ahref_pos->{$type});
        if ($flag == 1){
            my $umi_pos = $ahumi_pos->{$type} + ($read_end_pos - $bc_end_pos) + $ahref_pos->{$type} ;
            my $umi_seq = substr($seq, $umi_pos-1, $umi_len);
            return($bc, $cigar, "$bc_end_pos,$read_end_pos", $umi_pos, $umi_seq);
            #$final_umi_pos = $umi_pos ;
            #$final_umi_seq = $umi_seq ;
            #push @end_poses, "$bc_end_pos,$read_end_pos" ;
            #push @types, $type ;
            #push @bcs, $bc ;
            #push @cigars, $cigar ;
            #$bc_flag = 1 ;
            #last ;
        }
    }

    return("null");
}

#my ($bc1, $bc1_cigar) = &judge_bc1('bc1', $bc1_seq, $ahbcs, $ahbc_seq, $ahkmer, $kmerlen, $ahbc_len);
sub judge_bc1()
{
    my ($type, $seq, $ahbcs, $ahbc_seq, $ahkmer, $kmerlen, $ahbc_len) = @_ ;

    my $bc_len = length($seq) ;
    if (defined $ahbcs->{$seq}){
        return($ahbcs->{$seq}, "$bc_len=") ;
    }
    my %hcandi_bc = ();
    for (my $i=0; $i<$bc_len-$kmerlen+1; $i+=$kmerlen){
        my $kmer = substr($seq, $i, $kmerlen) ;
        if (defined $ahkmer->{$kmer}){
            for my $bc (@{$ahkmer->{$kmer}}){
                next if ($bc !~ m/^$type/);
                $hcandi_bc{$bc} = $ahbc_seq->{$bc} ;
            }
        }
    }
    for my $bc (keys %hcandi_bc){
        my $bc_seq = $hcandi_bc{$bc} ;
        my $ref_len = $ahbc_len->{$bc} ;
        if ($bc_len == $ref_len){
            my ($mismatch, $pos) = &dist_count($bc_seq, $seq) ;
            if ($mismatch == 1){
                my $cigar = ($pos-1)."=1X".($bc_len-$pos)."=" ;
                if ($pos == 1){
                    $cigar = "1X".($bc_len-$pos)."=" ;
                }
                return($bc, $cigar) ;
            }
        }
        else{
            my $minlen = $bc_len > $ref_len ? $ref_len : $bc_len ;
            my $pre_match_num = 0 ;
            for (my $i=0; $i<$minlen; $i++){
                if (substr($seq, $i, 1) eq substr($bc_seq, $i, 1)){
                    $pre_match_num++ ;
                }
                else{
                    last ;
                }
            }
            if ($bc_len > $ref_len){
                if ($ref_len - $pre_match_num > 0){
                    my $rest_seq = substr($seq, $pre_match_num+1, $bc_len - $pre_match_num - 1) ;
                    my $rest_bc_seq = substr($bc_seq, $pre_match_num, $bc_len - $pre_match_num - 1) ;
                    if ($rest_seq eq $rest_bc_seq){
                        return ($bc, "$pre_match_num=1I".($bc_len-$pre_match_num-1)."=") ;
                    }
                }
                else{
                    return($bc, "$ref_len=1I");
                }
            }
            else{
                if ($bc_len - $pre_match_num > 0){
                    my $rest_seq = substr($seq, $pre_match_num, $bc_len - $pre_match_num) ;
                    my $rest_bc_seq = substr($bc_seq, $pre_match_num+1, $bc_len - $pre_match_num) ;
                    if ($rest_seq eq $rest_bc_seq){
                        return ($bc, "$pre_match_num=1D".($bc_len-$pre_match_num)."=") ;
                    }
                }
                else{
                    return($bc, "$bc_len=1D");
                }
            }
        }
    }

    return("null");
}

#my ($bc1, $bc1_pos) = &judge_bc('bc1', $bc1_seq, $ahbcs, $ahbc_seq, $ahkmer, $kmerlen);
sub judge_bc()
{
    my ($type, $seq, $ahbcs, $ahbc_seq, $ahkmer, $kmerlen) = @_ ;

    my $bc_len = length($seq) ;
    my $ref_len = $hbc_type_len{$type} ;
    if (defined $ahbcs->{$seq}){
        return($ahbcs->{$seq}, "$bc_len=") ;
    }
    if ($bc_len == $hbc_type_len{$type}){
        for (my $i=0; $i<$bc_len-$kmerlen+1; $i+=$kmerlen){
            my $kmer = substr($seq, $i, $kmerlen) ;
            if (defined $ahkmer->{$kmer}){
                for my $bc (@{$ahkmer->{$kmer}}){
                    next if ($bc !~ m/^$type/);
                    my $bc_seq = $ahbc_seq->{$bc} ;
                    my ($mismatch, $pos) = &dist_count($bc_seq, $seq) ;
                    if ($mismatch == 1){
                        my $cigar = ($pos-1)."=1X".($bc_len-$pos)."=" ;
                        if ($pos == 1){
                            $cigar = "1X".($bc_len-$pos)."=" ;
                        }
                        return($bc, $cigar) ;
                    }
                }
            }
        }
    }
    else{
        my %hcandi_bc = ();
        for (my $i=0; $i<$bc_len-$kmerlen+1; $i+=$kmerlen){
            my $kmer = substr($seq, $i, $kmerlen) ;
            if (defined $ahkmer->{$kmer}){
                for my $bc (@{$ahkmer->{$kmer}}){
                    next if ($bc !~ m/^$type/);
                    $hcandi_bc{$bc} = $ahbc_seq->{$bc} ;
                }
            }
        }
        my $minlen = $bc_len > $ref_len ? $ref_len : $bc_len ;
        for my $bc (keys %hcandi_bc){
            my $bc_seq = $hcandi_bc{$bc} ;
            my $pre_match_num = 0 ;
            for (my $i=0; $i<$minlen; $i++){
                if (substr($seq, $i, 1) eq substr($bc_seq, $i, 1)){
                    $pre_match_num++ ;
                }
                else{
                    last ;
                }
            }
            if ($bc_len > $ref_len){
                if ($ref_len - $pre_match_num > 0){
                    my $rest_seq = substr($seq, $pre_match_num+1, $bc_len - $pre_match_num - 1) ;
                    my $rest_bc_seq = substr($bc_seq, $pre_match_num, $bc_len - $pre_match_num - 1) ;
                    if ($rest_seq eq $rest_bc_seq){
                        return ($bc, "$pre_match_num=1I".($bc_len-$pre_match_num-1)."=") ;
                    }
                }
                else{
                    return($bc, "$ref_len=1I");
                }
            }
            else{
                if ($bc_len - $pre_match_num > 0){
                    my $rest_seq = substr($seq, $pre_match_num, $bc_len - $pre_match_num) ;
                    my $rest_bc_seq = substr($bc_seq, $pre_match_num+1, $bc_len - $pre_match_num) ;
                    if ($rest_seq eq $rest_bc_seq){
                        return ($bc, "$pre_match_num=1D".($bc_len-$pre_match_num)."=") ;
                    }
                }
                else{
                    return($bc, "$bc_len=1D");
                }
            }
        }
    }

    return("null");
}

#my ($mismatch, $pos) = &dist_count($bc_seq, $seq) ;
sub dist_count()
{
    my ($seq1, $seq2) = @_ ;
    my $mask = $seq1 ^ $seq2 ;
    my $mismatch = 0 ;
    my @poses = ();
    while($mask =~ m/[^\0]/g){
        $mismatch++ ;
        push @poses, pos($mask) ;
    }
    return($mismatch, @poses) ;
}
#my $flag = &check_umi_seq($final_umi_seq, $umi_len);
sub check_umi_seq()
{
    my ($seq, $umi_len) = @_ ;

    my $flag = 1 ;
    if ($seq =~ m/N|A{$umi_len}|C{$umi_len}|G{$umi_len}|T{$umi_len}/){
        $flag = 0 ;
    }

    return($flag);
}

#my ($flag, $match_num, $mismatch_num, $cigar, $bc_end_pos, $read_end_pos) = &judge_match(\@bc_reg, \@read_reg, $bc_seq, $bc_len, $seq, $kmerlen, $ahref_pos->{$type});
sub judge_match()
{
    my ($abc_reg, $aread_reg, $bc_seq, $bc_len, $seq, $kmerlen, $ref_start) = @_ ;

    for (my $i=0; $i<@{$abc_reg}; $i++){
        my ($bc_start, $bc_end) = @{$abc_reg->[$i]} ;
        my ($read_start, $read_end) = @{$aread_reg->[$i]} ;
        $bc_end += $kmerlen-1 ;
        $read_end += $kmerlen-1 ;
        my $match = $bc_end - $bc_start + 1 ;
        if ($match >= $bc_len-1){
            my $cigar = $bc_len."=" ;
            if ($bc_start > $ref_start){
                $cigar = "1X".($bc_len-1)."=" ;
            }
            elsif ($bc_end < $bc_len+$ref_start-1){
                $cigar = ($bc_len-1)."=1X" ;
            }
            return (1, $match, $bc_len - $match, $cigar, $bc_end, $read_end) ;
        }
        if ($bc_start > $ref_start && $bc_end < $bc_len+$ref_start-1){
            next ;
        }
        my ($mismatch, $cigar, $bc_end_pos, $read_end_pos) = &reg_match($bc_start, $bc_end, $read_start, $read_end, $bc_seq, $bc_len, $seq, $ref_start);
        if ($mismatch == 1){
            return (1, $bc_len-1, 1, $cigar, $bc_end_pos, $read_end_pos);
        }
    }

    return(0);
}

#my ($mismatch, $cigar, $bc_end_pos, $read_end_pos) = &reg_match($bc_start, $bc_end, $read_start, $read_end, $bc_seq, $bc_len, $seq, $ref_start);
sub reg_match()
{
    my ($bc_start, $bc_end, $read_start, $read_end, $bc_seq, $bc_len, $seq, $ref_start) = @_ ;

    my $mismatch = -1 ;
    my $cigar = "" ;
    my $bc_end_pos = $bc_end ;
    my $read_end_pos = $read_end ;
    if ($bc_start > $ref_start){
        my $bc_unmap_len = $bc_start - $ref_start ;
        my $bc_cut_sq1 = substr($bc_seq, 0, $bc_unmap_len);
        my $bc_cut_sq2 = substr($bc_seq, 0, $bc_unmap_len-1);
        my $read_cut_sq1 = substr($seq, $read_start-$bc_unmap_len-2, $bc_unmap_len);
        my $read_cut_sq2 = substr($seq, $read_start-$bc_unmap_len, $bc_unmap_len-1);
        my $read_cut_sq3 = substr($seq, $read_start-$bc_unmap_len-1, $bc_unmap_len-1);
#
#                  <----->   read_cut_sq1
#                    <---->  read_cut_sq2
#                   <---->   read_cut_sq3
# read          ------------------------------
# match position           ||||||
# barcode           -------------
#                   <----->  bc_cut_sq1
#                   <---->   bc_cut_sq2
#
        if ($bc_cut_sq2 eq $read_cut_sq2){
            $mismatch = 1 ;
            $cigar = ($bc_unmap_len-1)."=1D".($bc_end-$bc_start+1)."=" ;
        }
        elsif ($bc_cut_sq2 eq $read_cut_sq3){
            $mismatch = 1 ;
            $cigar = ($bc_unmap_len-1)."=1X".($bc_end-$bc_start+1)."=" ;
        }
        elsif ($bc_cut_sq1 eq $read_cut_sq1){
            $mismatch = 1 ;
            $cigar = ($bc_unmap_len)."=1I".($bc_end-$bc_start+1)."=" ;
        }
        else{
            return(-1);
        }
    }
    elsif ($bc_end < $ref_start+$bc_len-1){
        my $bc_pos_end = $bc_end - $ref_start + 1 ;
        my $bc_cut_sq1 = substr($bc_seq, $bc_pos_end, $bc_len-$bc_pos_end);
        my $bc_cut_sq2 = substr($bc_seq, $bc_pos_end+1, $bc_len-$bc_pos_end-1);
        my $read_cut_sq1 = substr($seq, $read_end+1, $bc_len-$bc_pos_end-1);
        my $read_cut_sq2 = substr($seq, $read_end, $bc_len-$bc_pos_end-1);
        my $read_cut_sq3 = substr($seq, $read_end+1, $bc_len-$bc_pos_end);
#
#                          <---->  read_cut_sq1
#                         <---->   read_cut_sq2
#                          <-----> read_cut_sq3
# read          ------------------------------
# match pos         ||||||   
# barcode           -------------
#                         <----->  bc_cut_sq1
#                          <---->  bc_cut_sq2
#
        if ($bc_cut_sq2 eq $read_cut_sq1){
            $mismatch = 1 ;
            $cigar = "$bc_pos_end=1X".($bc_len-$bc_pos_end-1)."=" ;
            $bc_end_pos = $ref_start + $bc_len - 1 ;
            $read_end_pos = $read_end + $bc_len - $bc_pos_end ;
        }
        elsif ($bc_cut_sq1 eq $read_cut_sq3){
            $mismatch = 1 ;
            $cigar = "$bc_pos_end=1I".($bc_len-$bc_pos_end)."=" ;
            $bc_end_pos = $ref_start + $bc_len - 1 ;
            $read_end_pos = $read_end + $bc_len - $bc_pos_end + 1 ;
        }
        elsif ($bc_cut_sq2 eq $read_cut_sq2){
            $mismatch = 1 ;
            $cigar = "$bc_pos_end=1D".($bc_len-$bc_pos_end-1)."=" ;
            $bc_end_pos = $ref_start + $bc_len - 1 ;
            $read_end_pos = $read_end + $bc_len - $bc_pos_end - 1 ;
        }
        else{
            return(-1);
        }
    }

    return($mismatch, $cigar, $bc_end_pos, $read_end_pos);
}

#&output_bc_umi(\%hbc_umi, $outdir, $prefix);
sub output_bc_umi()
{
    my ($ahbc_umi, $outdir, $prefix) = @_ ;

    open (OUT, ">$outdir/$prefix.bc_umi") || die "$outdir/$prefix.bc_umi, $!\n" ;
    for my $bc (keys %{$ahbc_umi}){
        for my $umi (keys %{$ahbc_umi->{$bc}}){
            my $num = scalar (@{$ahbc_umi->{$bc}{$umi}});
            print OUT "$bc\t$umi\t$num" ;
            for my $id (@{$ahbc_umi->{$bc}{$umi}}){
                print OUT "\t$id" ;
            }
            print OUT "\n" ;
        }
    }
    close(OUT);

    return ;
}

