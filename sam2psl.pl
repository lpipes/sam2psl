#!/usr/bin/perl -w
#-----------------------------------------------------------------------------------------------------
#Author: Lenore Pipes, lp364@cornell.edu
#Copyright 2018

use strict;
use warnings;
if (scalar(@ARGV) < 2){
	print "sam2psl.pl - Convert sam file into psl format. Takes as input the target's chrom sizes and outputs alignment.psl file\n";
	print "usage:\n";
	print "sam2psl.pl target.chrom.sizes alignment.sam\n";
	print "\n";
	exit;
}
my $sam_file = $ARGV[1];
my $chrom_file = $ARGV[0];
my @splitname = split(/\./,$sam_file);
my $outfile = $splitname[0] . '.psl';
my %chrom_sizes = ();
open (OUTFILE, ">$outfile") || die ("Could not open $outfile!");
print OUTFILE "psLayout version 3\n";
print OUTFILE "\n";
print OUTFILE "match\tmis-\trep.\tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ\tQ\tQ\tQ\tT\tT\tT\tT\tblock\tblockSizes\tqStarts\ttStarts\n";
print OUTFILE "\tmatch\tmatch\t\tcount\tbases\tcount\tbases\t\tname\tsize\tstart\tend\tname\tsize\tstart\tend\tcount\n";
print OUTFILE "---------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
open (CHROMS, $chrom_file) || die ("Could not open chromosome sizes file!");
while (<CHROMS>){
	my $line = $_;
	chomp($line);
	my @spl = split('\t',$line);
	my $chrom_name = $spl[0];
	my $chrom_size = $spl[1];
	$chrom_sizes{$chrom_name} = $chrom_size;
}
close(CHROMS);
open (SAM, $sam_file) || die ("Could not open sam file!");
while (<SAM>){
        my $line = $_;
        chomp($line);
	if ($line !~ m/^\@/){
		my @spl = split('\t',$line);
		my $chrom = $spl[2];
		my $tSize = $chrom_sizes{$chrom};
		#if there is a match, process line
		if ($chrom ne "*"){
			my $qName = $spl[0];
			my $flag = $spl[1];
			my $tStart = $spl[3] - 1;
			my $cigar = $spl[5];
			my $sequence = $spl[9];
			my $misMatches = $spl[14];
			my @spl5 = split("nM:i:",$misMatches);
			$misMatches = $spl5[1];
			my $qSize = seqLength($sequence);
			my $totalmatches = 0;
			my $blockCount = 0;
			#process cigar string
			if ( $cigar =~ m/M/ ){
				my @counts = getCount($cigar,"M");
				$totalmatches = $counts[0];
				$blockCount = $counts[1];
			}
			my $nCount = 0;
			if ($cigar =~ m/N/ ){
				my @counts = getCount($cigar,"N");
				$nCount = $counts[0];
			}
			my $tBaseInsert = 0;
			my $tNumInsert = 0;
			if ($cigar =~ m/I/ ){
				my @counts = getCount($cigar,"I");
				$tBaseInsert = $counts[0];
				$tNumInsert = $counts[1];
			}
			my $qBaseInsert = 0;
			my $qNumInsert = 0;
			if ($cigar =~ m/D/ ){
				my @counts = getCount($cigar,"D");
				$qBaseInsert = $counts[0];
				$qNumInsert = $counts[1];
			}
			my $qStart = 0;
			if ($cigar =~ m/^[0-9]+S/ ){
				my @counts = getCount($cigar,"S");
				$qStart = $counts[2];
			}elsif( $cigar =~ m/^[0-9]+H/ ){
				my @counts = getCount($cigar,"H");
				$qStart = $counts[2];
			}
			my $qEnd = 0;
			my $matches = $totalmatches - $misMatches;
			my @qStart_blocks;
			$qStart_blocks[$blockCount] = undef;
			$qStart_blocks[0] = $qStart;
			my @blockSizes;
			my @qStarts;
			my $qStart_sum=$qStart;
			my @tStarts;
			my $tStart_sum = $tStart;
			if ($cigar =~ m/M/){
				my $tmp = $cigar;
				$tmp =~ s/[DNM=XISHP]/|/g;
				my @spl = split(/\Q|/,$tmp);	
				my $tmp2 = $cigar;
				$tmp2 =~ s/[0-9]+/|/g;
				my @spl2 = split(/\Q|/,$tmp2);
				for(my $i=0; $i<scalar(@spl2); $i++){
					if ($spl2[$i] eq 'M' or $spl2[$i] eq 'I'){
						if ($spl2[$i] eq 'M'){
							push(@blockSizes,$spl[$i-1]);
							my $scale = scalar(@blockSizes);
							if ( scalar(@blockSizes) > 1 ){
								push(@qStarts,$qStart_sum);
							}
							else{ push(@qStarts, $qStart); }
							$qStart_sum = $qStart_sum + $spl[$i-1];
						}else{
							$qStart_sum = $qStart_sum + $spl[$i-1];
						}	
					}
					if ($spl2[$i] eq 'M' or $spl2[$i] eq 'D' or $spl2[$i] eq 'N'){
						if ($spl2[$i] eq 'M'){
							if ( scalar(@tStarts) > 0 ){
								push(@tStarts,$tStart_sum);
							}else{ 
								push(@tStarts,$tStart); 
							}
							$tStart_sum = $tStart_sum + $spl[$i-1];
						}else{
							$tStart_sum = $tStart_sum + $spl[$i-1];
						}

	

					}	
				}
			
			}
			$qEnd = $qStarts[scalar(@qStarts)-1]+$blockSizes[scalar(@qStarts)-1];
			if ($cigar =~ m/[0-9]+H$/ ){
                                my @counts = getCount($cigar,"H");
                                $qEnd = $qSize - $counts[3];
                        }elsif ($cigar =~ m/[0-9]+S$/ ){
                                my @counts = getCount($cigar, "S");
                                $qEnd = $qSize - $counts[3];
                        }
			$cigar =~ s/([0-9]+[ISHP])+//g;
			my @spl2 = split('[MDN=X]',$cigar);
			my $tAlignment=0;
			for(my $i=0; $i<scalar(@spl2); $i++){
				$tAlignment = $tAlignment + $spl2[$i];
			}
			my $tEnd = $tStart + $tAlignment - 2;
			my @spl3 = split('',$flag);
			my $binary = dec2bin($flag);
			my @spl4 = split('',$binary);
			my $strand_code = $spl4[7];
			my $strand;
			if ($strand_code == 1){
				$strand = "-";
			}elsif($strand_code == 0){
				$strand = "+";
			}
			my $repMatches = 0;
			my $blockSizes_print = join(",",@blockSizes);
			my $qStarts_print = join(",",@qStarts);
			my $tStarts_print = join(",",@tStarts);
			print OUTFILE "$matches\t$misMatches\t$repMatches\t$nCount\t$qNumInsert\t$qBaseInsert\t$tNumInsert\t$tBaseInsert\t$strand\t$qName\t$qSize\t$qStart\t$qEnd\t$chrom\t$tSize\t$tStart\t$tEnd\t$blockCount\t$blockSizes_print,\t$qStarts_print,\t$tStarts_print,\n";
		}
	}
}
close(SAM);
#converts to binary
sub dec2bin {
	my $str = unpack("B32", pack("N", shift));
	my @spl = split('',$str);
	my $binary = $spl[20];
	for (my $i=21; $i<=31; $i++){
		$binary .= $spl[$i];
	}	
	return $binary;
}
#get length of sequence
sub seqLength {
	my @spl = split('',$_[0]);
	return scalar(@spl);
}
sub getCount {
	my $cigar = $_[0];
	my $symbol = $_[1];
	my @spl;
	my $beg=0;
	my $end=0;
	if ($symbol eq 'M'){
		$cigar =~ s/([0-9]+[ISHPDN=X])+//g;
		@spl = split('M', $cigar);
	}
	if ($symbol eq 'N'){
		$cigar =~ s/([0-9]+[ISHPDM=X])+//g;
		@spl = split('N',$cigar);
	}
	if ($symbol eq 'I'){
		$cigar =~ s/([0-9]+[NSHPDM=X])+//g;
		@spl = split('I',$cigar);
	}
	if ($symbol eq 'D'){
		$cigar =~ s/([0-9]+[NSHPIM=X])+//g;
		@spl = split('D',$cigar);
	}
	if ($symbol eq 'S'){
		if ( $cigar =~ m/^[0-9]+S/){
		     	$cigar =~ s/([0-9]+[NDHPIM=X])+//g;
			@spl = split('S',$cigar);
			if ( scalar(@spl) > 1){
                        	$beg = $spl[0];
				$end = $spl[1];
			}else{
				$beg = $spl[0];
			}
                }else{
			$cigar =~ s/([0-9]+[NDHPIM=X])+//g;
			@spl = split('S',$cigar);
			$end = $spl[0];
		}
	}
	if ($symbol eq 'H'){
		if ( $cigar =~ m/^[0-9]+H/){
			$cigar =~ s/([0-9]+[NDSPIM=X])+//g;
			@spl = split('H',$cigar);
			if (scalar(@spl) > 1){
				$beg = $spl[0];
				$end = $spl[1];
			}else{
				$beg=$spl[0];
			}
		}else{
			$cigar =~ s/([0-9]+[NDSPIM=X])+//g;
			@spl = split('H',$cigar);
			$end = $spl[0];
		}
	}
	my $count=0;
	for(my $i=0; $i<scalar(@spl); $i++){
		$count = $count + $spl[$i];
	}
	return $count,scalar(@spl),$beg,$end;
}
