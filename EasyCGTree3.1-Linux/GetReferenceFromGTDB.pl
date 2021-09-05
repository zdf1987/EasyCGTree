#！c:/Perl64/bin/perl
use strict;
use warnings;####
use List::MoreUtils qw/uniq/;
my @usage=qq(
====== EasyCGTree ======
     Version 3.1 by Dao-Feng Zhang (Hohai University, China)
	 Update 2021-06-22

Usage:\n perl_script + bac120_taxonomy_r95.tsv SeqDIR ID

);

die join("\n",@usage)."\n" unless $ARGV[2];  
print join("\n",@usage)."\n";
my $inputDir = $ARGV[1];
opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
my @query = readdir(DIR);
closedir(DIR);
die "Input directory $inputDir does not contain any files" unless scalar(@query);

my @inlist;
my $rrr=0;
foreach my $file (@query) {
	
	next unless $file =~ /\.faa/;
	$rrr++;
	print "Reading $inputDir/$file ($rrr of 120(122) files)...\n";#####
	

	open (FIL, "<$inputDir/$file") or die "Can't open '$file': $!";
	foreach my $in (<FIL>) {
		$in=~ s/\n// unless $in=~s/\r\n//;
		next unless $in;
		if ($in=~/\>/) {
			$in=~ s/\>//;
			push @inlist, $in;
		}
	}
	
}
@inlist= uniq @inlist;


my @list;
my $innu=0;
my $lacknum=0;
open (FILE, "<$ARGV[0]") or die "Can't open '$ARGV[0]': $!";
my (%out,%name,%spe,%num,%lnum,%lack,%list);

print "\n\nReadling $ARGV[0] ...";
foreach my $in (<FILE>) {
	$in=~ s/\n// unless $in=~s/\r\n//;


	next unless $in =~ /$ARGV[2]/i;
	$in=~ s/ /-/g;
	$innu++;
	
	my @in =split /\t/, $in;
	
	my @in1 =split /__/, $in;
	$lacknum++ unless $in[0] ~~ @inlist;
	push @list, $in[0] if $in[0] ~~ @inlist;
	$spe{$in[0]}="$in1[7],$in[0],$in[1]";####
	$name{$in[0]}="$in1[7].fas";
}
close FILE;
print "\n\nWARNNING: There are $innu items collected in the $ARGV[0] and belonging to $ARGV[2].\n I find no sequence for $lacknum items.\n\n";



my @genelist;
my $head = "Species,Assemble no.,Taxonomy,Present,Absent,Absent list,";
my $filenum=0;
foreach my $file (@query) {
	
	next unless $file =~ s/\.faa//;
	$filenum++;
	
	print "Handling $inputDir/$file ($filenum of 120/122 files)...\n";
	
	$head .="$file,";
	push @genelist, $file;
	open (FIL, "<$inputDir/$file.faa") or die "Can't open '$file': $!";
    my ($id,$sig,$sig0);
		my $num=0;
	my @list0;
	
	foreach my $in (<FIL>) {
		$in=~ s/\n// unless $in=~s/\r\n//;

		next unless $in;
		if ($in=~/\>/) {
			
			$in=~ s/\>//;
			if ($in~~ @list) {
				push @list0, $in;
				unless ($out{$in}) {
					$out{$in} =">$in\_$file\n";
					$num{$in} =1;
					if ($list{$in}) {
						$list{$in} .="+,";
					} else {
						$list{$in} ="+,";
					}
				} else {
					$out{$in} .=">$in\_$file\n";
					$num{$in}++;
					$list{$in} .="+,";
				} 
				
				$sig=1;
				$id =$in;
				$num++;
			} else {
				$sig=0;
	
			}
			last if $num == $#list+1;
		} elsif ($sig==1) {
			$out{$id} .="$in\n";
			
		}
		
	}
	close FIL; 
	foreach my $i (@list) {
		next if $i ~~ @list0;
		$num{$i} =0 unless $num{$i};
		
		if ($list{$i}){
			$list{$i} .="--,";
		} else {
			$list{$i} ="--,";
		}
		if ($lack{$i}) {
			$lack{$i} .="$file ";###
		} else {
			$lack{$i} ="$file ";####
		}
		$lnum{$i}++ if $lnum{$i};
		$lnum{$i} =1 unless $lnum{$i};
		
	}
	
	
}
open(OU, ">$ARGV[2]_log.csv") or die "Can't open '$ARGV[2]_log.csv': $!";
print OU "$head\n";
mkdir "$ARGV[2]" || die "Can't create blastDB directory '$ARGV[2]'\n";


my $innum=0;
foreach my $ii (keys %num) {
	$num{$ii} =0 unless $num{$ii};
	$lnum{$ii} =0 unless $lnum{$ii};
	$lack{$ii} ="--" unless $lack{$ii};#####
	next unless $out{$ii};
	
	print "I found $num{$ii} sequence totally for $name{$ii}!!\n\n";
	print OU "$spe{$ii},$num{$ii},$lnum{$ii},$lack{$ii},$list{$ii}\n";	
	$innum++;
	
	open(OUT, ">$ARGV[2]/$name{$ii}");
	print OUT "$out{$ii}";
	close OUT;	

}
close OU;
print "Job was finished: \n$innum gene sets were extracted, and please find the sequence in folder $ARGV[2] and more information in file $ARGV[2]_log.csv.\n";




