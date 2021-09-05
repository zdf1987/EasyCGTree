#!/usr/bin/perl
use warnings;
use strict;
use File::Copy qw(copy mv);
use Getopt::Long;
my @usage=qq(
====== EasyCGTree ======
     Version 3.1 by Dao-Feng Zhang (Hohai University, China)
	 Update 2021-06-22

Usage:\n
	perl GetRepRef.pl [Options]

Options:
-input <String>		
	Input data directory (essential option)
-gene_number <Int>		
	Number of genes randomly selected for distance calculation (6-24, default: 12; smaller number, less time cost)
-diverg_cutoff <Decimal, 0-0.5>	
	Minimum distance allowed to screen representative data sets. [default: 0.05]
-thread	<Int>		
	Number of threads to be used by 'clustalo'. [default: 2]
-local_dist				
	Distance matrix allready exists (no need to specify in the command line; this option can save time of doing alignment, and the file needs to be named after the input directory, for example: input.fas.dist).

-help					
	Display this message


);
my $inputDir;
my $genenum =12;
my $thread=2;
my $cutoff=0.05;


my %opt=qw();
GetOptions(\%opt,"input:s","gene_number:i","diverg_cutoff:f","thread:i","local_dist!","help|h!");
if (scalar(keys %opt )==0 || exists($opt{"help"})) {
	print join("\n",@usage)."\n\n";
	exit;
}
if (exists($opt{"input"})) {
	$inputDir = $opt{"input"};
} else {
	print join("\n",@usage)."\n\n";
	exit;	
}

$genenum = $opt{"gene_number"} if exists($opt{"gene_number"});
$thread=$opt{"thread"} if exists($opt{"thread"});
$cutoff=$opt{"diverg_cutoff"} if exists($opt{"diverg_cutoff"});
my $local_dist=1 if exists($opt{"local_dist"});

open (OUT, ">$inputDir.fas");
	
unless ($local_dist) {
	opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
	my @files = readdir(DIR);
	closedir(DIR);
	die "Input directory $inputDir does not contain any files" unless scalar(@files);
	my $time=0;
	my @geneloc;
	my @gene;	
	while (1) {
		my $nn =int(rand(116))-1;
		push @geneloc, $nn unless $nn ~~ @geneloc;
		last if $#geneloc >= $genenum;
	}
	@geneloc=sort @geneloc;
	my @outset;
	foreach (@files) {
		next unless /\.fas\Z/;
		my $seq;
		open (FILE, "./$inputDir/$_") or die "Can't open '$_': $!";
		unless (@gene) {
			my $count=0;
			my $sig=0;
			my @gene0;
			foreach my $in (<FILE>) {
				$in=~ s/\n// unless $in=~s/\r\n//;
				if ($in =~ /\>/) {
					if ($count ~~ @geneloc) {
						my @ii=split /_/, $in;
						push @gene0, $ii[3];
						$count++;
						$sig=1;
					} else {
						$count++;
						$sig=0;
					}
				} elsif ($sig) {

					unless ($seq) {
						$seq="$in";
						} else {
							$seq .="$in";
						}
				}
			}
			close FILE;
			if ($count >114) {
				@gene=@gene0;
				my $out=&formatSeq2fas($seq);
				print OUT ">$_\n$out";
			} else {
				push @outset, $_;
			}
		} else {
			my $count=0;
			my $sig=0;
			foreach my $in (<FILE>) {
				$in=~ s/\n// unless $in=~s/\r\n//;
				if ($in =~ /\>/) {
					my @ii=split /_/, $in;						
					$count++;
					if ($ii[3] ~~ @gene) {	
						$sig=1;
					} else {
						$sig=0;
					}
				} elsif ($sig) {
					unless ($seq) {
						$seq="$in";
						} else {
							$seq .="$in";
						}
				}
			}
			close FILE;
			if ($count >114) {
				my $out=&formatSeq2fas($seq);
				print OUT ">$_\n$out";
			} else {
				push @outset, $_;
			}
		}
	} 	
		
	close OUT;
	print "\n\n#######################################################\nThese $genenum genes were randomly selected for distance caculation: @gene\n###########################################################\n\n";
	my @para = ("-i", "./$inputDir.fas", "-o", "./$inputDir.aln.fas", "--outfmt=fasta","--distmat-out=./$inputDir.fas.dist", "--full","--output-order=tree-order", "--force","--threads=$thread","-v");
	system ('./bin/clustalo', @para);

	print "\n\n#######################################################\nThe alignment and distance caculation have been finished. These data sets were excluded from the following analysis for their gene number <115: @outset\n###########################################################\n\n";
} 
print "I am searching the most divergent taxa......\n";

open (FILE, "./$inputDir.fas.dist") or die "Can't open '$inputDir.fas.dist': $!";
my @taxa;
while (<FILE>) {
    chomp;
    my @in = split /\s+/;
	next unless $in[1];
     push @taxa, $in[0];
}
close FILE;
my (@taxa1,@taxa2,@dis);
my $line=0; 
open (FILE, "./$inputDir.fas.dist") or die "Can't open '$inputDir.fas.dist': $!";
while (<FILE>) {
     chomp;
     $line++;
     next if $line < 3;
     my @in = split /\s+/;
     my $end = $line-2;
     foreach (1..$end) {
        push @taxa1, $taxa[$end]; 
		push @taxa2, $taxa[$_-1]; 
		push @dis, $in[$_]; 
     }
  }
close FILE;

my @outRef;
while (1) {
	my $maxid;
	my $maxvalue=0;
	my @maxloc;
	foreach my $id (@taxa) {
		my $temid;
		my $temvalue=0;
		my @temloc;
		foreach my $num (0..$#taxa1) {###print "$taxa1[$num] eq $id || $taxa2[$num] eq $id \n";
			next unless $taxa1[$num] eq $id || $taxa2[$num] eq $id;
			$temvalue += $dis[$num];
			push @temloc, $num;
		}
		if ($maxvalue < $temvalue) {
			$maxvalue = $temvalue ;
			$maxid = $id;
			@maxloc=@temloc;
		}
	}
	push @outRef, $maxid;
	my $aver=$maxvalue/($#maxloc+1);
	print "$maxid was selected. The average distance with the rest taxa is: $aver.\n";
exit unless $maxid;
	my @outid;
	push @outid, $maxid;
	foreach my $num1 (@maxloc) {
		if ($dis[$num1]<= $cutoff) {
			push @outid, $taxa1[$num1] unless $taxa1[$num1] ~~ @outid;
			push @outid, $taxa2[$num1] unless $taxa2[$num1] ~~ @outid;
		}
	}
	my @outloc1;
	foreach my $num2 (0..$#taxa) {
		push @outloc1, $num2 if $taxa[$num2]~~@outid;
	}
	@outloc1=reverse @outloc1;

	foreach my $num21 (@outloc1) {
		splice @taxa, $num21, 1;
	}

	my @outloc2;
	foreach my $num22 (0..$#taxa1) {
		push @outloc2, $num22 if $taxa1[$num22]~~@outid ||$taxa2[$num22]~~@outid;
	}
	@outloc2=reverse @outloc2;
	foreach my $num3 (@outloc2) {
		splice @taxa1, $num3, 1;
		splice @taxa2, $num3, 1;
		splice @dis, $num3, 1;
	}
	last unless $taxa[1];
}
push @outRef, $taxa[0] if $taxa[0];
mkdir "Rep_$inputDir";	
foreach my $kk (@outRef) {
	copy ("./$inputDir/$kk","./Rep_$inputDir/$kk")|| die "Step 3: Can't move './$inputDir/$kk' into './Rep_$inputDir/$kk'!\n";
}	
my $number =$#outRef+1;
print "\nThese $number taxa were selected as presentatives: @outRef\n\n###########################################################\n\nA copy of related files had been written into directory: ./Rep_$inputDir/\n\n";

sub formatSeq2fas {
	my $i1=shift;
	$i1=~s/\*//g;
	my $subout;
	my $sub;
	my $start=0;
	while (1) {
		if (length($i1)-$start >60) {
			$sub = substr($i1,$start,60);
			unless ($subout) {
				$subout="$sub\n";
			} else {
				$subout .="$sub\n";
			}
		} else {
			unless ($subout) {
				$subout="$sub\n";
			} else {
				$subout .="$sub\n";
			}
			return $subout;
			last;
		}
		$start +=60;
	}
}
