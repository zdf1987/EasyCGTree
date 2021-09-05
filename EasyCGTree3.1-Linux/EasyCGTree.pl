#!/usr/bin/perl
use warnings;
use strict; 
no warnings 'experimental::smartmatch';
use File::Copy qw(move mv);
use Getopt::Long;
use List::MoreUtils qw/uniq/;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $rmon=$mon +1;
my $ryear=$year+1900;

print "\n##############################\nReading command line... checking options and input data...\n\n";
my @usage=qq(
====== EasyCGTree ======
     Version 3.1 by Dao-Feng Zhang (Hohai University, China)
	 Update 2021-08-02

Usage: perl EasyCGTree.pl [Options]

Options:
-input <String>	(Essential)
	Input data (genome/proteome) directory
-query <String>	(Essential)
	Directory of protein sequence data for gene calling
-seq_type <String, 'nucl', 'prot'> (Essential)
	Sequence type used for tree inference
	
-mode <String, 'A', 'A1', 'A2', 'A3'> (Optional)
	Set run mode. [default: A]
-tree <String, 'sm', 'st', 'snp'>(Optional)				
	Approach (sm, supermatrices; st, supertree; or snp, SNP) used for tree inference. [default: sm]
-reference <String> (Optional)
	Directory of protein sequence data used in tree inference
-thread <Int> (Optional)		
	Number of threads to be used by 'blast+'. [default: 2]
-blast_dir <String> (Optional)
	Directory for programs of blast+. [default: ./bin/]
-iden_cutoff <Int, 50..100> (Optional)
	Cutoff(%) for filtering the BLAST results. [default: 50]
-gene_cutoff <Decimal, 0.5..1> (Optional)
	Cutoff for omitting low-prevalence gene. [default: 0.8]
-genome_cutoff <Decimal, 0.5..1> (Optional)
	Cutoff for omitting low-quality genomes. [default: 0.8]
	
-help (Optional)				
	Display this message

);
my ($inputDir,$queryDir,$seqtype,$Reference,$blast,$check1,$check2,$check3,@query,@queryR);


########## Default parameters ################
my $tblastn="./bin/";
my $blastIdentitycutoff=50;
my $geneCutoff=0.8;
my $genomeCutoff=0.8;
my $num_threads=2;
my $tree_type ="sm";
my $mode = "A";
##########################################

my %opt=qw();
GetOptions(\%opt,"input:s","query:s","mode:s","seq_type:s","tree:s","reference:s","thread:i","blast_dir:s","iden_cutoff:f","gene_cutoff:f","genome_cutoff:f","help!");

if (scalar(keys %opt )==0 || exists($opt{"help"})) {
	print join("\n",@usage)."\n\n";
	exit;
}

unless (exists($opt{"seq_type"})+exists($opt{"query"})+exists($opt{"input"})==3) {
	print join("\n",@usage)."\n\nERROR: '-input', '-query' and '-seq_type' must be set to start.\n\n";
	exit;	
}	

{
	$inputDir = $opt{"input"};
	unless ($inputDir=~/\w/ || $inputDir=~/\d/) {
		print join("\n",@usage)."\n\nERROR: Argument '-input $inputDir': it should be a directory name, e.g., myGenome.\n\n";
		exit;
	}
}

{
	$queryDir = $opt{"query"};
	$check1= &checkDataType($queryDir);
	if ($check1 eq "nucl") {
		print join("\n",@usage)."\n\nERROR: Argument '-query $queryDir': directory '$queryDir' should contains file(s) of protein sequence, while only file(s) of nucleotide sequence were found.\n\n";
		exit;
	} elsif ($check1 eq "prot") {
	} elsif ($check1=~ /\.fas/ || $check1=~ /open input/) {
		print join("\n",@usage)."\n\nERROR: Argument '-query $queryDir': $check1\n\n";
		exit;
	} else {
		print join("\n",@usage)."\n\nERROR: Argument '-query $queryDir': directory '$queryDir' should only contains fasta-formated file(s) of protein sequence, while $check1\n\n";
		exit;		
	}
}

if (exists($opt{"reference"})) {
	$Reference = $opt{"reference"};
	$check3= &checkDataType($Reference);
	if ($check3 eq "nucl") {
		print join("\n",@usage)."\n\nERROR: Argument '-reference $Reference': directory '$Reference' should contains file(s) of protein sequence.\n\n";
		exit;
	} elsif ($check3 eq "prot") {
	} elsif ($check3=~ /\.fas/ || $check3=~ /open input/) {
		print join("\n",@usage)."\n\nERROR: Argument '-reference $Reference': $check3\n\n";
		exit;
	} else {
		print join("\n",@usage)."\n\nERROR: Argument '-reference $Reference': directory '$Reference' should contains fasta-formated file(s) of protein sequence, while $check1\n\n";
		exit;		
	}
	opendir(DIR, $Reference);
	@queryR = readdir(DIR);
	closedir(DIR);
	splice @queryR, 0, 1 unless $queryR[0]=~ /\.fas/;
	splice @queryR, 1, 1 unless $queryR[1]=~ /\.fas/;
}

$seqtype = $opt{"seq_type"};
if (exists($opt{"mode"})) {
	$mode = $opt{"mode"};
	if ($mode eq "A1") {
		if ($Reference) {
			print join("\n",@usage)."\n\nERROR: '-mode A/A1/A2' dose not allow the setting of '-reference '.\n\n";
			exit;
		}
	} elsif ($mode eq "A2") {
		if ($Reference) {
			print join("\n",@usage)."\n\nERROR: '-mode A/A1/A2' dose not allow the setting of '-reference '.\n\n";
			exit;
		}
	} elsif ($mode eq "A3") {
		
		unless ($Reference) {
			print join("\n",@usage)."\n\nERROR: '-mode A3' needs the setting of '-reference '.\n\n";
			exit;
		}
		unless  ($seqtype eq "prot") {
			print join("\n",@usage)."\n\nERROR: '-mode A3' needs the setting of '-seq_type prot'.\n\n" ; 
			exit;
		}
		
	} else {
		unless ($mode eq "A") {
			print LOG "\n\nERROR: Argument '-mode $mode': it should be 'A', 'A1', 'A2' or 'A3'.\n\n";
			exit;
		}
		
		if ($Reference) {
			print join("\n",@usage)."\n\nERROR: '-mode A/A1/A2' dose not allow the setting of '-reference '.\n\n";
			exit;
		} 
	}
} else {
		if ($Reference) {
		print join("\n",@usage)."\n\nERROR: Argument '-reference $Reference': no setting of '-mode' in command line means setting of '-mode A', which dose not allow the setting of '-reference '.\n\n";
		exit;
	} 
}

{
	$check2= &checkDataType($inputDir);
	if ($seqtype eq "nucl") {
		$blast = "tblastn";
		if ($check2 eq "nucl") {
		} elsif ($check2 eq "prot") {
			print join("\n",@usage)."\n\nERROR: '-seq_type nucl' only allows input directory '$inputDir' containing file(s) of nucleotide sequence.\n\n";
			exit;
		} elsif ($check2=~ /\.fas/ || $check2=~ /open input/) {
			print join("\n",@usage)."\n\nERROR: Argument '-input $inputDir': $check2\n\n";
			exit;
		} else {
			print join("\n",@usage)."\n\nERROR: '-seq_type nucl' only allows input directory '$inputDir' containing fasta-formated file(s) of nucleotide sequence, while $check2\n\n";
			exit;	
		}
	} elsif ($seqtype eq "prot") {
		$blast = "blastp";
		if ($check2 eq "nucl") {
			print join("\n",@usage)."\n\nERROR: '-seq_type prot' only allows input directory '$inputDir' containing file(s) of protein sequence.\n\n";
			exit;
		} elsif ($check2 eq "prot") {
		} elsif ($check2=~ /\.fas/|| $check2=~ /open input/) {
			print join("\n",@usage)."\n\nERROR: Argument '-input $inputDir': $check2\n\n";
			exit;
		} else {
			print join("\n",@usage)."\n\nERROR: '-seq_type prot' only allows input directory '$inputDir' containing fasta-formated file(s) of protein sequence, while $check2\n\n";
			exit;		
		}
	} else {
		print join("\n",@usage)."\n\nERROR: Argument '-seq_type $seqtype': it should be 'nucl' or 'prot'.\n\n";
		exit;
	}
	opendir(DIR, $inputDir);
	@query = readdir(DIR);
	closedir(DIR);
	splice @query, 1, 1 unless $query[1]=~ /\.fas/;
	splice @query, 0, 1 unless $query[0]=~ /\.fas/;
	my $fn1=$#query+1;
	my $fn2=1 + $#queryR;
	unless ($fn1 + $fn2 >=4) {
		if ($Reference) {
			print join("\n",@usage)."\n\nERROR: Number of files (taxa on a tree) in both Input directory '$inputDir' ($fn1) and Reference directory '$Reference' ($fn2) must be >=4 to start a run.\n\n";
		} else {
			print join("\n",@usage)."\n\nERROR: Number of files (taxa on a tree) in Input directory '$inputDir' ($fn1) must be >=4 to start a run.\n\n";
		}
		exit;
	}
}

open (LOG, ">$inputDir\_$ryear-$rmon-$mday-$hour-$min.log") or die "Can't open '$inputDir\_$ryear-$rmon-$mday-$hour-$min.log': $!\n";

my $treesig;
if (exists($opt{"tree"})) {
	$tree_type = $opt{"tree"};
	if ($tree_type eq "st") {
		unless ($seqtype eq "prot") {
			print join("\n",@usage)."\n\nERROR: '-tree st' only allows the setting '-seq_type prot'.\n\n";
			exit;
		} else {
		$treesig=1;
		}
	} elsif ($tree_type eq "snp") {
		unless ($seqtype eq "nucl") {
			print join("\n",@usage)."\n\nERROR: '-tree snp' only allows the setting '-seq_type nucl'.\n\n";
			exit;
		} else {
		$treesig=1;
		}
	} elsif ($seqtype eq "sm") {
	
	} else {
		print join("\n",@usage)."\nERROR: Argument '-tree $tree_type': the value '$tree_type' cannot be recognized (requested to be 'st', 'sm' or 'snp').\n\n";
		exit;
	}
} 

my $qwer;
$tblastn=$opt{"blast_dir"} if exists($opt{"blast_dir"});
if (exists($opt{"iden_cutoff"})) {
	$qwer=$opt{"iden_cutoff"};
	if ($opt{"iden_cutoff"} =~ /^\d+$/) {
		if ($opt{"iden_cutoff"}>=50 || 100>=$opt{"iden_cutoff"}) {
			$blastIdentitycutoff=$opt{"iden_cutoff"};
		} else {
			print join("\n",@usage)."\nERROR: Argument '-iden_cutoff $qwer': the value '$qwer' does not fall in recommended interval 50-100.\n\n";
			exit;
		}
	} else {
		print join("\n",@usage)."\nERROR: Argument '-iden_cutoff $qwer': '$qwer' is not an integer.\n\n";
		exit;
	}
}

my $rewq;
if (exists($opt{"gene_cutoff"}))	 {
	$rewq=$opt{"gene_cutoff"};
	if ($opt{"gene_cutoff"} =~ /^0?\.?\d+$/) {
		if ($opt{"gene_cutoff"}>=0.5 || 1>=$opt{"gene_cutoff"}) {
			$geneCutoff= $opt{"gene_cutoff"};
		} else {
			print join("\n",@usage)."\nERROR: Argument '-gene_cutoff $rewq': the value '$rewq' does not fall in recommended interval 0.5-1.\n\n";
			exit;
		}
	} else {
		print join("\n",@usage)."\n\nERROR: Argument '-gene_cutoff $rewq': '$rewq' is not a decimal.\n\n";
		exit;
	}
}
if ($treesig) {
	$geneCutoff=1;
	print "\nWARNING: '-tree $tree_type' requires '-gene_cutoff' to be set as '1'.\nThe value '1' will be used instead.\n";
	print LOG "\nWARNING: '-tree $tree_type' requires '-gene_cutoff' to be set as '1'.\nThe value '1' will be used instead.\n";		
} 

my $uyt;
if (exists($opt{"genome_cutoff"})) {
	$uyt=$opt{"genome_cutoff"};
if ($opt{"genome_cutoff"} =~ /^0?\.?\d+$/) {
		if ($opt{"genome_cutoff"}>=0.5 || 1>=$opt{"genome_cutoff"}) {
			$genomeCutoff= $opt{"genome_cutoff"};
		} else {
			print join("\n",@usage)."\nERROR: Argument '-ggenome_cutoff $uyt': the value '$uyt' does not fall in recommended interval 0.5-1.\n\n";
			exit;
		}
	} else {
		print join("\n",@usage)."\nERROR: Argument '-genome_cutoff $uyt': '$uyt' is not a decimal.\n\n";
		exit;
	}
}	
	
my $BNM;
if (exists($opt{"thread"})) {	
	$BNM=$opt{"thread"};
	if ($opt{"thread"} =~ /^\d+$/) {
			$num_threads=$opt{"thread"};
	} else {
		print join("\n",@usage)."\nERROR: Argument '-thread $BNM': '$BNM' is not an integer.\n\n";
		exit;
	}
}

print "\n\n###############Starting Informations ###############\n\n";
print LOG "\n\n###############Starting Informations ###############\n\n";

print LOG '====== EasyCGTree ======
	Version 3.1
by Dao-Feng Zhang (Hohai University, China)

';
print '====== EasyCGTree ======
	Version 3.1
by Dao-Feng Zhang (Hohai University, China)

';

print LOG "#########################\n\nOptions: \n";  
print "#########################\n\nOptions: \n";

my @np =("mode", "tree", "thread", "blast_dir", "iden_cutoff", "gene_cutoff", "genome_cutoff");

foreach my $kk (keys %opt) {
	next if $kk ~~ @np;
	print LOG "-$kk $opt{$kk}\n";  
	print "-$kk $opt{$kk}\n";
}
print "-mode $mode\n-tree $tree_type\n-thread $num_threads\n-blast_dir $tblastn\n-iden_cutoff $blastIdentitycutoff\n-gene_cutoff $geneCutoff\n-genome_cutoff $genomeCutoff\n";
print LOG "-mode $mode\n-tree $tree_type\n-thread $num_threads\n-blast_dir $tblastn\n-iden_cutoff $blastIdentitycutoff\n-gene_cutoff $geneCutoff\n-genome_cutoff $genomeCutoff\n";

print LOG "\n#########################\n\nJob Started at: $hour:$min:$sec,$ryear-$rmon-$mday\n\n";  
print "\n#########################\n\nJob Started at: $hour:$min:$sec,$ryear-$rmon-$mday\n\n";

if ($mode eq "A") {
	print "The user ordered a complete analysis (type A) of the pipeline!\n\n";
	print LOG "The user ordered a complete analysis (type A) of the pipeline!\n\n";
} elsif ($mode eq "A1") {
	print "The user ordered a first-part analysis (type A1) of the pipeline!\n\n";
	print LOG "The user ordered a first-part analysis (type A1) of the pipeline!\n\n";
} elsif ($mode eq "A2") {
	print "The user ordered a second-part analysis (type A2) of the pipeline!\n\n";
	print LOG "The user ordered a second-part analysis (type A2) of the pipeline!\n\n";
} elsif ($mode eq "A3") {
	print "The user ordered a complete analysis (type A3) with a pre-defined core-gene set of some genomes, which will be also used in current analysis!\n\n";
	print LOG "The user ordered a complete analysis (type A3) with a pre-defined core-gene set of some genomes, which will be also used in current analysis!\n\n";
}

my $TEMdir= $inputDir."_TEM";
unless ($mode eq "A2") {

###############Step 1: Make BlastDB #################
	print "\n############### Step 1: Make BlastDB ###############\n\n";
	print LOG "\n############### Step 1: Make BlastDB ###############\n\n";

	mkdir "$TEMdir" || die " Permission denied to create directories. Please contact the Administrator of your System\n"; 
	
	mkdir "$TEMdir/TEM1_blastDB";
	my $genomeNum = 0;
	
	print "Step 1: Making BLASTDB from directory '$inputDir': \n";
	foreach my $file (@query) {
		next unless $file =~ /\.fas$/;
		$genomeNum++;
		my $cmd;
		if ($tblastn =~ /\w/) {
			$cmd = "$tblastn"."makeblastdb -in ./$inputDir/$file -dbtype $seqtype -out ./$TEMdir/TEM1_blastDB/$file";
		} else {
			$cmd = "makeblastdb -in ./$inputDir/$file -dbtype $seqtype -out ./$TEMdir/TEM1_blastDB/$file";
		}
		my $response1= `$cmd` || die "Step 1: Can't execute 'makeblastdb': $!\nPlease make sure: 'makeblastdb.exe' is present in '$tblastn'; or the path was added to the environment variable.\n";  

		my @res1=split /\n/, $response1;
		foreach my $res1 (@res1) {
			next unless $res1 =~ /error/i;
			print LOG "Step 1 ERROR: Something goes wrong when execute 'makeblastdb':\n$res1\n\n\n\nEasyCGTree stopped because of problematic input data: ./$inputDir/$file\n";
			die  "Step 1 ERROR: Something goes wrong when execute 'makeblastdb':\n$res1\n\n\n\nEasyCGTree stopped because of problematic input data: ./$inputDir/$file\n";
			}
		
		print "$genomeNum ";
	}
	print "\n\n";
	print "=====$genomeNum BLAST Database were created from $inputDir!=====\n";
	print LOG "=====$genomeNum BLAST Database were created from $inputDir!=====\n";
	my $refN=0;
	if ($mode eq "A3") {
		print "Step 1: Making BLASTDB from directory '$Reference': \n";
		foreach my $file (@queryR) {
			next unless $file =~ /\.fas$/;
			$genomeNum++;
			$refN++;
			my $cmd;
			if ($tblastn =~ /\w/) {
				$cmd = "$tblastn"."makeblastdb -in ./$Reference/$file -dbtype $seqtype -out ./$TEMdir/TEM1_blastDB/$file";
			} else {
				$cmd = "makeblastdb -in ./$Reference/$file -dbtype $seqtype -out ./$TEMdir/TEM1_blastDB/$file";
			}

			my $response1= `$cmd` || die "Step 1: Can't execute 'makeblastdb': $!\nPlease make sure: 'makeblastdb.exe' is present in '$tblastn'; or its path was added to the environment variable.\n";  
			my @res1=split /\n/, $response1;
			foreach my $res1 (@res1) {
				next unless $res1 =~ /error/i;
				print LOG "Step 1 ERROR: Something goes wrong when execute 'makeblastdb':\n$res1\n\n\n\nEasyCGTree stopped because of problematic input data: ./$Reference/$file\n";
				die  "Step 1 ERROR: Something goes wrong when execute 'makeblastdb':\n$res1\n\n\n\nEasyCGTree stopped because of problematic input data: ./$Reference/$file\n";
			}
			print "$refN ";	
		}
	print "\n\n";
	print "===== $refN BLAST Database were created from '$Reference'!=====\n";
	print LOG "===== $refN BLAST Database were created from '$Reference'!=====\n";
	}
print "\n=====Totally, $genomeNum BLAST Database were created in '$TEMdir/TEM1_blastDB'!=====\n";
print LOG "\n=====Totally, $genomeNum BLAST Database were created in '$TEMdir/TEM1_blastDB'!=====\n";
	
############### Step 2: Blast Search #################
print "\n\n############### Step 2: Blast Search ###############\n\n";
print LOG "\n\n############### Step 2: Blast Search ###############\n\n";
	opendir(DIR2, $queryDir) || die "Step 1: Can't open input directory '$queryDir': $!\n";
	my @query2 = readdir(DIR2);
	closedir(DIR2);
	die "Step 2: Input Query directory $queryDir does not contain any files.\n" unless scalar(@query2);
	my $num2=0;
	open (QUERY, ">$queryDir.seq") or die "Step 1: Can't open input file '$queryDir.seq': $!\n";
	foreach my $file2 (@query2) {
		next unless $file2 =~ /\.fas$/;
		$num2++;
		open (FILE2, "<$queryDir/$file2") or die "Step 3: Can't open '$queryDir/$file2': $!\n";
		foreach my $in2 (<FILE2>) {
			next unless $in2;
			if ($in2 =~ /\>/) {
				print QUERY "$in2";
			} else {
				$in2 =~ s/\*//;
				while (1) {
					if (length($in2)>60) {
						my $sub = substr($in2,0,60);
						print QUERY "$sub\n";	
						$in2 =~ s/$sub//;
					} else {
						print QUERY "$in2\n";
						last;
					}
				}
			}
		}
		close FILE2;
	}
	close QUERY;
	print "The query file '$queryDir.seq' have been created!\n\n";
	print LOG "The query file '$queryDir.seq' have been created in current directory!\n\n";
	mkdir "$TEMdir/TEM2_blastOUT";
	
	my $query = "$queryDir.seq"; 
	my $number = 0;
	print "Doing BLAST search ($genomeNum in total) against '$inputDir': ";

	foreach my $file (@query) {
		next unless $file =~ /\.fas$/;
		my $cmd;
		if ($tblastn =~ /\w/) {
			$cmd="$tblastn"."$blast -query $query -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
		} else {
			$cmd="$blast -query $query -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
		}

		system $cmd || die "Step 2: Can't execute 'blastp/tblastn': $!\nPlease make sure: 'blastp/tblastn' is present in '$tblastn'; or their path was added to the environment variable.\n";
		$number++;
		print "$number ";
	}
	print "\n";
	print LOG "\n";
	
	print "===== $number BLAST searches against '$inputDir' has been performed!=====\n";
	print LOG "===== $number BLAST searches against '$inputDir' has been performed!=====\n";
	my $number0 = 0;
	print "\n";
	if ($mode eq "A3")	{
		print "Step 2: BLAST search ($genomeNum in total) against '$Reference': ";

		foreach my $file (@queryR) {
			next unless $file =~ /\.fas$/;
			my $cmd;
			if ($tblastn =~ /\w/) {
				$cmd="$tblastn"."$blast -query ./$Reference/$file -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/Reference_$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
			} else {
				$cmd="$blast -query ./$Reference/$file -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/Reference_$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
			}
			$number++;
			system $cmd || die "Step 2: Can't execute 'blastp/tblastn': $!\nPlease make sure: 'blastp/tblastn' is present in '$tblastn'; or their path was added to the environment variable.\n";
			$number0++;
			my $total= $#queryR + 1;
			print "$number0 ";
		}
		print "\n\n";
	print LOG "===== $number0 BLAST searches against '$Reference' has been performed!=====\n";
	print "===== $number0 BLAST searches against '$Reference' has been performed!=====\n";
	}

	print "=====Totally, $number BLAST search results has been written to '$TEMdir/TEM2_blastOUT'!=====\n";
	print LOG "=====Totally, $number BLAST search results has been written to '$TEMdir/TEM2_blastOUT'!=====\n";

############### Step 3: Filtrate the Blast Results ################# 	
	print "\n\n\n############### Step 3: Filtrate Blast Results ###############\n\n";
	print LOG "\n\n\n############### Step 3: Filtrate Blast Results ###############\n\n";
	mkdir "$TEMdir/TEM3_blastOUT_S";
	
	opendir(DIR3, "$TEMdir/TEM2_blastOUT") || die "Step 1: Can't open input directory '$TEMdir/TEM2_blastOUT': $!\n";
	my @query3 = readdir(DIR3);
	closedir(DIR3);
	die "Step 3: Input  directory $TEMdir/TEM2_blastOUT does not contain any files.\n" unless scalar(@query3);
	foreach my $file3 (@query3) {
		my (%out3,%score3);
		next unless $file3 =~ /\.fas$/;                    
		open (FIL3, "<$TEMdir/TEM2_blastOUT/$file3") or die "Can't open '$TEMdir/TEM2_blastOUT/$file3': $!\n";
		my (@key3,@result3);
		my $geneC=0;
		foreach my $in3 (<FIL3>) {
			$in3=~ s/\n// unless $in3=~s/\r\n//;
			next unless $in3;
			my @in3 =split /\s+/, $in3;                           
			next if $in3[4] < $blastIdentitycutoff;               
			my @in33 =split /_/, $in3[0];
			push @key3, $in33[3];									
			unless ($out3{$in33[3]}) {
				$out3{$in33[3]}=$in3;
				$score3{$in33[3]}=$in3[7];
				$geneC++;
			} else {
				$out3{$in33[3]}=$in3 if $score3{$in33[3]}<$in3[7];
				$score3{$in33[3]}=$in3[7] if $score3{$in33[3]}<$in3[7];
			}
		}
		close FIL3;
		my $outfile3="$geneC"."__$file3";
		
		open(OU3, ">$TEMdir/TEM3_blastOUT_S/$outfile3") or die "Step 3: Can't open '$TEMdir/TEM3_blastOUT_S/$outfile3': $!\n";
		foreach my $in (keys %out3) {
			print OU3 "$out3{$in}\n";
		}
		close OU3;
	}
	print "=====Screened BLAST search results has been written to '$TEMdir/TEM3_blastOUT_S'!=====\n";
	print LOG " =====Screened BLAST search results has been written to '$TEMdir/TEM3_blastOUT_S'!=====\n";	
}

if ($mode eq "A1") {
	print "\nThe user ordered first-part (A1) analysis, and it is finished!\n\n";
	print LOG "\nThe user ordered first-part (A1) analysis, and it is finished!\n\n";
} else {
opendir(DIR, "$TEMdir/TEM3_blastOUT_S") || die "Step 3: Can't open input directory '$TEMdir/TEM3_blastOUT_S': $!\n";
my @query3 = readdir(DIR);
closedir(DIR);
die "Step 3: Input directory TEM3_blastOUT_S does not contain any files.\n" unless scalar(@query3);

my (@querylist0,@querylistN);
my (@finalgenelist,@outgenelist,@finalgenome,@outgenome);
my (@finalfile);
my $geneNum=0;
my $genomeNum;

foreach my $kk (@query3) {
	next unless $kk =~ /.\.fas$/;
	$genomeNum++;
	my @yy=split /__/,$kk;
	$geneNum = $yy[0] if $geneNum < $yy[0]
}
my $cutoff2=int($geneNum*$genomeCutoff);     ####################
print "\n These genomes harboring more than $cutoff2 ($geneNum * $genomeCutoff) genes, which will be used in following analysis. This is the list (genome/gene_number):\n";
print LOG "\n These genomes harboring more than $cutoff2 ($geneNum * $genomeCutoff) genes, which will be used in following analysis. This is the list (genome/gene_number):\n";

foreach my $kk (@query3) {
	next unless $kk =~ /.\.fas$/;
	my @yy=split /__/,$kk;
	if ($yy[0]>=  $cutoff2) {
		push @finalfile, $kk ;
		push @finalgenome, $yy[1];
		print "$yy[1]/$yy[0] ";
		print LOG "$yy[1]/$yy[0] ";
	} else {
		push @outgenome, "$yy[1]/$yy[0]";
	}
}

my $dis=$#outgenome+1;
my $ingenomeN=$#finalgenome+1;	

print ".\n=====Totally, $ingenomeN of $genomeNum genomes were selected!=====\n\n";
print LOG ".\n=====Totally, $ingenomeN of $genomeNum genomes were selected!=====\n\n";
if (@outgenome) {
print "The following $dis genomes were excluded(genome/gene_number): @outgenome.\n\n";
print LOG "The following $dis genomes were excluded (genome/gene_number): @outgenome.\n\n";
}

my $num=0;
foreach my $file3 (@finalfile) {
	open (FIL3, "<$TEMdir/TEM3_blastOUT_S/$file3") or die "Step 3: Can't open '$TEMdir/TEM3_blastOUT_S/$file3': $!\n";
	foreach my $in3 (<FIL3>) {
		$in3=~ s/\n// unless $in3=~s/\r\n//;
		next unless $in3;
		my @in3 =split /\s/, $in3;
		my @in31 =split /_/, $in3[0];
		unless (@querylist0) {
			push @querylist0,$in31[3];
			$querylistN[$num]++;
			$num++;	
		} else {
			unless ($in31[3] ~~ @querylist0) {
				push @querylist0,$in31[3];
				$querylistN[$num]++;
					$num++;
			} else {
				foreach my $inn3 (0..$#querylist0) {
					if ($in31[3] eq $querylist0[$inn3]) {
					
						$querylistN[$inn3]++;	
					}
				}
			}
		}
	}
	close FIL;
}
my $cutoff=int($ingenomeN*$geneCutoff);           ##############

print "\nThese genes were present in >= $cutoff ($ingenomeN * $geneCutoff) of the $ingenomeN selected genomes, which will be used for tree inference (gene/prevalence):\n";
print LOG "\nThese genes were present in >= $cutoff ($ingenomeN * $geneCutoff) of the $ingenomeN selected genomes, which will be used for tree inference (gene/prevalence):\n";

foreach my $kk (0..$#querylistN) {
	if ($querylistN[$kk]<  $cutoff) {
		push @outgenelist, "$querylist0[$kk]/$querylistN[$kk]";
	} else {
		print "$querylist0[$kk]/$querylistN[$kk] ";
		print LOG "$querylist0[$kk]/$querylistN[$kk] ";
		push @finalgenelist, $querylist0[$kk];
	}
}
my $ingene =$#finalgenelist+1;
my $outgene =$#outgenelist+1;

print ".\n=====Totally, $ingene of $geneNum genes were selected!=====\n";
print LOG ".\n=====Totally, $ingene of $geneNum genes were selected!=====\n";
if (@outgenelist) {
print "\nThe following $outgene/$geneNum genes were excluded (gene/prevalence):\n@outgenelist\n";
print LOG "\nThe following $outgene/$geneNum genes were excluded (gene/prevalence):\n@outgenelist\n";
}

############### Step 4: Retrieve Sequences from Each Genome/proteome ################# 
print "\n\n\n############### Step 4: Retrieve Sequences from Each Genome/proteome ###############\n\n";
print LOG "\n\n\n############### Step 4: Retrieve Sequences from Each Genome/proteome ###############\n\n";

my $numg=0;
mkdir "$TEMdir/TEM4_GeneSeqs";

print "According to the Step 3, $ingenomeN genomes harboring >= $cutoff2 of the $ingene common genes, will be used in following analysis. The related gene sequences will be extracted.\n";
print LOG "According to the Step 3, $ingenomeN genomes harboring >= $cutoff2  of the $ingene common genes, will be used in following analysis. The related gene sequences will be extracted.\n******See more details in Step 3.******\n";

foreach my $i1 (0..$#finalfile) {
	$numg++;
	if ($seqtype eq "prot") {    
		my %geno;
		my ($seq,$seqID,$inname);
		$inname = $finalgenome[$i1];
		if ($finalgenome[$i1]=~ /^Reference_/) {
			$inname=~ s/^Reference_//;
			open (FILE, "<$Reference/$inname") or die "Step 4: Can't open '$Reference/$inname': $!\n";
		} else {
			open (FILE, "<$inputDir/$finalgenome[$i1]") or die "Step 4: Can't open '$inputDir/$finalgenome[$i1]': $!\n";########
		}
		foreach my $i3 (<FILE>) {
			$i3=~ s/\n// unless $i3=~ s/\r\n//;
			next unless $i3;		
			if ( $i3 =~ s/\>//) {
				if ($seq) {
					$geno{$seqID} = $seq;
					undef $seq;
				}
				my @ii = split /\s/, $i3;
				$seqID = $ii[0];			
			} else { 
				$seq .= $i3 if $seq;
				$seq = $i3 unless $seq;
			}
		}
		$geno{$seqID} = $seq;
		close FILE;
			
		my $seqnum=0;
		open(OUT, ">$TEMdir/TEM4_GeneSeqs/$finalgenome[$i1]") or die "Step 4: Can't open '$TEMdir/TEM4_GeneSeqs/$finalgenome[$i1]': $!\n";
		open (FIL, "<$TEMdir/TEM3_blastOUT_S/$finalfile[$i1]") or die "Step 4: Can't open 'TEM3_blastOUT_S/$finalfile[$i1]': $!\n";
		my $seqTag="1_1_$inname";		
		$seqTag=~ s/\.fas/_/;		
		foreach my $i2 (<FIL>) {
			$i2=~ s/\n// unless $i2=~s/\r\n//;
			next unless $i2;
			my @in0 = split /\s/, $i2;
			my @ff =split /_/, $in0[0];	
			my $geneSeq = $geno{$in0[1]};			
			print OUT ">$seqTag"."$ff[3]\n$geneSeq\n";
			$seqnum++;
		}
		close FIL;	
			
		close OUT;
		print "$numg: Got $seqnum/$geneNum sequences from genome: $finalgenome[$i1]\n";
	} else {
		my %geno;
		my ($seq,$seqID);
		open (FILE, "<$inputDir/$finalgenome[$i1]") or die "Step 4: Can't open '$inputDir/$finalgenome[$i1]': $!\n";########
		foreach my $i3 (<FILE>) {
			$i3=~ s/\n// unless $i3=~ s/\r\n//;
			next unless $i3;		
			if ( $i3 =~ s/\>//) {
				if ($seq) {
					$geno{$seqID} = $seq;
					undef $seq;
				}
				my @ii = split /\s/, $i3;
				$seqID = $ii[0];			
			} else { 
				$seq .= $i3 if $seq;
				$seq = $i3 unless $seq;
			}
		}
		$geno{$seqID} = $seq;
		close FILE;
			
		my $seqnum=0;
		open(OUT, ">$TEMdir/TEM4_GeneSeqs/$finalgenome[$i1]") or die "Step 4: Can't open '$TEMdir/TEM4_GeneSeqs/$finalgenome[$i1]': $!\n";
		open (FIL, "<$TEMdir/TEM3_blastOUT_S/$finalfile[$i1]") or die "Step 4: Can't open 'TEM3_blastOUT_S/$finalfile[$i1]': $!\n";
		my $seqTag="1_1_$finalgenome[$i1]";		
		$seqTag=~ s/\.fas/_/;					
		foreach my $i2 (<FIL>) {
			$i2=~ s/\n// unless $i2=~s/\r\n//;
			next unless $i2;
			my @in0 = split /\s/, $i2;
			my@ff =split /_/, $in0[0];
			my ($start, $end,$length);
			if ($in0[2] < $in0[3]) {
				if ($in0[2]-1>0) {
					$start=$in0[2]-1;
				} else {
					$start=0;
				}
				if (length($geno{$in0[1]}) > $in0[3]+1) {
					$length =$in0[3]-$in0[2]+1;
				} else {
					$length =length($geno{$in0[1]})-$start;
				}	
			} else {
				if ($in0[3]-1>0) {
					$start=$in0[3]-1;
				} else {
					$start=0;
				}
				if (length($geno{$in0[1]}) > $in0[2]+1) {
					$length =$in0[2]-$in0[3]+1;
				} else {
					$length =length($geno{$in0[1]})-$start;
				}
			}	
			my $geneSeq = substr($geno{$in0[1]},$start, $length);	
			$geneSeq = reverse_complement ($geneSeq) if ($in0[2] > $in0[3]);				
			print OUT ">$seqTag"."$ff[3]\n$geneSeq\n";
			$seqnum++;
		}
		close FIL;	
			
		close OUT;
		print "	$numg: Got $seqnum/$geneNum sequences from genome:$finalgenome[$i1]!\n";
	}
}
print "=====\nTotally, $numg sequence files has been written to '$TEMdir/TEM4_GeneSeqs'!=====\n";
print LOG "=====\nTotally, $numg sequence files has been written to '$TEMdir/TEM4_GeneSeqs'!=====\n";

############### Step 5: Gather the Orthologs in Single File #################
print "\n\n\n############### Step 5: Gather the Orthologs in Single File ###############\n\n";
print LOG "\n\n\n############### Step 5: Gather the Orthologs in Single File ###############\n\n";

opendir(DIR1, "$TEMdir/TEM4_GeneSeqs") || die "Step 5: Can't open input directory '$TEMdir/TEM4_GeneSeqs': $!\n";
my @in1 = readdir(DIR1);
closedir(DIR1);
die "Step 5: Input directory '$TEMdir/TEM4_GeneSeqs' does not contain any files.\n" unless scalar(@in1);

mkdir "$TEMdir/TEM5_GeneCluster";
my %outClu;
my $num5=0;
print "Step 5: Reading the";
foreach my $in5 (@finalgenome) {	
	my $tt =$in5;
	$tt =~ s/\.fas$//;
	open (FILE, "<$TEMdir/TEM4_GeneSeqs/$in5") or die "Step 5: Can't open '$TEMdir/TEM4_GeneSeqs/$in5': $!\n";
	my ($id,$sig);
	foreach my $i1 (<FILE>) {
		$i1=~ s/\n// unless $i1=~s/\r\n//;
		next unless $i1;
		if ( $i1 =~ /\>/ ) {
			my @seqids=split /_/, $i1;
			$id =$seqids[3];
			$sig=1 if $id ~~ @finalgenelist;     
		} elsif ($sig) {
			unless ($outClu{$id}) {
				$outClu{$id}= ">$tt\n$i1\n";
			} else {
				$outClu{$id} .= ">$tt\n$i1\n"; 
			}
			$sig=0;
		}
	}
	close FILE;
	$num5++;
	print " $num5";
}

my $n=0;
print " files.\nStep 5: Writing the files (file no./sequence no.):";
foreach my $in (keys %outClu) {
	my $cnt;
	$n++;
	open (OUT, ">./$TEMdir/TEM5_GeneCluster/$in.fas") or die "Step 5: Can't open '$TEMdir/TEM5_GeneCluster/$in.fas': $!\n";
	print OUT "$outClu{$in}";
	$outClu{$in} =~ s/>/"X".$cnt++/ge;
	print " $n/$cnt";
	close OUT;
}
print ".\n";	
print LOG "Totally, $n gene clusters were written to '$TEMdir/TEM5_GeneCluster'.\n";

############### Step 6: Deduplicate Sequences #################
print "\n\n\n############### Step 6: Deduplicate Sequences ###############\n\n";
print LOG "\n\n\n############### Step 6: Deduplicate Sequences ###############\n\n";
mkdir "$TEMdir/TEM6_DedupGeneCluster";
mkdir "$TEMdir/TEM60_DedupList";

opendir(DIR1, "$TEMdir/TEM5_GeneCluster") || die "Step 6: Can't open input directory '$TEMdir/TEM5_GeneCluster': $!\n";
my @in6 = readdir(DIR1);
closedir(DIR1);
die "Step 6: Input directory $TEMdir/TEM5_GeneCluster does not contain any files.\n" unless scalar(@in6);
print "Step 6: Reading the files (file no./sequence types):";

my %outDedup;
my $num6=0;
foreach my $inx (@finalgenelist) {
	my $in =$inx;
	my $genome=$in;
	 $in .="\.fas";
	open (FILE, "<$TEMdir/TEM5_GeneCluster/$in") or die "Step 6: Can't open '$TEMdir/TEM5_GeneCluster/$in': $!\n";
	my (@id,@seq,$refID,$refSeq);
	foreach my $i1 (<FILE>) {
		$i1=~ s/\n// unless $i1=~s/\r\n//;
		next unless $i1;
		if ( $i1 =~ /\>/ ) {
			push @id, $i1;
		} else {
			push @seq, $i1;
		}
	}
	close FILE;
	die "Step 6: There must be something wrong in file $in!\n" unless $#id==$#seq;
	open (OU1, ">./$TEMdir/TEM6_DedupGeneCluster/$genome.fas") or die "Step 6: Can't open '$TEMdir/TEM6_DedupGeneCluster/$genome.fas': $!\n";
	open (OU2, ">./$TEMdir/TEM60_DedupList/$genome.txt") or die "Step 6: Can't open '$TEMdir/TEM60_DedupList/$genome.txt': $!\n";
	my $loop;
	
	while (1) {
		last unless @id;
		$refID= shift @id;
		$refSeq=shift @seq;
		my $id=$refID;
		my @splice;
		foreach my $i1 (0..$#id) {
			if ($seq[$i1] eq $refSeq) {
				$id .="$id[$i1]";
				push @splice, $i1;
			}
		}
		@splice=reverse @splice;		
		foreach my $in4 (@splice) {
			splice @seq, $in4, 1;
			splice @id, $in4, 1;
		}
		print OU1 "$refID\n";
		while (1) {
			$refSeq=~s/\*//g;
			if (length($refSeq)>60) {
				my $sub = substr($refSeq,0,60);
				print OU1 "$sub\n";	
				$refSeq =~ s/$sub//;		
			} else {
				print OU1 "$refSeq\n";last;
				}
		}
		print OU2 "$id	$refID\n";
		$loop++;
	}
	$num6++;
	print " $num6/$loop";
}
print ".\n";
print LOG "\nTotally, $num6 files has been processed successfully.\n Related records has been written to '$TEMdir/TEM60_DedupList' and '$TEMdir/TEM6_DedupGeneCluster'.\n";

############### Step 7: Do Alignment #################
print "\n\n\n############### Step 7: Do Alignment ###############\n\n"; 
print LOG "\n\n\n############### Step 7: Do Alignment ###############\n\n";
mkdir "$TEMdir/TEM7_Alignment";
my $num7 = 0;

foreach (@finalgenelist) {
	my $program = "clustalo";##Different from Windows version
	my $outname = $_ . ".fas.fasta";
	my @para = ("-i", "./$TEMdir/TEM6_DedupGeneCluster/$_.fas", "-o", "./$TEMdir/TEM7_Alignment/$outname", "--outfmt=fasta", "--output-order=tree-order","--threads=$num_threads","-v","--force");###Different from Windows version
	unless (system './bin/clustalo', @para) {##Different from Windows version
		$num7++;
		print "\nStep 7: $_ has been completed! It is the $num7/$ingene file.\n\n";	
	} else {
		die "Step 7: Can't execute './bin/clustalo': $!.\n";   ##Different from Windows version
	}
}
print LOG "\n$num7 alignments has been created in '$TEMdir/TEM7_Alignment'!\n\n";

############### Step 8: Trim Sequence Algnments #################
print "\n\n\n############### Step 8: Trim Sequence Algnments ###############\n\n";
print LOG "\n\n\n############### Step 8: Trim Sequence Algnments ###############\n\n";
mkdir "$TEMdir/TEM8_AlnTrimmed";

opendir(DIR, "$TEMdir/TEM7_Alignment") || die "Step 8: Can't open input directory '$TEMdir/TEM7_Alignment':$!\n";
my @in8 = readdir(DIR);
closedir(DIR);
die "Step 8: Input directory $TEMdir/TEM7_Alignment does not contain any files.\n" unless scalar(@in8);

my (@problem,@error);
my $num8 = 0;
print "Step 8: Reading the files:";
foreach my $filex (@finalgenelist) {
	my $file = $filex;
	$file .= "\.fas\.fasta";
	(my $outname = $file) =~ s/\.fas//;
	open (FILE, "$TEMdir/TEM7_Alignment/$file") or die "Step 8: Can't open '$TEMdir/TEM7_Alignment/$file': $!\n";
	my ($seq, @seqs);
	while (<FILE>) {
		if (/\A\>/) { 
			if ($seq) {
				$seq .= "\n";
				push @seqs, ($seq,$_);
				$seq = undef;
			} else {
				push @seqs, $_;
			}
		} else {
			chomp;
			$seq .= $_;}
	}
	$seq .="\n";
	push @seqs, $seq;
	close FILE;

	my %seqs = @seqs;
	my $length =0;
	my (@posF,@posE);
	foreach (sort keys %seqs) {
		$seqs{$_} =~ /\A-*(?<start>[A-Z])/;
		my $fpos = index($seqs{$_}, $+{start});
		push @posF, $fpos;
		$seqs{$_} =~ /(?<end>[A-Z])\*?-*\Z/;
		my $epos = rindex($seqs{$_}, $+{end});
		push @posE, $epos;
		next if $length;
		$length = length($seqs{$_}) if $length < length($seqs{$_});
	}
	@posF = sort @posF;
	@posE=sort @posE;
	my $t8=int($#posF/2+1);
	my $fore_pos = $posF[$t8];
	my $end_pos = $posE[$t8];

	open (OUT, ">$TEMdir/TEM8_AlnTrimmed/$outname") or die "Step 8: Can't open '$TEMdir/TEM8_AlnTrimmed/$outname': $!\n";
	foreach (sort keys %seqs) {
		$seqs{$_} = substr($seqs{$_},$fore_pos, $end_pos - $fore_pos + 1);
		print OUT "$_$seqs{$_}\n";
	}

	close OUT;
	$num8++;
	print " $num8";
	my $ratio =($end_pos - $fore_pos + 1)/$length;
	push @problem, "$file $ratio\n" if $ratio < 0.5;
	push @error, "$file $ratio\n" if $ratio < 0;
}
print LOG "Totally, $num8 alignments has been trimmed successfully and written to '$TEMdir/TEM8_AlnTrimmed'!\n\n";
if (@problem) {
	my $pro=$#problem +1;

	print ".\nStep 8: These $pro cluster(s) are problemic (less than half length was retained for following analysis), please check related file(s) in 'TEM7_Alignment':\n*****Indicated by the filenames and the coverages*****\n@problem\n";
	print LOG ".\nStep 8: These $pro cluster(s) are problemic (less than half length was retained for following analysis), please check related file(s) in 'TEM7_Alignment':\n*****Indicated by the filenames and the coverages*****\n@problem\n";
}
if (@error) {
	my $err=$#error +1;
	print "These $err +1 clusters are ERROR:\n@error\n";
	print LOG "These $err +1 clusters are ERROR:\n@error\n";
}

############### Step 9: ReCall Sequences #################
print "\n\n\n############### Step 9: ReCall Sequences ###############\n\n";
print LOG "\n\n\n############### Step 9: ReCall Sequences ###############\n\n";
mkdir "$TEMdir/TEM9_ReCallAln";

opendir(DIR1, "$TEMdir/TEM8_AlnTrimmed") || die "Step 9: Can't open input directory '$TEMdir/TTEM8_AlnTrimmed': $!\n";
my @in9 = readdir(DIR1);
closedir(DIR1);

die "Step 9: Input directory $TEMdir/TEM7_Alignment does not contain any files.\n" unless scalar(@in9);
print "Step 9: Reading the files(file no./sequence no.):";
my $num9=0;
foreach my $in (@finalgenelist) {
	my $genome=$in;
	open (FILE, "<$TEMdir/TEM8_AlnTrimmed/$in.fasta") or die "Step 9: Can't open '$TEMdir/TEM8_AlnTrimmed/$in.fasta': $!\n";
	my (@id,@seq,$seq,$seq9);
	my $num99=0;
	foreach my $i1 (<FILE>) {
		$i1=~ s/\n// unless $i1=~s/\r\n//;
		next unless $i1;
		if ( $i1 =~ /\>/ ) {
			if ($seq) {push @seq, $seq; $seq=undef;}
			push @id, $i1;
		} else {
			$seq .="$i1\n";
		}
	}
	push @seq, $seq;
	close FILE;
	open (FIL, "<$TEMdir/TEM60_DedupList/$genome.txt") or die "Step 9: Can't open '$TEMdir/TEM60_DedupList/$genome.txt': $!\n";
	open (OU, ">./$TEMdir/TEM9_ReCallAln/$genome.fasta") or die "Step 9: Can't open '$TEMdir/TEM9_ReCallAln/$genome.fasta': $!\n";
	my @current;
	foreach my $i (<FIL>)	{
		$i=~ s/\n// unless $i=~s/\r\n//;
		next unless $i;
		my @in=split /\s/, $i;
		my @in1=split /\>/,$in[0];
		foreach my $i1 (0..$#id) {
			if ($in[1] eq $id[$i1]) {
				foreach my $i2 (1..$#in1) {
					print OU ">$in1[$i2]\n$seq[$i1]";
					my $ll=$in1[$i2];
					push @current, $ll;
					$num99++;
					next if $seq9;
					$seq9=$seq[$i1];
					$seq9=~ s/[A-Z]/-/g;   ############################
				}
			last;
			}
		}
	}
	close FIL;
	$num9++;
	foreach my $i9 (@finalgenome) {
		$i9 =~ s/\.fas//;
		next if $i9~~@current;
		print OU ">$i9\n$seq9";
	}
	close OU;	
	print " $num9/$num99";
}
print ".\n";
print LOG "Totally, $num9 gene clusters has been processed successfully and written to '$TEMdir/TEM9_ReCallAln'!\n\n";

############### Step 10: Concatenate Sequences #################
print "\n\n\n############### Step 10: Concatenate Sequences ###############\n\n";
print LOG "\n\n\n############### Step 10: Concatenate Sequences ###############\n\n";

opendir(DIR, "$TEMdir/TEM9_ReCallAln") || die "Step 10: Can't open input directory '$TEMdir/TEM9_ReCallAln': $!\n";
my @in10 = readdir(DIR);
closedir(DIR);
die "Step 10: Input directory $TEMdir/TEM9_ReCallAln does not contain any files.\n" unless scalar(@in10);

print "Step 10: Reading the files:";
my $num10 = 0;
my (%out,@id,@seq,$id,$seq);
foreach my $filex ( @finalgenelist) {
	my $file ="$filex.fasta";
	unless (@id) {
		open (FILE, "$TEMdir/TEM9_ReCallAln/$file") or die "Step 10: Can't open '$TEMdir/TEM9_ReCallAln/$file': $!\n";
		foreach my $i1 (<FILE>) {
			$i1=~ s/\n// unless $i1=~s/\r\n//;
			next unless $i1;
			if ( $i1 =~ /\>/ ) {
				if ($seq) {
					push @seq, $seq; 
					$seq=undef;
				}
				push @id, $i1;
			} else {
				$seq .="$i1";
			}
		}
		close FILE;
		push @seq,$seq;
		foreach my $i2 (0..$#id) {
			$out{$id[$i2]}=$seq[$i2];
		}
	} else {
		open (FILE, "$TEMdir/TEM9_ReCallAln/$file") or die "Step 10: Can't open '$TEMdir/TEM9_ReCallAln/$file': $!\n";
		foreach my $i1 (<FILE>) {
			$i1=~ s/\n// unless $i1=~s/\r\n//;
			next unless $i1;
			if ( $i1 =~ /\>/ ) {
				$id= $i1;
			} else {
				$out{$id} .="$i1";
			}
		}
	}
 	$num10++;
	print " $num10";
}
print ".\n";

open (OUT, ">$inputDir.concatenation.fas") or die "Can't open '$inputDir.concatenation.fas': $!\n";
my $nn=0;
foreach (sort keys %out) {
	$nn++;
	print OUT "$_\n";
	my $i1=$out{$_};
	while (1) {
		if (length($i1)>60) {
			my $sub = substr($i1,0,60);
			print OUT "$sub\n";	
			$i1 =~ s/$sub//;	
		} else {
			print OUT "$i1\n";
			last;
		}
	}
}
close OUT;
print "Got $nn concatenated sequences written to '$inputDir.concatenation.fas'!\n";
print LOG "Got $nn concatenated sequences written to '$inputDir.concatenation.fas'!\n";

############### Step 11: Construct a Tree #################
print "\n\n\n############### Step 11: Construct a Tree ###############\n\n";
print LOG "\n\n\n############### Step 11: Construct a Tree ###############\n\n";

if ($tree_type eq "sm") {
	print "Infering phylogeny of the supermatrices approach......\n\n";
	my $cmd11;
	if ($seqtype eq "nucl") {
		$cmd11="./bin/FastTreeMP -gtr -nt $inputDir.concatenation.fas > $inputDir.supermatrix.tree";##Different from Windiws version
	} elsif ($seqtype eq "prot") {
		$cmd11="./bin/FastTreeMP $inputDir.concatenation.fas > $inputDir.supermatrix.tree";	##Different from Windiws version
	}
	&FastreeResponse ($cmd11);
	print LOG "The supermatrix tree has been written successfully in '$inputDir.supermatrix.tree'.\n\n";
} elsif ($tree_type eq "st") {
	mkdir "$TEMdir/TEM11_GeneTrees";
	opendir(DIR, "$TEMdir/TEM9_ReCallAln") || die "Step 10: Can't open input directory '$TEMdir/TEM9_ReCallAln': $!\n";
	my @in11 = readdir(DIR);
	closedir(DIR);
	die "Step 11: Input directory $TEMdir/TEM9_ReCallAln does not contain any files.\n" unless scalar(@in11);
	my $num11 = 0;
	foreach my $file11x ( @finalgenelist) {
		my $file11 ="$file11x.fasta";
		$num11++;
		print "Step 11: Infering phylogeny from '$TEMdir/TEM9_ReCallAln/$file11': $num11/$num10.\n";
		my $cmd11;
		if ($seqtype eq "nucl") {
			$cmd11="./bin/FastTreeMP -gtr -nt $TEMdir/TEM9_ReCallAln/$file11 > $TEMdir/TEM11_GeneTrees/$file11.tree";##Different from Windiws version
		} elsif ($seqtype eq "prot") {
			$cmd11="./bin/FastTreeMP $TEMdir/TEM9_ReCallAln/$file11 > $TEMdir/TEM11_GeneTrees/$file11.tree";	##Different from Windiws version
		}
		&FastreeResponse ($cmd11);	
	}
	print LOG "Totally, $num11 phylogenies has been written to '$TEMdir/TEM11_GeneTrees' and 'intree'.\n\n";
	print  "Totally, $num11 phylogenies has been written to '$TEMdir/TEM11_GeneTrees' and 'intree'.\n\n";
	open (OU11, ">intree") or die "Step 9: Can't open 'intree': $!\n";
	foreach my $file11 ( @finalgenelist) {
		$file11 .="\.fasta\.tree";
		open (FIL11, "<$TEMdir/TEM11_GeneTrees/$file11") or die "Step 9: Can't open '$TEMdir/TEM11_GeneTrees/$file11': $!\n";
		foreach my $in11 (<FIL11>) {
			next unless $in11=~ /\d/;
			print OU11 "$in11";
		}
		close FIL11;
	}
	close OU11;
	my $cmd11t;
	$cmd11t="./bin/consenseM";##Different from Windows version
	my $response11= `$cmd11t`;  
	print "$response11\n";
	rename "outtree","$inputDir.supertree.tree" || die "Step 3: Can't rename 'outtree' to '$inputDir.supertree.tree'!\n";
	rename "outfile","$inputDir.supertree.outfile" || die "Step 3: Can't rename 'outfile' to '$inputDir.supertree.outfile'!\n";
	print LOG "The supertree has been written successfully in '$inputDir.supertree.tree'.\n\n";
	
} elsif ($tree_type eq "snp") {
	print "Extracting SNP from '$inputDir.concatenation.fas'.\n\n";
	###### ExtractSNP
	{					
		my $outfile="$inputDir.concatenationSNP.fas";
		my %gene;
		my ($i,$name);
		my (@snp,@n,@gap);
		foreach  (keys %out) {
			my @genes;
			$name = $_;
				$gene{$name} = $out{$_};
				@genes = split (//,$out{$_});
				if (@snp) {
					for ($i=0 ;$i<$#genes ;$i++) {
						if ($genes[$i] eq "-") {
							push @gap, $i;
						}
						if ($snp[$i] ne $genes[$i]) {
								push @n, $i;
						}
					}
				}
				@snp = split (//,$out{$_});
		}
		close FILE;
		my @uniq_n = uniq @n;
		my @uniq_gap = uniq @gap;
		my @loc;
		foreach my $uniq_n (@uniq_n) {
			unless ($uniq_n ~~ @uniq_gap) {
					push @loc, $uniq_n;
			}
		}
		open OUTSNP, '>', $outfile or die "Can't open '$outfile': $!";
		foreach my $key (keys %gene) {
			print OUTSNP "$key\n";
			my @genes = split (//,$gene{$key});
			my $nSNP=0;
			foreach my $loc (@loc) {
				$nSNP++;
				print OUTSNP "$genes[$loc]";
				if ($nSNP==60) {
					print OUTSNP "\n";
					$nSNP=0;
				}
			}
			print OUTSNP "\n";
		}
	}
	close OUTSNP;
	print LOG "The SNPs has been extracted and written to '$inputDir.concatenationSNP.fas'\n\n";
	
	my $cmd11="./bin/FastTreeMP -gtr -nt $inputDir.concatenationSNP.fas > $inputDir.SNP.tree";##Different from Windows version
	print "Infering phylogeny based on SNPs......\n";
	&FastreeResponse ($cmd11);
	print LOG "The phylogeny based on SNPs has been written to '$inputDir.SNP.tree'.\n\n";
}
}

############### Ending Information ###############
print "\n\n\n############### Ending Information ###############\n\n";
print LOG "\n\n\n############### Ending Information ###############\n\n";
my ($sec1,$min1,$hour1,$mday1,$mon1,$year1,$wday1,$yday1,$isdst1) = localtime(time);
my $rmon1=$mon1 +1;
my $ryear1=$year1+1900;

if ($mode eq "A2") {
	print "The user ordered the second-part (A2) analysis.\nJob was finished at: $sec1:$min1:$hour1,$mday1-$rmon1-$ryear1.\n";
	print LOG "The user ordered the second-part (A2) analysis.\nJob was finished at: $sec1:$min1:$hour1,$mday1-$rmon1-$ryear1.\n";
} else {
	print LOG "Job was finished at: $hour1:$min1:$sec1,$ryear1-$rmon1-$mday1.\n";
	print "Job was finished at: $hour1:$min1:$sec1,$ryear1-$rmon1-$mday1.\n";
}

my $sig=0;
my ($sec0,$min0,$hour0,$yday0);
if ($sec1>=$sec) {
	$sec0=$sec1-$sec;
	$sig=0;
} else {
	$sec0=$sec1-$sec+60;
	$sig=1;
}
	$min1=$min1-$sig;
if ($min1>=$min) {
	$min0=$min1-$min;
	$sig=0;
} else {
	$min0=$min1-$min+60;
	$sig=1;
}
	$hour1=$hour1-$sig;
if ($hour1>=$hour) {
	$hour0=$hour1-$hour;
	$sig=0;
} else {
	$hour0=$hour1-$hour+24;
	$sig=1;
}
	$yday1=$yday1-$sig;
if ($yday1>=$yday) {
	$yday0=$yday1-$yday;
	$sig=0;
} else {
	$yday0=$yday1-$yday+365;
	$sig=1;
}	
print "Running time: $yday0 d $hour0 h $min0 min $sec0 sec.\n\n\nEasyCGTree Version 2.0 by Dao-Feng Zhang\n\n"; 
print LOG "Running time: $yday0 d $hour0 h $min0 min $sec0 sec.\n\n\nEasyCGTree Version 2.0 by Dao-Feng Zhang\n\n";
close LOG;

sub FastreeResponse {
	my $cmdsub =shift;
	my $response11= `$cmdsub`;  
	my @res11=split /\n/, $response11;
	print "$response11\n";
	foreach my $res1 (@res11) {
		if  ($res1 =~ /truncated/i && $res1 =~ /be too long/i) {
			print LOG "Step 11 ERROR: EasyCGTree stopped when execute 'FastTree':\n$response11\n\n\n\nPlease check whether the genome/proteome files were formated correctly. \nPlease find more informations in section 3.1 of the Manual.\n";##Different from Linux version
			die  "Step 11 ERROR: EasyCGTree stopped when execute 'FastTree':\n$response11\n\n\n\nPlease check whether the genome/proteome files were formated correctly. \nPlease find more informations in section 3.1 of the Manual.\n\n";##Different from Linux version
		} elsif ($res1 =~ /out of memory/i) {
			print LOG "Step 11 ERROR: EasyCGTree stopped when execute 'FastTree':\n$response11\n\n\n\nThis ERROR occured probably because 'TastTree' needs larger memory space than that of your computer. Please find a more powerful PC or server to complete the run.\n";##Different from Linux version
			die  "Step 11 ERROR: EasyCGTree stopped when execute 'FastTree':\n$response11\n\n\n\nThis ERROR occured probably because 'TastTree' needs larger memory space than that of your computer. Please find a more powerful PC or server to complete the run.\n";##Different from Linux version
		}
	}
}

sub reverse_complement {
        my $dna = shift;
        # reverse the DNA sequence
        my $revcomp = reverse($dna);
        $revcomp =~tr/tgca/TGCA/;
        # complement the reversed DNA sequence
        $revcomp =~tr/ACGT/TGCA/;
		return $revcomp;
}

sub checkDataType {
	my $dna = shift;
	my $out;	
	unless (opendir(DIR, $dna)) {
		$out= "Can't open input directory '$dna': $!\n";
		return $out;
		last;
	}
	my @querysub = readdir(DIR);
	closedir(DIR);
	my @name;
	foreach my $in00 (@querysub) {
		next unless $in00 =~ /\w/;	
		unless ($in00 =~ /\.fas$/) {
			$out="it contain file(s) of which the name ends without '.fas'. Please formate them by using script 'FormatNames.pl'.\n";
			last;
		} else {
			push @name, $in00;
		}
	} 
	if ($out) {
		return $out;
	} else {
		unless (@name) {
			$out="Directory $dna does not contain any files.\n" ;
			return $out;
		}else {
			my @geneloc;
			while (1) {
				my $nn =int(rand($#name+1));
				push @geneloc, $nn unless $nn ~~ @geneloc;
				last if $#geneloc == $#name || $#geneloc == 12;
			}
			my @gene;	
			foreach my $file (@geneloc) {
				open (FILE, "./$dna/$name[$file]");
				my $type;
				my $count=0;	
				my $sig=0;
				my $seq;
				foreach my $in (<FILE>) {
					$in=~ s/\n// unless $in=~s/\r\n//;
					if ($in =~ /\>/) {
						$count++;
						$sig=1;
					} else {
						unless ($seq) {
						$seq="$in";
						} else {
							$seq .="$in";
						}
					}
					last if $count ==10;
				}
				close FILE;
				unless ($sig) {
					$out = "./$dna/$name[$file] seems not a fasta-formated file. Please check others.\n";
					last;
				} else {
					my $len = length($seq);
					my ($cntA,$cntG,$cntC,$cntT);
					$seq =~ s/A/"X".$cntA++/ge;
					$seq =~ s/G/"X".$cntG++/ge;
					$seq =~ s/C/"X".$cntC++/ge;
					$seq =~ s/T/"X".$cntT++/ge;
					my $rio = ($cntA+$cntG+$cntC+$cntT)/$len;
					if ($rio>0.9) {
						$type = "nucl";	
					} else {
						$type ="prot";
					}
					$out =$type unless $out;
				}
				unless ($type eq $out) {
					$out = "it containing files of two types of sequence (both DNA and protein).\n" ;
					last; 
				}
			}
		}		
	}
	return $out;
}
