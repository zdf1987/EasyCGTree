#!/usr/bin/perl
use warnings;
use strict; 
no warnings 'experimental::smartmatch';
use File::Copy qw(move mv);
use Getopt::Long;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $rmon=$mon +1;
my $ryear=$year+1900;

print "##############################\nReading command line... checking options and input data...\n\n";
my @usage=qq(
====== EasyCGTree ======
     Version 2.0 by Dao-Feng Zhang (Hohai University, China)
	 Update 2021-06-22

Usage: perl EasyCGTree.pl [Options]

Options:
-input <String>	(Essential)
	Input data (genome/proteome) directory
-query <String>	(Essential)
	Directory of protein sequence data for gene calling
-mode <String, 'A', 'A1', 'A2', 'A3'> (Essential)
	Set run mode
-seq_type <String, 'nucl', 'prot'> (Essential)
	Sequence type used for tree inference
	
-reference <String> (Essential)
	Directory of protein sequence data used in tree inference
-thread <Int> (Optional)		
	Number of threads to be used by 'blast+'. [default: 2]
-blast_dir <String> (Optional)
	Directory for programs of blast+, e.g. /share/bin/
-iden_cutoff <Int, 50..100> (Optional)
	Cutoff(%) for filtering the BLAST results. [default: 50]
-gene_cutoff <Decimal, 0.5..1> (Optional)
	Cutoff for omitting low-prevalence gene. [default: 0.8]
-genome_cutoff <Decimal, 0.5..1> (Optional)
	Cutoff for omitting low-quality genomes. [default: 0.8]
	
-help (Optional)				
	Display this message

);
my ($inputDir,$queryDir,$mode,$seqtype,$Reference,$blast,$check1,$check2,$check3,@query,@queryR);

my %opt=qw();
GetOptions(\%opt,"input:s","query:s","mode:s","seq_type:s","reference:s","thread:i","blast_dir:s","iden_cutoff:f","gene_cutoff:f","genome_cutoff:f","help!");


if (scalar(keys %opt )==0 || exists($opt{"help"})) {
	print join("\n",@usage)."\n\n";
	exit;
}

unless (exists($opt{"mode"})+exists($opt{"seq_type"})+exists($opt{"query"})+exists($opt{"input"})==4) {
	print join("\n",@usage)."##############\n\nERROR: '-input', '-query', '-seq_type' and '-mode' must be specified to start.\n\n";
	exit;
		
}	

{
	$inputDir = $opt{"input"};
	unless ($inputDir=~/\w/ || $inputDir=~/\d/) {
		print join("\n",@usage)."##############\n\nERROR: Argument '-input': it should be a folder name, e.g., myGenome.\n\n";
		exit;
	}

}

{
	$queryDir = $opt{"query"};
	$check1= &checkDataType($queryDir);
	if ($check1 eq "nucl") {
		print join("\n",@usage)."##############\n\nERROR: Argument '-query $queryDir': directory '$queryDir' should contain file(s) of protein sequence, while only file(s) of nucleotide sequence were found.\n\n";
		exit;
	} elsif ($check1 eq "prot") {
	} elsif ($check1=~ /\.fas/ || $check1=~ /open input/) {
		print join("\n",@usage)."##############\n\nERROR: Argument '-query $queryDir': $check1\n\n";
		exit;
	} else {
		print join("\n",@usage)."##############\n\nERROR: Argument '-query $queryDir': directory '$queryDir' should only contain fasta-formated file(s) of protein sequence, while $check1\n\n";
		exit;		
	}
}


if (exists($opt{"reference"})) {
	$Reference = $opt{"reference"};
	$check3= &checkDataType($Reference);
	if ($check3 eq "nucl") {
		print join("\n",@usage)."##############\n\nERROR: '-reference $Reference' only allows the directory '$Reference' containing file(s) of protein sequence.\n\n";
		exit;
	} elsif ($check1 eq "prot") {
	} elsif ($check1=~ /\.fas/) {
		print join("\n",@usage)."##############\n\nERROR: $check3\n\n";
		exit;
	} else {
		print join("\n",@usage)."##############\n\nERROR: '-reference $Reference' only allows the directory '$Reference' containing fasta-formated file(s) of protein sequence, while $check1\n\n";
		exit;		
	}
	opendir(DIR, $Reference);
	@queryR = readdir(DIR);
	closedir(DIR);
	splice @queryR, 0, 1 if $queryR[0]=~ /\.fas/;
	splice @queryR, 1, 1 if $queryR[1]=~ /\.fas/;
}


$seqtype = $opt{"seq_type"};

{
	$mode = $opt{"mode"};
	if ($mode eq "A") {
		if ($Reference) {
			print join("\n",@usage)."##############\n\nERROR: '-mode A/A1/A2' dose not allow the setting of '-reference '.\n\n";
			exit;
		}
	} elsif ($mode eq "A1") {
		if ($Reference) {
			print join("\n",@usage)."##############\n\nERROR: '-mode A/A1/A2' dose not allow the setting of '-reference '.\n\n";
			exit;
		}
	} elsif ($mode eq "A2") {
		if ($Reference) {
			print join("\n",@usage)."##############\n\nERROR: '-mode A/A1/A2' dose not allow the setting of '-reference '.\n\n";
			exit;
		}
	} elsif ($mode eq "A3") {
		
		unless ($Reference) {
			print join("\n",@usage)."##############\n\nERROR: '-mode A3' needs the setting of '-reference '.\n\n";
			exit;
		}
		unless  ($seqtype eq "prot") {
			print join("\n",@usage)."##############\n\nERROR: '-mode A3' needs the setting of '-seq_type prot'.\n\n" ; 
			exit;
		}
		
	} else {
		print join("\n",@usage)."##############\n\nERROR: Argument '-mode': it should be 'A', 'A1', 'A2' or 'A3'.\n\n"; 

		exit;
	}

} 




{
#####$seqtype = $opt{"seq_type"};
	$check2= &checkDataType($inputDir);
	
	opendir(DIR, $inputDir);
	@query = readdir(DIR);
	closedir(DIR);
	
	splice @query, 0, 1 if $query[0]=~ /\.fas/;
	splice @query, 1, 1 if $query[1]=~ /\.fas/;
	unless ($#query + $#queryR >=3) {
		print join("\n",@usage)."##############\n\nERROR: $#query + $#queryR Number of files in both Input directory '$inputDir' and Reference directory '$Reference' must be >5 to start a run.\n\n";
		exit;
	}
	
	if ($seqtype eq "nucl") {
		$blast = "tblastn";
		if ($check2 eq "nucl") {

		} elsif ($check2 eq "prot") {
			print join("\n",@usage)."##############\n\nERROR: '-seq_type nucl' only allows Input directory '$inputDir' containing file(s) of nucleotide sequence.\n\n";
			exit;
		} elsif ($check2=~ /\.fas/) {
			print join("\n",@usage)."##############\n\nERROR: $check2\n\n";
			exit;
		} else {
			print join("\n",@usage)."##############\n\nERROR: '-seq_type nucl' only allows Input directory '$inputDir' containing fasta-formated file(s) of nucleotide sequence, while $check2\n\n";
			exit;	
		}
	} elsif ($seqtype eq "prot") {
		$blast = "blastp";
		if ($check2 eq "nucl") {
			print join("\n",@usage)."##############\n\nERROR: '-seq_type prot' only allows Input directory '$inputDir' containing file(s) of protein sequence.\n\n";
			exit;
		} elsif ($check2 eq "prot") {
		} elsif ($check2=~ /\.fas/) {
			print join("\n",@usage)."##############\n\nERROR: $check2\n\n";
			exit;
		} else {
			print join("\n",@usage)."##############\n\nERROR: '-seq_type prot' only allows Input directory '$inputDir' containing fasta-formated file(s) of protein sequence, while $check2\n\n";
			exit;		
		}
	} else {
		print join("\n",@usage)."##############\n\nERROR: Argument '-seq_type': it should be 'nucl' or 'prot'.\n\n";
		exit;
	}
}



########## Default parameters ################
my $tblastn=" ";
my $blastIdentitycutoff=50;
my $geneCutoff=0.8;
my $genomeCutoff=0.8;
my $num_threads=2;
##########################################

open (LOG, ">$inputDir\_$ryear-$rmon-$mday-$hour-$min.log") or die "Can't open '$inputDir\_$ryear-$rmon-$mday-$hour-$min.log': $!\n";

print "\n\n\n###############Starting Informations ###############\n\n";
print LOG "\n\n\n###############Starting Informations ###############\n\n";

my $qwer;
$tblastn=$opt{"blast_dir"} if exists($opt{"blast_dir"});
if (exists($opt{"iden_cutoff"})) {
	$qwer=$opt{"iden_cutoff"};
	if ($opt{"iden_cutoff"} =~ /^\d+$/) {
		if ($opt{"iden_cutoff"}>=50 || 100>=$opt{"iden_cutoff"}) {
			$blastIdentitycutoff=$opt{"iden_cutoff"};
		} else {
			print join("\n",@usage)."##############\n\nWARNING: Argument '-iden_cutoff': the value '$qwer' is out of recommended interval 50-100.\nThe default value '50' will be used instead.\n";##$opt{"iden_cutoff"}
			print LOG join("\n",@usage)."##############\n\nWARNING: Argument '-iden_cutoff': the value '$qwer' is out of recommended interval 50-100.\nThe default value '50' will be used instead.\n";
		}
	} else {
		print join("\n",@usage)."##############\n\nWARNING: Argument '-iden_cutoff': '$qwer' is not an integer.\nThe default value '50' will be used instead.\n";
		print LOG join("\n",@usage)."##############\n\nWARNING: Argument '-iden_cutoff': '$qwer' is not an integer.\nThe default value '50' will be used instead.\n";
	}
}	
my $rewq;
if (exists($opt{"gene_cutoff"}))	 {
	$rewq=$opt{"gene_cutoff"};
	if ($opt{"gene_cutoff"} =~ /^0?\.?\d+$/) {
		if ($opt{"gene_cutoff"}>=0.5 || 1>=$opt{"gene_cutoff"}) {
			$geneCutoff= $opt{"gene_cutoff"};
		} else {
			print join("\n",@usage)."##############\n\nWARNING: Argument '-gene_cutoff': the value '$rewq' is out of recommended interval 0.5-1.\nThe default value '0.8' will be used instead.\n";
			print LOG join("\n",@usage)."##############\n\nWARNING: Argument '-gene_cutoff': the value '$rewq' is out of recommended interval 0.5-1.\nThe default value '0.8' will be used instead.\n";
		}
	} else {
		print join("\n",@usage)."##############\n\nWARNING: Argument '-gene_cutoff': '$rewq' is not a decimal.\nThe default value '0.8' will be used instead.\n";
		print LOG join("\n",@usage)."##############\n\nWARNING: Argument '-gene_cutoff': '$rewq' is not a decimal.\nThe default value '0.8' will be used instead.\n";
	}
}
my $uyt;
if (exists($opt{"genome_cutoff"})) {
	$uyt=$opt{"genome_cutoff"};
	if ($opt{"genome_cutoff"} =~ /^0?\.?\d+$/) {
		if ($opt{"genome_cutoff"}>=0.5 || 1>=$opt{"genome_cutoff"}) {
			$genomeCutoff= $opt{"genome_cutoff"};
		} else {
			print join("\n",@usage)."##############\n\nWARNING: Argument '-ggenome_cutoff': the value '$uyt' is out of recommended interval 0.5-1.\nThe default value '0.8' will be used instead.\n";
			print LOG join("\n",@usage)."##############\n\nWARNING: Argument '-genome_cutoff': the value '$uyt' is out of recommended interval 0.5-1.\nThe default value '0.8' will be used instead.\n";
		}
	} else {
		print join("\n",@usage)."##############\n\nWARNING: Argument '-genome_cutoff': '$uyt' is not a decimal.\nThe default value '0.8' will be used instead.\n";
		print LOG join("\n",@usage)."##############\n\nWARNING: Argument '-genome_cutoff': '$uyt' is not a decimal.\nThe default value '0.8' will be used instead.\n";
	}
}	
my $BNM;
if (exists($opt{"thread"})) {	
	$BNM=$opt{"thread"};
	if ($opt{"thread"} =~ /^\d+$/) {
			$num_threads=$opt{"thread"};
		
	} else {
		print join("\n",@usage)."##############\n\nWARNING: Argument '-thread': '$BNM' is not an integer.\nThe default value '2' will be used instead.\n";
		print LOG join("\n",@usage)."##############\n\nWARNING: Argument '-thread': '$BNM' is not an integer.\nThe default value '2' will be used instead.\n";
	}
}


print LOG '====== EasyCGTree ======
	Version 2.0
by Dao-Feng Zhang (Hohai University, China)
';
print '====== EasyCGTree ======
	Version 2.0
by Dao-Feng Zhang (Hohai University, China)
';


print LOG "#########################\n\nOptions: \n";  
print "#########################\n\nOptions: \n";

foreach my $kk (keys %opt) {
	print LOG "-$kk $opt{$kk}\n";  
	print "-$kk $opt{$kk}\n";
}
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

	mkdir "$TEMdir" || die " Permission denied to create directories. Please contact the Administrator of your System\n"; ##########21.06.19
	
	
	opendir(DIR, $inputDir) || die "Step 1: Can't open input directory '$inputDir': $!\n";
	my @query = readdir(DIR);
	closedir(DIR);
	die "Step 1: Input directory $inputDir does not contain any files.\n" unless scalar(@query);
	mkdir "$TEMdir/TEM1_blastDB";##########21.06.19
	my $genomeNum = 0;
	
	foreach my $file (@query) {
		next unless $file =~ /\.fas/;
		$genomeNum++;
		my $cmd;
		if ($tblastn) {
			$cmd = "$tblastn"."makeblastdb -in ./$inputDir/$file -dbtype $seqtype -out ./$TEMdir/TEM1_blastDB/$file";
		} else {
			$cmd = "makeblastdb -in ./$inputDir/$file -dbtype $seqtype -out ./$TEMdir/TEM1_blastDB/$file";
		}
		my $response1= `$cmd` || die "Step 1: Can't execute 'makeblastdb': $!\nPlease make sure: blast+ package was added to the environment variable of yous system; or the path was set correctly in 'parameters.txt'.\n";  ##########21.06.19

		my @res1=split /\n/, $response1;##########21.06.19
		foreach my $res1 (@res1) {##########21.06.19
			
			next unless $res1 =~ /error/i;##########21.06.19
			print LOG "Step 1 ERROR: Something goes wrong when execute 'makeblastdb':\n$res1\n\n\n\nEasyCGTree stopped because of problematic input data: ./$inputDir/$file\n";##########21.06.19
			die  "Step 1 ERROR: Something goes wrong when execute 'makeblastdb':\n$res1\n\n\n\nEasyCGTree stopped because of problematic input data: ./$inputDir/$file\n";##########21.06.19
			}
		
		print "Step 1: Making BLASTDB: It's the $genomeNum task of directory '$inputDir'!\n";
		print LOG "Step 1: Making BLASTDB: It's the $genomeNum task of directory '$inputDir'!\n";
	}
	my $refN=0;
	my @queryR;
	if ($mode eq "A3") {
		opendir(DIR, "Reference") || die "Step 1: Can't open input directory 'Reference': $!\n";
		@queryR = readdir(DIR);
		closedir(DIR);
		die "Step 1: Input directory 'Reference' does not contain any files.\n" unless scalar(@queryR);
		foreach my $file (@queryR) {
			next unless $file =~ /\.fas/;
			$genomeNum++;
			$refN++;
			my $cmd;
			if ($tblastn) {
				$cmd = "$tblastn"."makeblastdb -in ./$inputDir/$file -dbtype $seqtype -out ./$TEMdir/TEM1_blastDB/$file";
			} else {
				$cmd = "makeblastdb -in ./$inputDir/$file -dbtype $seqtype -out ./$TEMdir/TEM1_blastDB/$file";
			}
			my $response1= `$cmd` || die "Step 1: Can't execute 'makeblastdb': $!\nPlease check whether blast+ package was installed correctly!\n";  ##########21.06.19
			my @res1=split /\n/, $response1;##########21.06.19
			foreach my $res1 (@res1) {##########21.06.19
				next unless $res1 =~ /error/i;##########21.06.19
				print LOG "Step 1 ERROR: Something goes wrong when execute 'makeblastdb':\n$res1\n\n\n\nEasyCGTree stopped because of problematic input data: ./Reference/$file\n";##########21.06.19
				die  "Step 1 ERROR: Something goes wrong when execute 'makeblastdb':\n$res1\n\n\n\nEasyCGTree stopped because of problematic input data: ./Reference/$file\n";##########21.06.19
			}
			print "Step 1: Making BLASTDB: It's the $refN task of directory 'Reference'!\n";
			print LOG "Step 1: Making BLASTDB: It's the $refN task of directory 'Reference'!\n";
		}
	print "Step 1: ##############Totally, $genomeNum BLAST Database were created!##############\n";
	print LOG "Step 1: ##############Totally, $genomeNum BLAST Database were created!##############\n";
	}

	
############### Step 2: Blast Search #################
print "\n\n\n############### Step 2: Blast Search ###############\n\n";
print LOG "\n\n\n############### Step 2: Blast Search ###############\n\n";
	opendir(DIR2, $queryDir) || die "Step 1: Can't open input directory '$queryDir': $!\n";
	my @query2 = readdir(DIR2);
	closedir(DIR2);
	die "Step 2: Input Query directory $queryDir does not contain any files.\n" unless scalar(@query2);
	my $num2=0;
	open (QUERY, ">$queryDir.seq") or die "Step 1: Can't open input file '$queryDir.seq': $!\n";
	foreach my $file2 (@query2) {
		next unless $file2 =~ /\.fas/;
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
	print "Step 2: The query file $queryDir.seq have been created!\n\n";
	print LOG "Step 2: The query file $queryDir.seq have been created!\n\n";
	mkdir "$TEMdir/TEM2_blastOUT";############################################21.06.19
	
	my $query = "$queryDir.seq"; 

	my $number = 0;
	foreach my $file (@query) {
		next unless $file =~ /\.fas/;
		my @in = split (/_/, $file);
		my $cmd;
		if ($tblastn) {
			$cmd="$tblastn"."$blast -query $query -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
		} else {
			$cmd="$blast -query $query -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
		}

		my $response2 = `$cmd`;
		print LOG "Step 2 ERROR: Something goes wrong when execute BLAST search:$response2\nPlease check the files in folder '$queryDir'. They should be fasta-formated amino acid sequence and not empty!\n\n\n\nEasyCGTree stopped because of problematic input file(s) in '$queryDir'\n" if $response2 =~ /empty/i || $response2 =~ /warning/i;##########21.06.19
		die "Step 2 ERROR: Something goes wrong when execute BLAST search:$response2\nPlease check the files in folder '$queryDir'. They should be fasta-formated amino acid sequence and not empty!\n\n\n\nEasyCGTree stopped because of problematic input file(s) in '$queryDir'\n" if $response2 =~ /empty/i || $response2 =~ /warning/i;##########21.06.19

		$number++;
		
		print "Step 2: BLAST: It's the $number/$genomeNum task!\n";
		print LOG "Step 2: BLAST: It's the $number/$genomeNum task!\n";
	}
	my $number0 = 0;
	if ($mode eq "A3")	{
		foreach my $file (@queryR) {
			next unless $file =~ /\.fas/;
			my @in = split (/_/, $file);
			my $cmd;
			if ($tblastn) {
				$cmd="$tblastn"."$blast -query $query -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
			} else {
				$cmd="$blast -query $query -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
			}
			$number++;
			system $cmd || die "Step 2: Can't execute '$tblastn': $!\n";
			$number0++;
			my $total= $#queryR + 1;
			print "Step 2: BLAST: It's the $number0/$refN task of directory 'Reference'!\n";
			print LOG "Step 2: BLAST: It's the $number0/$refN task of directory 'Reference'!\n";
		}
	}
	print "Step 2: ##############Totally, $number BLAST searches were performed!##############\n";
	print LOG "Step 2: ##############Totally, $number BLAST BLAST searches were performed!##############\n";


############### Step 3: Filtrate the Blast Results ################# 符合标准的blast结果都输出，不考虑所有基因组共有情况	
print "\n\n\n############### Step 3: Filtrate Blast Results ###############\n\n";
print LOG "\n\n\n############### Step 3: Filtrate Blast Results ###############\n\n";

	mkdir "$TEMdir/TEM3_blastOUT_S";##########21.06.19
	
	opendir(DIR3, "$TEMdir/TEM2_blastOUT") || die "Step 1: Can't open input directory '$TEMdir/TEM2_blastOUT': $!\n";
	my @query3 = readdir(DIR3);
	closedir(DIR3);
	die "Step 3: Input  directory $TEMdir/TEM2_blastOUT does not contain any files.\n" unless scalar(@query3);

	foreach my $file3 (@query3) {

		my (%out3,%score3);
		next unless $file3 =~ /\.fas/;                    

		open (FIL3, "<$TEMdir/TEM2_blastOUT/$file3") or die "Can't open '$TEMdir/TEM2_blastOUT/$file3': $!\n";
		my (@key3,@result3);

		foreach my $in3 (<FIL3>) {
			$in3=~ s/\n// unless $in3=~s/\r\n//;
			next unless $in3;
		
			my @in3 =split /\s+/, $in3;                           
			next if $in3[4] < $blastIdentitycutoff;               ###############
			my @in33 =split /_/, $in3[0];
			push @key3, $in33[3];									################$in33[3]还是$in33[2]根据query基因标签决定，r89和r95的faa序列标签格式不同。
			unless ($out3{$in33[3]}) {
				$out3{$in33[3]}=$in3;
				$score3{$in33[3]}=$in3[7];
				
			} else {
				$out3{$in33[3]}=$in3 if $score3{$in33[3]}<$in3[7];
				$score3{$in33[3]}=$in3[7] if $score3{$in33[3]}<$in3[7];
			}
		}
		close FIL3;
	
		open(OU3, ">$TEMdir/TEM3_blastOUT_S/$file3") or die "Step 3: Can't open '$TEMdir/TEM3_blastOUT_S/$file3': $!\n";

		foreach my $in (keys %out3) {
			print OU3 "$out3{$in}\n";
		}
		close OU3;
	}
		
	

}

if ($mode eq "A1") {
	print "The user ordered first-part (A1) analysis, and it is finished!\n\n";
	print LOG "The user ordered first-part (A1) analysis, and it is finished!\n\n";

} else {

opendir(DIR, "$TEMdir/TEM3_blastOUT_S") || die "Step 3: Can't open input directory '$TEMdir/TEM3_blastOUT_S': $!\n";
my @query3 = readdir(DIR);
closedir(DIR);
die "Step 3: Input directory TEM3_blastOUT_S does not contain any files.\n" unless scalar(@query3);
my (@querylist0,@querylistN);

my (@filename,@genelist);
my $num=0;
foreach my $file3 (@query3) {
	
	next unless $file3 =~ /.\.fas/;
	push @filename, $file3;
	open (FIL3, "<$TEMdir/TEM3_blastOUT_S/$file3") or die "Step 3: Can't open '$TEMdir/TEM3_blastOUT_S/$file3': $!\n";
     	my $putin;
	
	foreach my $in3 (<FIL3>) {
		$in3=~ s/\n// unless $in3=~s/\r\n//;
		next unless $in3;

		my @in3 =split /\s/, $in3;

		my @in31 =split /_/, $in3[0];
		if ($putin) {
			die "Step 3: This file $file3 contains reduplicative sequences at $in31[3] locus.\n" if $putin =~ /$in31[3]	/ ;
			$putin .="$in31[3]	";################$in31[3]还是$in31[2]根据query基因标签决定，r89和r95的faa序列标签格式不同。
			
		} else {
			$putin ="$in31[3]	";
			
		}
		unless (@querylist0) {
			push @querylist0,$in31[3];
			$querylistN[$num]++;
			$num++;
			next;
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

	push @genelist, $putin;
	close FIL;
}


my (@finalgenelist,@outgenelist,@finalgenome,@outgenome);
my $cutoff=int(($#filename+1)*$geneCutoff);           ###############
foreach my $kk (0..$#querylistN) {
	push @outgenelist, $querylist0[$kk] if $querylistN[$kk]<  $cutoff;
	push @finalgenelist, $querylist0[$kk] if $querylistN[$kk]>=  $cutoff;
}
my $cutoff2=int(($#finalgenelist+1)*$genomeCutoff);     ####################
foreach my $kk (0..$#filename) {
	my @yy=split /\t/,$genelist[$kk];
	my $nummv=0;
	foreach my $vv (@yy) {
		$nummv++ if $vv~~@finalgenelist;
	}
	
	push @finalgenome, $filename[$kk] if $nummv>=  $cutoff2;
	push @outgenome, $filename[$kk] if $nummv <  $cutoff2;
}
my $ingene =$#finalgenelist+1;
my $outgene =$#outgenelist+1;
my $cmd33="$TEMdir/TEM3_blastOUT_S/Discard";##########################
mkdir "$cmd33" || die "Step 3: Can't create blastOUT directory '$TEMdir/TEM3_blastOUT_S/Discard'!\n";
my $dis=0;
foreach my $kk (@outgenome) {
	move ("./$TEMdir/TEM3_blastOUT_S/$kk","./$TEMdir/TEM3_blastOUT_S/Discard/$kk")|| die "Step 3: Can't move './$TEMdir/TEM3_blastOUT_S/$kk' into './$TEMdir/TEM3_blastOUT_S/Discard/$kk'!\n";
	
	$dis++;
}
my $ingenomeN=$#finalgenome+1;

my $genomeNum=$#filename+1;
print "Step 3: ========\n There are $ingene/$num common genes present in more than $cutoff ($genomeNum * $geneCutoff [optional setting]) genomes, will be used in following analysis. This is the list:\n@finalgenelist;\n";
print LOG "Step 3: ========\n There are $ingene/$num common genes present in more than $cutoff ($genomeNum * $geneCutoff [optional setting]) genomes, will be used in following analysis. This is the list:\n@finalgenelist;\n";
if ($outgene) {
print "\n$outgene/$num genes were excluded.  This is the list:\n@outgenelist\n";
print LOG "\n$outgene/$num genes were excluded.  This is the list:\n@outgenelist\n";
}

print "========\n$ingenomeN genomes harboring >= $cutoff2 ($ingene * $genomeCutoff [optional setting]) of the $ingene common genes, will be used in following analysis.This is the list:\n@finalgenome.\n";
print LOG "========\n$ingenomeN genomes harboring >= $cutoff2 ($ingene * $genomeCutoff [optional setting]) of the $ingene common genes, will be used in following analysis.This is the list:\n@finalgenome.\n";

if ($dis) {
print "\n$dis genomes were excluded and moved into $TEMdir/TEM3_blastOUT_S/Discard/. This is the list:\n:@outgenome.\n";
print LOG "\n$dis genomes were excluded and moved into $TEMdir/TEM3_blastOUT_S/Discard/. This is the list:\n:@outgenome.\n";
}


############### Step 4: Retrieve Sequences from Each Genome/proteome ################# 提取所有符合标准的基因序列
print "\n\n\n############### Step 4: Retrieve Sequences from Each Genome/proteome ###############\n\n";
print LOG "\n\n\n############### Step 4: Retrieve Sequences from Each Genome/proteome ###############\n\n";

my $numg=0;
mkdir "$TEMdir/TEM4_GeneSeqs";##########21.06.19

opendir(DIR, "$TEMdir/TEM3_blastOUT_S") || die "Step 4: Can't open input directory '$TEMdir/TEM3_blastOUT_S':$!\n";
my @query4 = readdir(DIR);
closedir(DIR);
die "Step 4: Input directory $TEMdir/TEM3_blastOUT_S does not contain any files.\n" unless scalar(@query4);

print "Step 4: According to the Step 3, $ingenomeN genomes harboring >= $cutoff2 ($ingene * $genomeCutoff [optional setting]) of the $ingene common genes, will be used in following analysis. The related gene sequences will be extracted.\n";
print LOG "Step 4: According to the Step 3, $ingenomeN genomes harboring >= $cutoff2 ($ingene * $genomeCutoff [optional setting]) of the $ingene common genes, will be used in following analysis. The related gene sequences will be extracted.\n";

foreach my $i1 (@query4) {
	next unless $i1=~ /.\.fas/;
	next unless $i1 ~~ @finalgenome;	##################21.02.24
	$numg++;
	if ($seqtype eq "prot") {    ##################21.02.24
	

	
	my %geno;
	my ($seq,$seqID,$inname);
	$inname = $i1;
	if ($i1=~ /Reference_/) {
	$inname=~ s/Reference_//;
		
		open (FILE, "<Reference/$inname") or die "Step 4: Can't open 'Reference/$inname': $!\n";
	} else {
		open (FILE, "<$inputDir/$i1") or die "Step 4: Can't open '$inputDir/$i1': $!\n";########
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
	
	open(OUT, ">$TEMdir/TEM4_GeneSeqs/$inname") or die "Step 4: Can't open '$TEMdir/TEM4_GeneSeqs/$i1': $!\n";
	open (FIL, "<$TEMdir/TEM3_blastOUT_S/$i1") or die "Step 4: Can't open 'TEM3_blastOUT_S/$i1': $!\n";
	my $seqTag="1_1_$inname";		
	$seqTag=~ s/\.fas/_/;		
	foreach my $i2 (<FIL>) {
		$i2=~ s/\n// unless $i2=~s/\r\n//;
		next unless $i2;
		my @in0 = split /\s/, $i2;
		my@ff =split /_/, $in0[0];

				
		my $geneSeq = $geno{$in0[1]};	
						
		print OUT ">$seqTag"."$ff[3]\n$geneSeq\n";################$ff[3]还是$ff[2]根据query基因标签决定，r89和r95的faa序列标签格式不同。##################21.06.20
		$seqnum++;
	}
	close FIL;	
			
	close OUT;
	print "Step 4: I got $seqnum/$num sequences from genome $numg/$ingenomeN:$i1!\n";
	print LOG "Step 4: I got $seqnum/$num sequences from genome $numg/$ingenomeN:$i1!\n";
	
	} else {
	
	
	my %geno;
	my ($seq,$seqID);
	open (FILE, "<$inputDir/$i1") or die "Step 4: Can't open '$inputDir/$i1': $!\n";########
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
	
	open(OUT, ">$TEMdir/TEM4_GeneSeqs/$i1") or die "Step 4: Can't open '$TEMdir/TEM4_GeneSeqs/$i1': $!\n";
	open (FIL, "<$TEMdir/TEM3_blastOUT_S/$i1") or die "Step 4: Can't open 'TEM3_blastOUT_S/$i1': $!\n";
	my $seqTag="1_1_$i1";		
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
		print OUT ">$seqTag"."$ff[3]\n$geneSeq\n";################$ff[3]还是$ff[2]根据query基因标签决定，r89和r95的faa序列标签格式不同。
		$seqnum++;
	}
	close FIL;	
			
	close OUT;
	print "Step 4: I got $seqnum/$num sequences for genome $numg/$ingenomeN:$i1!\n";
	print LOG "Step 4: I got $seqnum/$num sequences for genome $numg/$ingenomeN:$i1!\n";
	
	
	}
	
	
	
	
}
		


############### Step 5: Gather the Orthologs in Single File #################
print "\n\n\n############### Step 5: Gather the Orthologs in Single File ###############\n\n";
print LOG "\n\n\n############### Step 5: Gather the Orthologs in Single File ###############\n\n";

opendir(DIR1, "$TEMdir/TEM4_GeneSeqs") || die "Step 5: Can't open input directory '$TEMdir/TEM4_GeneSeqs': $!\n";
my @in1 = readdir(DIR1);
closedir(DIR1);
die "Step 5: Input directory '$TEMdir/TEM4_GeneSeqs' does not contain any files.\n" unless scalar(@in1);

mkdir "$TEMdir/TEM5_GeneCluster";###############21.06.19

my %outClu;
my $num5=0;
print "Step 5: I am handling the";
print LOG "Step 5: I am handling the";
foreach my $in5 (@in1) {
	
	next unless $in5=~ /.\.fas/;
	my $ad0id="Reference_"."$in5";
	next unless $in5 ~~ @finalgenome || $ad0id ~~ @finalgenome;	##################21.05.17
	
	
	my $tt =$in5;
	$tt=~ s/\.fas//;
	open (FILE, "<$TEMdir/TEM4_GeneSeqs/$in5") or die "Step 5: Can't open '$TEMdir/TEM4_GeneSeqs/$in5': $!\n";
	my ($id,$sig);
	foreach my $i1 (<FILE>) {
	
		$i1=~ s/\n// unless $i1=~s/\r\n//;
		next unless $i1;
		if ( $i1 =~ /\>/ ) {
			my @seqids=split /_/, $i1;##################21.06.20
			$id =$seqids[3];##################21.06.20
			
			
			my $adid="Reference_"."$id"; #######################注意这里有改动
			$sig=1 if $id ~~ @finalgenelist;     #######################注意这里有改动
		} elsif ($sig) {
			unless ($outClu{$id}) {
				$outClu{$id}= ">$id\_$tt\n$i1\n";
			} else {
				$outClu{$id} .= ">$id\_$tt\n$i1\n"; ###
			}
			$sig=0;
		}
	}
	close FILE;
	$num5++;
	print " $num5/$ingenomeN";
	print LOG " $num5/$ingenomeN";
}

my $n=0;
print " files.\nStep 5: I am handling the files (file no./sequence no.):";
print LOG " files.\nStep 5: I am handling the files (file no./sequence no.):";
foreach my $in (keys %outClu) {
	my $cnt;
	$n++;
	
	open (OUT, ">./$TEMdir/TEM5_GeneCluster/$in.fas") or die "Step 5: Can't open '$TEMdir/TEM5_GeneCluster/$in.fas': $!\n";
	print OUT "$outClu{$in}";
	$outClu{$in} =~ s/>/"X".$cnt++/ge;
	print " $n/$cnt";
	print LOG " $n/$cnt";
	close OUT;
}
print ".\n";	
print LOG ".\n";


############### Step 6: Deduplicate Sequences #################
print "\n\n\n############### Step 6: Deduplicate Sequences ###############\n\n";
print LOG "\n\n\n############### Step 6: Deduplicate Sequences ###############\n\n";

mkdir "$TEMdir/TEM6_DedupGeneCluster";##########21.06.19

mkdir "$TEMdir/TEM60_DedupList";##########21.06.19

opendir(DIR1, "$TEMdir/TEM5_GeneCluster") || die "Step 6: Can't open input directory '$TEMdir/TEM5_GeneCluster': $!\n";
my @in6 = readdir(DIR1);
closedir(DIR1);
die "Step 6: Input directory $TEMdir/TEM5_GeneCluster does not contain any files.\n" unless scalar(@in6);
print "Step 6: I am handling the files (file no./sequence types):";
print LOG "Step 6: I am handling the files (file no./sequence types):";
my %outDedup;
my $num6=0;
foreach my $in (@in6) {
	my $genome=$in;
	next unless $genome=~ s/\.fas//;
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
	print LOG " $num6/$loop";
}
print ".\n";
print LOG ".\n";


############### Step 7: Do Alignment #################
print "\n\n\n############### Step 7: Do Alignment ###############\n\n"; 
print LOG "\n\n\n############### Step 7: Do Alignment ###############\n\n";
mkdir "$TEMdir/TEM7_Alignment";##########21.06.19


opendir(DIR1, "$TEMdir/TEM6_DedupGeneCluster") || die "Step 7: Can't open input directory '$TEMdir/TEM6_DedupGeneCluster': $!.\n";
my @in7 = readdir(DIR1);
closedir(DIR1);
die "Step 7: Input directory $TEMdir/TEM6_DedupGeneCluster does not contain any files.\n" unless scalar(@in7);


my $num7 = 0;

foreach (@in7) {
	next unless /\w/;
	next unless /\.fas/;   #######21.02.24
	my $program = "clustalo";##########################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	my $outname = $_ . ".fasta";
	my @para = ("-i", "./$TEMdir/TEM6_DedupGeneCluster/$_", "-o", "./$TEMdir/TEM7_Alignment/$outname", "--outfmt=fasta", "--output-order=tree-order","--threads=$num_threads","-v");##########################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	unless (system './bin/clustalo', @para) {##########################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		$num7++;
		print "\nStep 7: $_ has been completed! It is the $num7/$ingene file.\n\n";
		print LOG "\nStep 7: $_ has been completed! It is the $num7/$ingene file.\n\n";
	} else {
		die "Step 7: Can't execute './bin/clustalo': $!.\n";   ##########################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	}
}


############### Step 8: Trim Sequence Algnments #################
print "\n\n\n############### Step 8: Trim Sequence Algnments ###############\n\n";
print LOG "\n\n\n############### Step 8: Trim Sequence Algnments ###############\n\n";

mkdir "$TEMdir/TEM8_AlnTrimmed";##########21.06.19


opendir(DIR, "$TEMdir/TEM7_Alignment") || die "Step 8: Can't open input directory '$TEMdir/TEM7_Alignment':$!\n";
my @in8 = readdir(DIR);
closedir(DIR);
die "Step 8: Input directory $TEMdir/TEM7_Alignment does not contain any files.\n" unless scalar(@in8);
my (@problem,@error);
my $num8 = 0;
print "Step 8: I am handling the files:";
print LOG "Step 8: I am handling the files:";
foreach my $file (@in8) {
	next unless $file =~ /\.fasta$/;
  
	(my $outname = $file) =~ s/\.fas//;
	open (FILE, "$TEMdir/TEM7_Alignment/$file") or die "Step 8: Can't open '$file': $!\n";
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
	print LOG " $num8";
	my $ratio =($end_pos - $fore_pos + 1)/$length;
	push @problem, "$file $ratio\n" if $ratio < 0.5;
	push @error, "$file $ratio\n" if $ratio < 0;
}
my $pro=$#problem +1;
my $err=$#error +1;
print ".\nStep 8: These $pro cluster(s) are problemic (less than half length was retained for following analysis), please check related file(s) in 'TEM7_Alignment':\n*****Indicated by the filenames and the coverages*****\n@problem\n";
print LOG ".\nStep 8: These $pro cluster(s) are problemic (less than half length was retained for following analysis), please check related file(s) in 'TEM7_Alignment':\n*****Indicated by the filenames and the coverages*****\n@problem\n";

if (@error) {
print "These $err +1 clusters are ERROR:\n@error\n";
print LOG "These $err +1 clusters are ERROR:\n@error\n";
}


############### Step 9: ReCall Sequences #################
print "\n\n\n############### Step 9: ReCall Sequences ###############\n\n";
print LOG "\n\n\n############### Step 9: ReCall Sequences ###############\n\n";

mkdir "$TEMdir/TEM9_ReCallAln";##########21.06.19


opendir(DIR1, "$TEMdir/TEM8_AlnTrimmed") || die "Step 9: Can't open input directory '$TEMdir/TTEM8_AlnTrimmed': $!\n";
my @in9 = readdir(DIR1);
closedir(DIR1);
die "Step 9: Input directory $TEMdir/TEM7_Alignment does not contain any files.\n" unless scalar(@in9);
print "Step 9: I am handling the files(file no./sequence no.):";
print LOG "Step 9: I am handling the files(file no./sequence no.):";
my $num9=0;
foreach my $in (@in9) {
	my $genome=$in;
	next unless $genome=~ s/\.fasta//;
	open (FILE, "<$TEMdir/TEM8_AlnTrimmed/$in") or die "Step 9: Can't open '$TEMdir/TEM8_AlnTrimmed/$in': $!\n";
	
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
	die "Step 9: There must be something wrong in file $in!\n" unless $#id==$#seq;
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
					$ll=~ s/$genome\_//;
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
		$i9 =~ s/Reference_//;
		next if $i9~~@current;
		print OU ">$genome\_$i9\n$seq9";
	}


	close OU;	
	
	print " $num9/$num99";
	print LOG " $num9/$num99";
}
print ".\n";
print LOG ".\n";


############### Step 10: Concatenate Sequences #################
print "\n\n\n############### Step 10: Concatenate Sequences ###############\n\n";
print LOG "\n\n\n############### Step 10: Concatenate Sequences ###############\n\n";

opendir(DIR, "$TEMdir/TEM9_ReCallAln") || die "Step 10: Can't open input directory '$TEMdir/TEM9_ReCallAln': $!\n";
my @in10 = readdir(DIR);
closedir(DIR);
die "Step 10: Input directory $TEMdir/TEM9_ReCallAln does not contain any files.\n" unless scalar(@in10);

print "Step 10: I am handling the files:";
print LOG "Step 10: I am handling the files:";
my $num10 = 0;
my (%out,@id,@seq,$id,$seq);
foreach my $file ( sort @in10) {
	next unless $file =~ /\.fasta$/;
  	$file =~ /\.fasta$/ || die "Step 10: File '$file' does not have a name in xxxx.fasta format.\n";;
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
				my @kk =split /_/,$i1;
				push @id, $kk[1];
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
				my @kk =split /_/,$i1;
				$id=$kk[1];
			} else {
				$out{$id} .="$i1";
			}
		}
	}
 	$num10++;
	
	print " $num10";
        print LOG " $num10";

}
print ".\n";
print LOG ".\n";
open (OUT, ">$inputDir.concatenation.fas") or die "Can't open '$inputDir.concatenation.fas': $!\n";
my $nn=0;
foreach (sort keys %out) {

	$nn++;
	
	print OUT ">$_\n";
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

print "Step 10: I got $nn concatenated sequences!\n";
print LOG "Step 10: I got $nn concatenated sequences!\n";


############### Step 11: Construct a Tree #################
print "\n\n\n############### Step 11: Construct a Tree ###############\n\n";
print LOG "\n\n\n############### Step 11: Construct a Tree ###############\n\n";

print "Step 11: I am making the tree......";
print LOG "Step 10: I am making the tree......";
my $cmd11;
if ($seqtype eq "nucl") {
	$cmd11="./bin/FastTreeMP -gtr -nt $inputDir.concatenation.fas > $inputDir.tree";##########################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
} elsif ($seqtype eq "prot") {
	$cmd11="./bin/FastTreeMP $inputDir.concatenation.fas > $inputDir.tree";	##########################%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
}

my $response11= `$cmd11`;  ##########21.06.19
my @res11=split /\n/, $response11;##########21.06.19
print LOG "$response11\n";
print "$response11\n";
foreach my $res1 (@res11) {##########21.06.19
	if  ($res1 =~ /truncated/i && $res1 =~ /be too long/i) {##########21.06.19
		print LOG "Step 11 ERROR: EasyCGTree stopped when execute 'FastTreeMP':\n$response11\n\n\n\nPlease check whether the genome/proteome files were formated correctly. \nPlease find more informations in section 3.1 of the Manual.\n";##########21.06.19
		die  "Step 11 ERROR: EasyCGTree stopped when execute 'FastTreeMP':\n$response11\n\n\n\nPlease check whether the genome/proteome files were formated correctly. \nPlease find more informations in section 3.1 of the Manual.\n\n";##########21.06.19
	} elsif ($res1 =~ /out of memory/i) {
		print LOG "Step 11 ERROR: EasyCGTree stopped when execute 'FastTreeMP':\n$response11\n\n\n\nThis ERROR occured probably because 'TastTreeMP' needs larger memory space than that of your computer. Please find a more powerful PC or server to complete the run.\n";##########21.06.19
		die  "Step 11 ERROR: EasyCGTree stopped when execute 'FastTreeMP':\n$response11\n\n\n\nThis ERROR occured probably because 'TastTreeMP' needs larger memory space than that of your computer. Please find a more powerful PC or server to complete the run.\n";##########21.06.19
	}
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
	my @query = readdir(DIR);
	closedir(DIR);
	
	my @name;
	foreach my $in00 (@query) {
	
		next unless $in00 =~ /\w/;	
		unless ($in00 =~ /\.fas/) {
			$out="it contain file(s) of which the name ends out of '.fas'. Please formate them by using script 'FormatNames.pl'.\n";
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