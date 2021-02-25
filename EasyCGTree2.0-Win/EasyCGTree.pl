#!/usr/bin/perl
use warnings;
use strict; 
my @usage=qq(
====== EasyCGTree ======
     Version 2.0 by Dao-Feng Zhang (Hohai University, China)
	 Update 2020-12-22

Usage: perl EasyCGTree.pl genomeFolder queryFolder A/A1/A2/A3 nucl/prot parameters.txt(Optional)

*************** Set in parameters.txt (Optional) ***************

tblastn=					specify the location of the program tblastn, e.g. /share/bin/;  

num_threads=2				specify the number of threads used by the program tblastn; 

blastIdentitycutoff=50		specify the cutoff for filtering the tblastn results;

geneCutoff=0.8				specify the cutoff for omitting low-prevalence gene;

genomeCutoff=0.8			specify the cutoff for omitting low-quality genomes;
);

print join("\n",@usage)."\n" unless $ARGV[3];       
die unless $ARGV[3];
my @querylist;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $inputDir = $ARGV[0];
my $rmon=$mon +1;
my $ryear=$year+1900;
open (LOG, ">$inputDir\_$ryear-$rmon-$mday-$hour-$min.log") or die "Can't open '$inputDir\_$ryear-$rmon-$mday-$hour-$min.log': $!\n";
print LOG '====== EasyCGTree ======
Version 2.0
by Dao-Feng Zhang (Hohai University, China)
';
print '====== EasyCGTree ======
Version 2.0
by Dao-Feng Zhang (Hohai University, China)
';

print "\n\n\n###############Starting Informations ###############\n\n";
print LOG "\n\n\n###############Starting Informations ###############\n\n";

print LOG "Command:  perl EasyCGTree.pl @ARGV\n\nJob Started at: $hour:$min:$sec,$mday-$rmon-$ryear\n\n";  
print "Command:  perl EasyCGTree.pl @ARGV\n\nJob Started at: $hour:$min:$sec,$mday-$rmon-$ryear\n\n";

if ($ARGV[2] eq "A") {
	print "The user ordered a complete analysis (type A) of the pipeline!\n\n";
	print LOG "The user ordered a complete analysis (type A) of the pipeline!\n\n";
} elsif ($ARGV[2] eq "A1") {
	print "The user ordered a first-part analysis (type A1) of the pipeline!\n\n";
	print LOG "The user ordered a first-part analysis (type A1) of the pipeline!\n\n";
} elsif ($ARGV[2] eq "A2") {
	print "The user ordered a second-part analysis (type A2) of the pipeline!\n\n";
	print LOG "The user ordered a second-part analysis (type A2) of the pipeline!\n\n";
} elsif ($ARGV[2] eq "A3") {
	print "The user ordered a complete analysis (type A3) with a pre-defined core-gene set of some genomes, which will be also used in current analysis!\n\n";
	print LOG "The user ordered a complete analysis (type A3) with a pre-defined core-gene set of some genomes, which will be also used in current analysis!\n\n";
	die "ERROR: Parameter 'A3' is conflicting to parameter 'nucl'. !\n\n" if $ARGV[3] eq "nucl";
	
} else {
	print join("\n",@usage)."\n"; 
	print LOG "WARNING: Problematic command line!\n\n";
	die "WARNING: Problematic command line!\n\n";
}
my ($dbtype,$blast);
if ($ARGV[3] eq "nucl") {
	$blast = "tblastn";
	print "NOTE: The input directory $ARGV[1] must include genomes of fasta-formated!\n\n";
	print LOG "NOTE: The input directory $ARGV[1] must include genomes of fasta-formated!\n\n";
} elsif ($ARGV[3] eq "prot") {
	$blast = "blastp";
	print "NOTE: The input directory $ARGV[1] must include proteomes of fasta-formated!\n\n";
	print LOG "NOTE: The input directory $ARGV[1] must include proteomes of fasta-formated!\n\n";
} else {
	print join("\n",@usage)."\n"; 
	print LOG "WARNING: Problematic command line!\n\n";
	die "WARNING: Problematic command line!\n\n";
}


########## Default parameters ################
my $tblastn=" ";
my $blastIdentitycutoff=50;
my $geneCutoff=0.8;
my $genomeCutoff=0.8;
my $num_threads=2;
##########################################
if ($ARGV[4]) {
	print "################# Specified parameters s in parameters.txt #################\n";
	print LOG "################# Specified parameters in parameters.txt #################\n";
	open (PARA, "<$ARGV[4]") or die "Can't open '$ARGV[4]': $!\n";
	foreach my $para (<PARA>) {
		next if $para =~ /\#/;
		next if $para =~ /\*/;
		print "$para";
		print LOG "$para";
		next unless $para =~ /\w/;;
		my @p1=split /	/, $para;
		my @p2=split /\=/,$p1[0];
		if ($p2[0] eq "tblastn") {
			$tblastn=$p2[1];
		} elsif ($p2[0] eq "num_threads") {
			$num_threads =$p2[1];
		} elsif ($p2[0] eq "geneCutoff") {
			$geneCutoff=$p2[1];
		} elsif ($p2[0] eq "blastIdentitycutoff") {
			$blastIdentitycutoff=$p2[1];
		} elsif ($p2[0] eq "genomeCutoff") {
			$genomeCutoff=$p2[1];
		}
	}

} else {
print '################# Parameters (Default) ##############

tblastn=				specify the location of the program tblastn (use double quotation marks if space exists in the path), e.g. D:\"Program Files"\NCBI\blast-2.2.29+\bin\;  

num_threads=2				specify the number of threads used by the program tblastn; 

blastIdentitycutoff=70			specify the cutoff for filterring the tblastn results;

geneCutoff=0.8				specify the cutoff for omitting low-prevalence gene;

genomeCutoff=0.8			specify the cutoff for omitting low-quality genomes;

';

print LOG '################# Parameters (Default) ##############

tblastn=				specify the location of the program tblastn (use double quotation marks if space exists in the path), e.g. D:\"Program Files"\NCBI\blast-2.2.29+\bin\;  

num_threads=2				specify the number of threads used by the program tblastn; 

blastIdentitycutoff=50			specify the cutoff for filter the tblastn results;

geneCutoff=0.8				specify the cutoff for drop low-prevalence gene;

genomeCutoff=0.8			specify the cutoff for drop low-quality genomes;

';
}


my $TEMdir= $inputDir."_TEM";
unless ($ARGV[2] eq "A2") {
     ##################

###############Step 1: Make BlastDB #################
	print "\n############### Step 1: Make BlastDB ###############\n\n";
	print LOG "\n############### Step 1: Make BlastDB ###############\n\n";

	my $cmd0="md $TEMdir";
	system $cmd0 || die "Step 1: Can't create blastDB directory '$TEMdir'\n";
	
	
	
	opendir(DIR, $inputDir) || die "Step 1: Can't open input directory '$inputDir': $!\n";
	my @query = readdir(DIR);
	closedir(DIR);
	die "Step 1: Input directory $inputDir does not contain any files.\n" unless scalar(@query);
	my $cmd1="md $TEMdir\\TEM1_blastDB";
	system $cmd1 || die "Step 1: Can't create blastDB directory '$TEMdir\\TEM1_blastDB': $!\n";
	my $genomeNum = 0;
	
	foreach my $file (@query) {
		next unless $file =~ /\.fas/;
		$genomeNum++;

		my $cmd = "$tblastn"."makeblastdb -in ./$inputDir/$file -dbtype $ARGV[3] -out ./$TEMdir/TEM1_blastDB/$file";
		print LOG "$tblastn"."makeblastdb -in ./$inputDir/$file -dbtype $ARGV[3] -out ./$TEMdir/TEM1_blastDB/$file \n\n";
		
		system "$cmd" || die "Step 1: Can't execute 'makeblastdb': $!\n";
		print "Step 1: Creating BLASTDB: **************It's the $genomeNum task in directory $inputDir!**************\n";
		print LOG "Step 1: Creating BLASTDB: **************It's the $genomeNum task in directory $inputDir!**************\n";
	}
	my $refN=0;
	my @queryR;
	if ($ARGV[2] eq "A3") {
		opendir(DIR, "Reference") || die "Step 1: Can't open input directory 'Reference': $!\n";
		@queryR = readdir(DIR);
		closedir(DIR);
		die "Step 1: Input directory 'Reference' does not contain any files.\n" unless scalar(@queryR);
			foreach my $file (@queryR) {
			next unless $file =~ /\.fas/;
			$genomeNum++;
			$refN++;
			my $cmd = "$tblastn"."makeblastdb -in ./Reference/$file -dbtype $ARGV[3] -out ./$TEMdir/TEM1_blastDB/$file";
		
		
			system "$cmd" || die "Step 1: Can't execute 'makeblastdb': $!\n";
			print "Step 1: Creating BLASTDB: **************It's the $refN task in directory 'Reference'!**************\n";
			print LOG "Step 1: Creating BLASTDB: **************It's the $refN task in directory 'Reference'!**************\n";
		}
	print "Step 1: ##############Totally, $genomeNum BLAST DatabaseCreating were created!##############\n";
	print LOG "Step 1: ##############Totally, $genomeNum BLAST DatabaseCreating were created!##############\n";
	}
	
	
############### Step 2: Blast Search #################
print "\n\n\n############### Step 2: Blast Search ###############\n\n";
print LOG "\n\n\n############### Step 2: Blast Search ###############\n\n";
	opendir(DIR2, $ARGV[1]) || die "Step 1: Can't open input directory '$ARGV[1]': $!\n";
	my @query2 = readdir(DIR2);
	closedir(DIR2);
	die "Step 2: Input Query directory $ARGV[1] does not contain any files.\n" unless scalar(@query2);
	my $num2=0;
	open (QUERY, ">$ARGV[1].seq") or die "Step 1: Can't open input file '$ARGV[1].seq': $!\n";
	foreach my $file2 (@query2) {
		next unless $file2 =~ /\w/;
		$num2++;
		open (FILE2, "<$ARGV[1]/$file2") or die "Step 3: Can't open '$ARGV[1]/$file2': $!\n";
		foreach my $in2 (<FILE2>) {
			next unless $in2;
			if ($in2 =~ /\>/) {
				print QUERY "$in2";
				my @aa2=split /_/, $in2;
		
				push @querylist, $aa2[2];
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
	print "Step 2: The query file $ARGV[1].seq have been created!\n\n";
	print LOG "Step 2: The query file $ARGV[1].seq have been created!\n\n";
	my $cmd2="md $TEMdir\\TEM2_blastOUT";##################################
	system $cmd2 || die "Step 2: Can't create blastOUT directory '$TEMdir/TEM2_blastOUT'\n";
	my $query = "$ARGV[1].seq"; 

	my $number = 0;
	foreach my $file (@query) {
		next unless $file =~ /\.fas/;
		my @in = split (/_/, $file);
		my $cmd = "$tblastn"."$blast -query $query -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';

		system $cmd || die "Step 2: Can't execute '$tblastn': $!\n";
		$number++;
		
		print "Step 2: BLAST: It's the $number/$genomeNum task!\n";
		print LOG "Step 2: BLAST: It's the $number/$genomeNum task!\n";
	}
	my $number0 = 0;
	if ($ARGV[2] eq "A3")	{
		foreach my $file (@queryR) {
			next unless $file =~ /\.fas/;
			my @in = split (/_/, $file);
			my $cmd = "$tblastn"."$blast -query ./Reference/$file -db ./$TEMdir/TEM1_blastDB/$file -out ./$TEMdir/TEM2_blastOUT/Reference_$file -evalue 1e-5 -num_threads $num_threads".' -outfmt "6 qseqid sseqid sstart send pident qcovs evalue bitscore"';
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


############### Step 3: Filtrate Blast Results ################# 符合标准的blast结果都输出，不考虑所有基因组共有情况	
print "\n\n\n############### Step 3: Filtrate Blast Results ###############\n\n";
print LOG "\n\n\n############### Step 3: Filtrate Blast Results ###############\n\n";

	my $cmd3="md $TEMdir\\TEM3_blastOUT_S";
	system $cmd3 || die "Step 3: Can't create blastOUT directory '$TEMdir/TEM3_blastOUT_S'\n";
	opendir(DIR3, "$TEMdir\\TEM2_blastOUT") || die "Step 1: Can't open input directory '$TEMdir\\TEM2_blastOUT': $!\n";
	my @query3 = readdir(DIR3);
	closedir(DIR3);
	die "Step 3: Input  directory $TEMdir\\TEM2_blastOUT does not contain any files.\n" unless scalar(@query3);

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

if ($ARGV[2] eq "A1") {
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
my $cmd33="md $TEMdir\\TEM3_blastOUT_S\\Discard";##########################
system $cmd33 || die "Step 3: Can't create blastOUT directory '$TEMdir\\TEM3_blastOUT_S\\Discard'!\n";
my $dis=0;
foreach my $kk (@outgenome) {
	my $cmd333="move .\\$TEMdir\\TEM3_blastOUT_S\\$kk .\\$TEMdir\\TEM3_blastOUT_S\\Discard\\$kk";
	system $cmd333 || die "Step 3: Can't move '.\\$TEMdir\\TEM3_blastOUT_S\\$kk' into '.\\$TEMdir\\TEM3_blastOUT_S\\Discard\\$kk'!\n";##################21.02.24
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
my $cmd4="md $TEMdir\\TEM4_GeneSeqs";
system $cmd4 || die "Step 4: Can't create blastOUT directory '$TEMdir/TEM4_GeneSeqs'\n";;

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
	if ($ARGV[3] eq "prot") {    ##################21.02.24
	

	
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
			
			
	foreach my $i2 (<FIL>) {
		$i2=~ s/\n// unless $i2=~s/\r\n//;
		next unless $i2;
		my @in0 = split /\s/, $i2;
		my@ff =split /_/, $in0[0];

				
		my $geneSeq = $geno{$in0[1]};	
						
		print OUT ">$ff[3]\n$geneSeq\n";################$ff[3]还是$ff[2]根据query基因标签决定，r89和r95的faa序列标签格式不同。
		$seqnum++;
	}
	close FIL;	
			
	close OUT;
	print "Step 4: I got $seqnum/$num sequences for genome $numg/$ingenomeN:$i1!\n";
	print LOG "Step 4: I got $seqnum/$num sequences for genome $numg/$ingenomeN:$i1!\n";
	
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
		print OUT ">$ff[3]\n$geneSeq\n";################$ff[3]还是$ff[2]根据query基因标签决定，r89和r95的faa序列标签格式不同。
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

my $cmd5="md $TEMdir\\TEM5_GeneCluster";
system $cmd5 || die "Step 5: Can't create blastOUT directory '$TEMdir/TEM5_GeneCluster': $!\n";
my %outClu;
my $num5=0;
print "Step 5: I am handling the";
print LOG "Step 5: I am handling the";
foreach my $in5 (@in1) {
	
	next unless $in5=~ /.\.fas/;
	next unless $in5 ~~ @finalgenome;	##################21.02.24
	my $tt =$in5;
	$tt=~ s/\.fas//;
	open (FILE, "<$TEMdir/TEM4_GeneSeqs/$in5") or die "Step 5: Can't open '$TEMdir/TEM4_GeneSeqs/$in5': $!\n";
	my ($id,$sig);
	foreach my $i1 (<FILE>) {
	
		$i1=~ s/\n// unless $i1=~s/\r\n//;
		next unless $i1;
		if ( $i1 =~ /\>/ ) {
			
			$id =$i1;
			$id=~ s/\>//;
			$sig=1 if $id ~~ @finalgenelist;
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

my $cmd6="md $TEMdir\\TEM6_DedupGeneCluster";
system $cmd6 || die "Step 6: Can't create Deduplicated Seq directory '$TEMdir\\TEM6_DedupGeneCluster': $!\n";
my $cmd60="md $TEMdir\\TEM60_DedupList";
system $cmd60 || die "Step 6: Can't create Deduplication List directory '$TEMdir\\TEM60_DedupList': $!\n";

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
my $cmd7="md $TEMdir\\TEM7_Alignment";###################################
system $cmd7 || die "Step 7: Can't create Alignment directory '$TEMdir\\TEM7_Alignment'\n";

opendir(DIR1, "$TEMdir/TEM6_DedupGeneCluster") || die "Step 7: Can't open input directory '$TEMdir\\TEM6_DedupGeneCluster': $!.\n";
my @in7 = readdir(DIR1);
closedir(DIR1);
die "Step 7: Input directory $TEMdir\\TEM6_DedupGeneCluster does not contain any files.\n" unless scalar(@in7);


my $num7 = 0;

foreach (@in7) {
	next unless /\w/;
	next unless /\.fas/;   #######21.02.24
	my $program = ".\\bin\\muscle.exe";##############################
	my $outname = $_ . ".fasta";
	my @para = ("-in", "./$TEMdir/TEM6_DedupGeneCluster/$_", "-out", "./$TEMdir/TEM7_Alignment/$outname");
	unless (system "$program", @para) {#############################
		$num7++;
		print "\nStep 7: $_ has completed! It is the $num7/$ingene file.\n";
		print LOG "\nStep 7: $_ has completed! It is the $num7/$ingene file.\n";
	} else {
		die "Step 7: Can't execute '.\\bin\\muscle.exe': $!.\n";
	}
}


############### Step 8: Trim Sequence Algnments #################
print "\n\n\n############### Step 8: Trim Sequence Algnments ###############\n\n";
print LOG "\n\n\n############### Step 8: Trim Sequence Algnments ###############\n\n";

my $cmd8="md $TEMdir\\TEM8_AlnTrimmed";
system $cmd8 || die "Step 8: Can't create Alignment directory '$TEMdir\\TEM8_AlnTrimmed'\n";

opendir(DIR, "$TEMdir/TEM7_Alignment") || die "Step 8: Can't open input directory '$TEMdir\\TEM7_Alignment':$!\n";
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
print ".\nStep 8: These $pro clusters are problemic (less than half length was retained for following analysis), please check related file(s) in 'TEM7_Alignment':\n*****Indicated by the filenames and the coverages*****\n@problem\n";
print LOG ".\nStep 8: These $pro clusters are problemic (less than half length was retained for following analysis), please check related file(s) in 'TEM7_Alignment':\n*****Indicated by the filenames and the coverages*****\n@problem\n";

if (@error) {
print "These $err +1 clusters are ERROR:\n@error\n";
print LOG "These $err +1 clusters are ERROR:\n@error\n";
}


############### Step 9: ReCall Sequences #################
print "\n\n\n############### Step 9: ReCall Sequences ###############\n\n";
print LOG "\n\n\n############### Step 9: ReCall Sequences ###############\n\n";

my $cmd9="md $TEMdir\\TEM9_ReCallAln";
system $cmd9 || die "Step 9: Can't create Alignment directory '$TEMdir\\TEM9_ReCallAln': $!\n";

opendir(DIR1, "$TEMdir/TEM8_AlnTrimmed") || die "Step 9: Can't open input directory '$TEMdir\\TTEM8_AlnTrimmed': $!\n";
my @in9 = readdir(DIR1);
closedir(DIR1);
die "Step 9: Input directory $TEMdir\\TEM7_Alignment does not contain any files.\n" unless scalar(@in9);
print "Step 9: I am handling the files(file no./sequence no.):";
print LOG "Step 9: I am handling the files(file no./sequence no.):";
my $num9=0;
foreach my $in (@in9) {
	my $genome=$in;
	next unless $genome=~ s/\.fasta//;
	open (FILE, "<$TEMdir/TEM8_AlnTrimmed/$in") or die "Step 9: Can't open '$TEMdir\\TEM8_AlnTrimmed/$in': $!\n";
	
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
	open (FIL, "<$TEMdir/TEM60_DedupList/$genome.txt") or die "Step 9: Can't open '$TEMdir\\TEM60_DedupList\\$genome.txt': $!\n";
	open (OU, ">./$TEMdir/TEM9_ReCallAln/$genome.fasta") or die "Step 9: Can't open '$TEMdir\\TEM9_ReCallAln\\$genome.fasta': $!\n";
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

opendir(DIR, "$TEMdir/TEM9_ReCallAln") || die "Step 10: Can't open input directory '$TEMdir\\TEM9_ReCallAln': $!\n";
my @in10 = readdir(DIR);
closedir(DIR);
die "Step 10: Input directory $TEMdir\\TEM9_ReCallAln does not contain any files.\n" unless scalar(@in10);

print "Step 10: I am handling the files:";
print LOG "Step 10: I am handling the files:";
my $num10 = 0;
my (%out,@id,@seq,$id,$seq);
foreach my $file ( sort @in10) {
	next unless $file =~ /\.fasta$/;
  	$file =~ /\.fasta$/ || die "Step 10: File '$file' does not have a name in xxxx.fasta format.\n";;
	unless (@id) {
 	
		open (FILE, "$TEMdir/TEM9_ReCallAln/$file") or die "Step 10: Can't open '$TEMdir\\TEM9_ReCallAln/$file': $!\n";
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
		###my @tem;
		foreach my $i2 (0..$#id) {
			$out{$id[$i2]}=$seq[$i2];
		}
	} else {
		open (FILE, "$TEMdir/TEM9_ReCallAln/$file") or die "Step 10: Can't open '$TEMdir\\TEM9_ReCallAln/$file': $!\n";
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
open (OUT, ">$ARGV[0].concatenation.fas") or die "Can't open '$ARGV[0].concatenation.fas': $!\n";
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
if ($ARGV[3] eq "nucl") {
	$cmd11=".\\bin\\FastTree.exe -gtr -nt $ARGV[0].concatenation.fas > $ARGV[0].tree";
} elsif ($ARGV[3] eq "prot") {
	$cmd11=".\\bin\\FastTree.exe $ARGV[0].concatenation.fas > $ARGV[0].tree";
}
system $cmd11 || "Step 10: Can't execute '.\\bin\\FastTree.exe': $!";
}


print "\n\n\n############### Ending Information ###############\n\n";
print LOG "\n\n\n############### Ending Information ###############\n\n";
my ($sec1,$min1,$hour1,$mday1,$mon1,$year1,$wday1,$yday1,$isdst1) = localtime(time);

my $rmon1=$mon1 +1;
my $ryear1=$year1+1900;


if ($ARGV[2] eq "A2") {
	print "The user ordered the second-part (A2) analysis.\nJob was finished at: $sec1:$min1:$hour1,$mday1-$rmon1-$ryear1.\n";
	print LOG "The user ordered the second-part (A2) analysis.\nJob was finished at: $sec1:$min1:$hour1,$mday1-$rmon1-$ryear1.\n";
} else {

	print LOG "Job was finished at: $hour1:$min1:$sec1,$mday1-$rmon1-$ryear1.\n";

	print "Job was finished at: $hour1:$min1:$sec1,$mday1-$rmon1-$ryear1.\n";
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