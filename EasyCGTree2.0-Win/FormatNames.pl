#!/usr/bin/perl  
use warnings;
use strict;
my @usage=qq(
====== EasyCGTree ======
     Version 2.0

Usage: perl formatGenomes.pl genomeFolder

);

die join("\n",@usage)."\n" unless $ARGV[0];      
my $inputDir = $ARGV[0];
opendir(DIR, $inputDir) || die "Can't open input directory '$inputDir'\n";
my @files = readdir(DIR);
closedir(DIR);
die "Input directory $inputDir does not contain any files" unless scalar(@files);

my $number1 = 0;
foreach my $file (@files) {
next unless $file =~ /\w/;
my $file1=$file;
$file1 =~ s/\s/-/g;
$file1 =~ s/_/-/g;
$file1=~ s/\.\w+\z/\.fas/;
$file1 .=".fas" unless $file1=~ /\./;
$file1=~ s/\(//;
$file1=~ s/\)//;

die "Can't rename './$inputDir/$file': $!" unless rename("./$inputDir/$file","./$inputDir/$file1"); 
$number1++;

}
print "\nTotally, $number1 files have been renamed successfully!\n";