#!/usr/bin/perl
## CONVERTS PHASE (CHROMOPAINTER) FORMAT TO BEAGLE FORMAT
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

sub help {
print("EXTRACTS SNP RANGE FROM PHASE (CHROMOPAINTER) FORMAT\n");

print("usage:   perl phasesubsample.pl <options> <from> <to> <phasefile> <outputphasefile>\n");

print("where:\n");
print("<from>:		First SNP to retain (1 is the first snp)\n");
print("<to>:		Final SNP to retain (L is the last snp, the number on the 3rd line of the phase file)\n");
print("<phasefile>:		ChromoPainter/PHASE style (http://www.hapmap.org/downloads/phasing/2007-08_rel22/phased/00README.txt) SNP file\n");
print("<outputphasefile>:       Output phase file\n\n");

print("<options>:\n");
print("-v: Verbose mode\n");
die "\n";
}

###############################
## ARGUMENT PROCESSING

my $from=0;
my $to=-1;
my $newnumsnps=-1;
my $phasefile="";
my $outfile="";
my $verbose=0;

GetOptions ('v|verbose' => \$verbose);
if(@ARGV != 4) {help();}

$from=$ARGV[0];
$to=$ARGV[1];
$phasefile=$ARGV[2];
$outfile=$ARGV[3];

if(!looks_like_number($from) || !looks_like_number($to) || $from<0 || $to<0 || $to<$from){
    die("Invalid arguments: from and to must be specified and must be in the range 1 to L\n");
}
$newnumsnps=$to-$from+1;

####################################
## Define global variables
my @snplocs; # location of the SNPS

my $numsnps=0; # number of SNPS defined in the file
my $numinds=0; # number of individuals defined in the file
my $numhaps=0; # number of haplotypes observed
my $ploidy=-1; # number of haps per ind

####################################
## File IO

## Check we can read the input files
open PHASEFILE, $phasefile or die $!;

## Create output files
open OUTFILE, ">", $outfile or die $!;

####################################
## Functions we need
sub trim($){  # remove whitespace from beginning and end of the argument
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

####################################
## Read the phasefile 

## read the PHASEFILE header
my $skip=1;
my @tmarr;
while ($skip) {
	my $tmp=<PHASEFILE>;
	my @tmpvals = split(/ /, $tmp);
	if($tmpvals[0] eq "P"){ # found the line with all the SNP locations
		shift @tmpvals;
		@snplocs= split(/ /, $tmp);
		$tmp=<PHASEFILE>; # read the line of S's
		$numsnps=trim(pop @tmarr);
		$numinds=trim(pop @tmarr);
		$skip=0;
	}else {
		push @tmarr, $tmpvals[0];
	}
}
print "Detected $numinds individuals\n";
print "And $numsnps SNPs\n";

if($from>$numsnps || $to>$numsnps){
    die("from or to is larger than the number of SNPs detected ($numsnps) in the file\n");
}

## Print the new header
print OUTFILE "0\n$numinds\n$newnumsnps\nP ";
for(my $i=$from;$i <=$to;++$i){
    print OUTFILE "$snplocs[$i] ";
}
print OUTFILE "\n";
for(my $i=$from;$i <=$to;++$i){
    print OUTFILE "S";
}
print OUTFILE "\n";

# remaining lines are SNPs
while (my $tmp=<PHASEFILE>) {
    if($verbose){print "Reading haplotype $numhaps\n";}
    ++$numhaps;
    my @tarr=split(//,trim($tmp));
    if(scalar(@tarr)!=$numsnps){
	my $tmp=scalar(@tarr);
	die "Expected $numsnps SNPs on haplotype $numhaps, but received $tmp\n";
    }
    for(my $i=$from;$i <=$to;++$i){
	print OUTFILE $tarr[$i];
    }
    print OUTFILE "\n";
}
close PHASEFILE;
close OUTFILE;
