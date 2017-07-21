#!/usr/local/perl/bin/perl
use warnings;
use strict 'vars';
use Getopt::Long;
use Pod::Usage;
use Cwd;


###############################################################################
# PARSE COMMAND LINE OPTIONS
###############################################################################

=head1 NAME
gspan_format.pl

=head2 DESCRIPTION
This script will make a gspan file for graph kernel input. The output
file order is NOT correspodance to the FLAG options. The order would 
be always as following: 
"peptide-sequence, peptide-polarity, peptide-charge, peptide-hydrophobicity
domain-sequence, domain-polarity, domain-charge, domain-hydrophobicity".
N.B. There is a file called by this script which contains all the 
domains name and the sequence in two column.

e.g: perl gspan_format.pl -i input.tsv -o output.gspan -pep -pp -dom -dc

AUTHOR: Kousik Kundu

=head1 OPTIONS

        -help   This message.
	-i      Input file should be in tab seperated format with 1st and 2nd column 
		with domain name and interacting peptide, respectively.
	-o	<FILE_NAME> e.g. GRB2.gspan. Name of the output file. 
**
FLAGS:	-pep	<default no>	Tag for the peptide analysis. if no agrument is given, 
				it will identify only the amino acids in a peptide sequence.
	-pp 	<default no>	Indicates the polarity of each amino acids in a peptide sequence.
	-pc 	<default no>	Indicates the charge (postitive, negative, neutral) of each amino 
				acids in a peptide sequence.
	-pseq	<default no>	Identify only the amino acids in a peptide sequence.
	-phy 	<default no>	Indicates the hydrophobicity of each amino acids in a peptide sequence.
	-dom 	<default no>	Tag for the peptide analysis. if no agrument is given, 
				it will identify only the amino acids in a domain sequence.
	-dp 	<default no>	Indicates the polarity of each amino acids in a domain sequence.
	-dc 	<default no>	Indicates the charge (postitive, negative, neutral) of each amino 
				acids in a domain sequence.
	-dhy 	<default no>	Indicates the hydrophobicity of each amino acids in a domain sequence.
	-dseq	<default no>	Identify only the amino acids in a domain sequence.
	-bothway <defauly no>	Indicates the edges would be both way.. otherwise single way



=cut

#### COMMAND LINE OPTIONS 

my ($i_help, $in_file, $out_file, $pep, $dom, $pseq, $dseq, $pp, $dp, $pc, $dc, $phy, $dhy, $bway);  


GetOptions('help'               =>  \$i_help,
	   'i|input=s'          =>  \$in_file,
           'o|output=s'         =>  \$out_file,
	   'pep|peptide'	=>  \$pep,
	   'dom|domain'		=>  \$dom,
	   'pseq|peptide'	=>  \$pseq,
	   'dseq|domain'	=>  \$dseq,
	   'pp|pep_pol'      	=>  \$pp,
	   'dp|dom_pol'		=>  \$dp,
	   'pc|pep_charge'	=>  \$pc,
	   'dc|dom_charge'	=>  \$dc,
	   'phy|pep_hydro'	=>  \$phy,
	   'dhy|dom_hydro'	=>  \$dhy,
	   'bway|bothway'	=>  \$bway);
	
           



check input
pod2usage(-exitstatus => 1, -verbose => 1) if $i_help;

&usage unless $in_file && $out_file;
die "ERROR: File $in_file does not exist.\n" if (!(-e $in_file));
#die "ERROR: File $out_file does already exist.\n" if (-e $out_file);

open (OUT, ">$out_file");
open(FH, "$in_file");
my @file=<FH>; 

my ($i, $dom_name, $pep_seq, $dom_seq, $count);
my $x =1;my $gap = 0;
foreach(@file) {
	$count = 1; 
	my $line= $_;
	chomp($line);
	my @info= split(/\t/, $line);
	$dom_name = $info[0];
	$pep_seq= $info[1]; 
	$pep_seq=uc($pep_seq);
	print OUT "t\t\# id/instance $x\n"; $x++;
	open(OP, "/home/kousik/Projects/SH3_Domains/data/SH3_domain_MUSCLE.out");
	while(<OP>) {
	
	my $f_line=$_;
	chomp($f_line);
	my @f_info= split(/\t/, $f_line);
	if($f_info[0]=~ /$dom_name/) {
	$dom_seq = $f_info[1]; 
	$dom_seq=uc($dom_seq); last;
	}
	
	}
	$gap = 0;
	if($pep) {
	pep_seq_amino_acids_any($pep_seq);
	}
	
	if($pp) {
	die "ERROR: For -pp, the -pep option is madatory.\n" if (!($pep));
	pep_seq_polarity($pep_seq);
	}
	if($pc) {
	die "ERROR: For -pc, the -pep option is madatory.\n" if (!($pep));
	pep_seq_charge($pep_seq);
	}
	if($phy) {
	die "ERROR: For -phy, the -pep option is madatory.\n" if (!($pep));
	pep_seq_hydrophobicity($pep_seq);
	}
	if($pseq) {
	die "ERROR: For -pseq, the -pep option is madatory.\n" if (!($pep));
	pep_seq_amino_acids($pep_seq);
	}
	
	if($dom) {
	dom_seq_amino_acids_any($dom_seq);
	}
	
	if($dp) {
	die "ERROR: For -dp, the -dom option is madatory.\n" if (!($dom));
	dom_seq_polarity($dom_seq);
	}
	if($dc) {
	die "ERROR: For -dc, the -dom option is madatory.\n" if (!($dom));
	dom_seq_charge($dom_seq);
	}
	if($dhy) {
	die "ERROR: For -dhy, the -dom option is madatory.\n" if (!($dom));
	dom_seq_hydrophobicity($dom_seq);
	}
	if($dseq) {
	die "ERROR: For -dseq, the -dom option is madatory.\n" if (!($dom));
	dom_seq_amino_acids($dom_seq);
	}


}


sub pep_seq_amino_acids_any {
my $p_seq = shift;
chomp($p_seq); 
my $len=length($p_seq); 
my $j=$count;
for(my $i=0; $i<$len; $i++) {
	#my $aa = substr($p_seq,$i,1);
	print OUT "v\t$count\tX\n"; $count++;
	}

while($j<$count-1) {
my $m=$j+1;
print OUT "e\t$j\t$m\n"; $j++;
}
}

#############################


#############################

sub pep_seq_amino_acids {
my $p_seq = shift;
chomp($p_seq); 
my $len=length($p_seq); 
my $j=$count;
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($p_seq,$i,1);
	print OUT "v\t$count\tPA$aa\n"; $count++;
	}

while($j<$count) {
my $m=$j-$len;


if($bway) {
print OUT "e\t$j\t$m\n";
}
print OUT "e\t$m\t$j\n"; $j++;
}
}

#############################



sub pep_seq_polarity {
my $p_seq = shift;
chomp($p_seq); 
my $len=length($p_seq); 
my $j=$count;
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($p_seq,$i,1);
	if($aa eq "A" || $aa eq "G" || $aa eq "I" || $aa eq "L" || $aa eq "M" || $aa eq "F" || $aa eq "P" || $aa eq "W" || $aa eq "V") {
	my $pol = "PPN";
	print OUT "v\t$count\t$pol\n"; $count++;
	}
	else {
	my $pol = "PPP";
	print OUT "v\t$count\t$pol\n"; $count++;
	}
}

while($j<$count) {
my $m=$j-$len;

if($bway) {
print OUT "e\t$j\t$m\n";
}

print OUT "e\t$m\t$j\n"; $j++;
}
}

#############################

sub pep_seq_charge {
my $p_seq = shift;
chomp($p_seq); 
my $len=length($p_seq); 
my $j=$count;
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($p_seq,$i,1);
	if($aa eq "D" || $aa eq "E") {
	my $ch= "PCA";
	print OUT "v\t$count\t$ch\n"; $count++;
	}
	elsif($aa eq "R" || $aa eq "K" || $aa eq "H" ) {
	my $ch= "PCB";
	print OUT "v\t$count\t$ch\n"; $count++;
	}
	else {
	my $ch= "PCN";
	print OUT "v\t$count\t$ch\n"; $count++;
	}
}
while($j<$count) {
my $m=$j-$len;

if($bway) {
print OUT "e\t$j\t$m\n";
}

print OUT "e\t$m\t$j\n"; $j++;
}
}

#############################

sub pep_seq_hydrophobicity {
my $p_seq = shift;
chomp($p_seq); 
my $len=length($p_seq); 
my $j=$count;
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($p_seq,$i,1);
	if($aa eq "L" || $aa eq "V" || $aa eq "I") {
	my $hy = "PHVH";
	print OUT "v\t$count\t$hy\n"; $count++;
	}
	elsif($aa eq "A" || $aa eq "M" || $aa eq "C" || $aa eq "F") {
	my $hy = "PHH";
	print OUT "v\t$count\t$hy\n"; $count++;
	}
	elsif($aa eq "G" || $aa eq "T" || $aa eq "S" || $aa eq "W" || $aa eq "Y" || $aa eq "P") {
	my $hy = "PHL";
	print OUT "v\t$count\t$hy\n"; $count++;
	}
	else {
	my $hy = "PHVL";
	print OUT "v\t$count\t$hy\n"; $count++;
	}
}

while($j<$count) {
my $m=$j-$len;

if($bway) {
print OUT "e\t$j\t$m\n";
}

print OUT "e\t$m\t$j\n"; $j++;
}
}


#############################

	
sub dom_seq_amino_acids_any {
my $d_seq = shift;
chomp($d_seq); 
my $len=length($d_seq); 
my $j=$count;
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($d_seq,$i,1);
	my $pos=$i+1;
	print OUT "v\t$count\tD$pos\n"; $count++;
	}

while($j<$count-1) {
my $m=$j+1;
print OUT "e\t$j\t$m\n"; $j++;
}
}

#############################

	
sub dom_seq_amino_acids {
my $d_seq = shift;
chomp($d_seq); 
my $len=length($d_seq); 
my $j=$count; my @gap_c=""; my $a=0;
my $gap_t=0;
if($gap == 1) {
	$gap_t=1;
}
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($d_seq,$i,1);
	if($aa eq "-") {
		if($gap_t == 1) {
			$gap_c[$a]= "$count"; $a++; $count++; next;
		}
		$gap = 1;
		print OUT "v\t$count\tGAP\n"; $count++;
	}
	else {
		print OUT "v\t$count\tDA$aa\n"; $count++;
	}
}

label: while($j<$count) {
my $m=$j-$len;

foreach(@gap_c) {
	my $line =$_; 
	chomp($line); 
	if($j eq $line) { 
		$j++; next label;
	}
	}
if($bway) {
print OUT "e\t$j\t$m\n";
}

print OUT "e\t$m\t$j\n"; $j++;
}
}

#############################


sub dom_seq_polarity {
my $d_seq = shift;
chomp($d_seq); 
my $len=length($d_seq); 
my $j=$count; my @gap_c=""; my $a=0;
my $gap_t=0;
if($gap == 1) {
	$gap_t=1;
}
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($d_seq,$i,1);
	if($aa eq "A" || $aa eq "G" || $aa eq "I" || $aa eq "L" || $aa eq "M" || $aa eq "F" || $aa eq "P" || $aa eq "W" || $aa eq "V") {
	my $pol = "DPN";
	print OUT "v\t$count\t$pol\n"; $count++;
	}
	elsif ($aa eq "-") {
	my $pol = "GAP";
	if($gap_t == 1) {
		$gap_c[$a]= "$count"; $a++; $count++; next;
	}
	$gap = 1;
	print OUT "v\t$count\t$pol\n"; $count++;
	}
	else {
	my $pol = "DPP";
	print OUT "v\t$count\t$pol\n"; $count++;
	}
}

label: while($j<$count) {
my $m=$j-$len;
print @gap_c;
foreach(@gap_c) {
	my $line =$_; 
	chomp($line); 
	if($j eq $line) { 
		$j++; next label;
	}
	}
if($bway) {
print OUT "e\t$j\t$m\n";
}

print OUT "e\t$m\t$j\n"; $j++;
}
}

#############################


sub dom_seq_charge {
my $d_seq = shift;
chomp($d_seq); 
my $len=length($d_seq); 
my $j=$count;  my @gap_c=""; my $a=0;
my $gap_t=0;
if($gap == 1) {
	$gap_t=1;
}
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($d_seq,$i,1);
	if($aa eq "D" || $aa eq "E") {
	my $ch= "DCA";
	print OUT "v\t$count\t$ch\n"; $count++;
	}
	elsif($aa eq "R" || $aa eq "K" || $aa eq "H" ) {
	my $ch= "DCB";
	print OUT "v\t$count\t$ch\n"; $count++;
	}
	elsif($aa eq "-") { 
	my $ch= "GAP";
	if($gap_t == 1) { 
		$gap_c[$a]= "$count";$a++; $count++; next;
	}
	$gap = 1;
	print OUT "v\t$count\t$ch\n"; $count++;
	}
	else {
	my $ch= "DCN";
	print OUT "v\t$count\t$ch\n"; $count++;
	}
}

label: while($j<$count) {
my $m=$j-$len;
foreach(@gap_c) {
	my $line =$_; 
	chomp($line); 
	if($j eq $line) { 
		$j++; next label;
	}
	}
if($bway) {
print OUT "e\t$j\t$m\n";
}

print OUT "e\t$m\t$j\n"; $j++;
}
}

#############################


sub dom_seq_hydrophobicity {
my $d_seq = shift;
chomp($d_seq); 
my $len=length($d_seq); 
my $j=$count;  my @gap_c=""; my $a=0;
my $gap_t=0;
if($gap == 1) {
	$gap_t=1;
}
for(my $i=0; $i<$len; $i++) {
	my $aa = substr($d_seq,$i,1);
	if($aa eq "L" || $aa eq "V" || $aa eq "I") {
	my $hy = "DHVH";
	print OUT "v\t$count\t$hy\n"; $count++;
	}
	elsif($aa eq "A" || $aa eq "M" || $aa eq "C" || $aa eq "F") {
	my $hy = "DHH";
	print OUT "v\t$count\t$hy\n"; $count++;
	}
	elsif($aa eq "G" || $aa eq "T" || $aa eq "S" || $aa eq "W" || $aa eq "Y" || $aa eq "P") {
	my $hy = "DHL";
	print OUT "v\t$count\t$hy\n"; $count++;
	}
	elsif($aa eq "-") {
	my $hy= "GAP";
	if($gap_t == 1) {
		$gap_c[$a]= "$count"; $a++; $count++; next;
	}
	$gap = 1;
	print OUT "v\t$count\t$hy\n"; $count++;
	}
	else {
	my $hy = "DHVL";
	print OUT "v\t$count\t$hy\n"; $count++;
	}
}

label: while($j<$count) {
my $m=$j-$len;
foreach(@gap_c) {
	my $line =$_; 
	chomp($line); 
	if($j eq $line) { 
		$j++; next label;
	}
	}
if($bway) {
print OUT "e\t$j\t$m\n";
}

print OUT "e\t$m\t$j\n"; $j++;
}

}












































