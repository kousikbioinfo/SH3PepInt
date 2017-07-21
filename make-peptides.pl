#!/usr/local/perl/bin/perl
use warnings;
use strict 'vars';
use Getopt::Long;
use Pod::Usage;
use Cwd;


###############################################################################
# PARSE COMMAND LINE OPTIONS
###############################################################################


my ($i_help, $in_file, $out_file, $len, $w_size);   # in- and output files


GetOptions('help'               =>  \$i_help,
	   'i|input=s'          =>  \$in_file,
           'o|output=s'         =>  \$out_file,
           'l|length=s'		=>  \$len,
           'w|window=s'		=>  \$w_size);
           

=head1 OPTIONS

        -help   This message.
	-i      Input file should be in fasta format with one line sequence (NOT multiple lines).
	-o	<FILE_NAME> e.g. whole-proteome.fasta. Name of the output file. The file would be tab seperated
		1 col: protein uniprot-id with location (P31946-1-15)
		2 col: 15 amino acids length sequence.
		If this is not given, the output is written to STDOUT.
	-l	Length of the sequence
	-w	window size 

=cut


check input
pod2usage(-exitstatus => 1, -verbose => 1) if $i_help;

&usage unless $in_file && $out_file;
die "ERROR: File $in_file does not exist.\n" if (!(-e $in_file));

die "ERROR: -l (Length) option is madatory.\n" if (!($len));

die "ERROR: For -w (Window size) option is madatory.\n" if (!($w_size));

	

open (OUT, ">$out_file");
open(FH, "$in_file");
my @file=<FH>;
my $i =0; my $j=1; my ($line, $uniprot_id, $max, $seq, $length, $start, $end);
foreach(@file) {
	$line= $_;
	chomp($line);
	if($line=~ />(\S+)/) {
		$uniprot_id = $1; 
		#open (OUT1, ">/scratch/0/kousik/SH3_Project/data/whole-proteome/proteins/$j.txt"); $j++;
		next;
	}
	$length= length($line);
	$max= $length-$len;
	for ($i=0; $i<=$max; $i=$i+$w_size) {
		$seq= substr($line, $i, $len);
		$start= $i+1;
		$end= $i+$len;
		print OUT "$uniprot_id"."-"."$start"."-"."$end\t"."$seq\n";
		#print OUT ">$uniprot_id"."-"."$start"."-"."$end\n$seq\n";
	}

}
























