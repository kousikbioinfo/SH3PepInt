#!/bin/bash

TEMPDIR=binding-results-$RANDOM;
mkdir $TEMPDIR;
FASTA=$1;
MODEL=$2;
perl fastapl -p1 $FASTA > $TEMPDIR/input_seq_file.txt; ## Reformat each sequence into one line
perl make-peptides.pl -i $TEMPDIR/input_seq_file.txt -o $TEMPDIR/output_seq_file.txt -l 15 -w 5; ## Make peptide sequences having Tyrosine residue at 3rd position along with the sequence id and the position
perl gspan_format_new_GAP.pl -i $TEMPDIR/output_seq_file.txt -o $TEMPDIR/Test-file.gspan -pep -pc -phy -pseq;

echo -e "Seq-ID:\tSequence:" > $TEMPDIR/header1.txt;
#echo $MODEL;
if [ "$MODEL" != "" ] ; then
	basename ${MODEL%.*}|tr '\n' '\t' >> $TEMPDIR/header2.txt;
	svmsgdnspdk_dir/svmsgdnspdk -a TEST -d $TEMPDIR/Test-file.gspan -m $MODEL -gt DIRECTED > $TEMPDIR/out.txt ;
	awk '{print $2}' output.predictions  > $TEMPDIR/output.txt;
	rm output.predictions ;
	paste $TEMPDIR/header1.txt $TEMPDIR/header2.txt > $TEMPDIR/header.txt;
	paste $TEMPDIR/output_seq_file.txt $TEMPDIR/output.txt > $TEMPDIR/values.txt;
	cat $TEMPDIR/header.txt $TEMPDIR/values.txt > $TEMPDIR.txt;
	#awk '{OFS="\t";if ($4>0) print $1,$2,$3,$4}' $TEMPDIR/results.txt > $TEMPDIR.txt;
	
else

 	for i in models/*; do basename ${i%.*}|tr '\n' '\t' >> $TEMPDIR/header2.txt ; 
 	svmsgdnspdk_dir/svmsgdnspdk -a TEST -d $TEMPDIR/Test-file.gspan -m $i -gt DIRECTED > $TEMPDIR/out.txt;
	awk '{print $2}' output.predictions  > $TEMPDIR/output.txt;
	rm output.predictions ; done
 	paste $TEMPDIR/header1.txt $TEMPDIR/header2.txt > $TEMPDIR/header.txt;
 	paste $TEMPDIR/output_seq_file.txt models/*.out > $TEMPDIR/values.txt;
 	cat $TEMPDIR/header.txt $TEMPDIR/values.txt > $TEMPDIR.txt;
 	#awk '{OFS="\t";if ($4>0) print $0}' $TEMPDIR/results.txt > $TEMPDIR.txt;
	rm models/*.out;
fi

echo 'results in '$TEMPDIR'.txt';

rm -r $TEMPDIR;
