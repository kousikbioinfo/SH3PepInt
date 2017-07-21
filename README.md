SH3PepInt
===


SH3PepInt is a prediction tool for SH3-peptide interactions, based on Graph Kernel using SVM. 
Current version contains a total 70 human SH3 domain models. 

SH3PepInt 1.0
January, 2013 
Authors: Kousik Kundu, Fabrizio Costa and Rolf Backofen

Platform:
------------

Unix and Linux


Installation:
------------

A current version of SH2PepInt you can get from:

http://www.bioinf.uni-freiburg.de/Software/SH3PepInt
 
or please write an email to Kousik Kundu <kk8@sanger.ac.uk>

To install the tool, please extract the src archive somewhere. Then change
into that directory and type

  bash COMPILE.sh

the script compiles the svmsgdnspdk and  
creates the master script, namely SH3PepInt.sh. 




Dependency:
-------------

In order to compile SH3PepInt correctly, you need "PERL" already installed.



Usage:
--------------

SH3PepInt <protein/peptide fasta file> <model file>

e.g. SH3PepInt sample.fasta models/SRC-84-145.model

The model file is OPTIONAL, by default the program uses all the 70 models for the predictions. 




Files:
-----------------
sample.fasta: sample file for run the tool.
peptides.fasta: sample peptide file with 15 amino acids length.






