SH3PepInt
===

SH3PepInt is a prediction tool for SH3-peptide interactions, based on Graph Kernel using SVM. 
Current version contains a total 70 human SH3 domain models. 

SH3PepInt 1.0, January, 2013 

Authors: Kousik Kundu, Fabrizio Costa and Rolf Backofen

Platform:
------------

Unix and Linux


Installation:
------------

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


Webserver:
-----------------
http://modpepint.informatik.uni-freiburg.de/SH3PepInt/Input.jsp


Contact:
-----------------
Kousik Kundu (kk8@sanger.ac.uk)


Citation:
-----------------
* Kousik Kundu, Martin Mann, Fabrizio Costa, and Rolf Backofen.
[MoDPepInt: An interactive webserver for prediction of modular domain-peptide interactions
Bioinformatics, 2014.](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu350)

* Kousik Kundu, Fabrizio Costa, and Rolf Backofen.
[A graph kernel approach for alignment-free domain-peptide interaction prediction with an application to human SH3 domains
Bioinformatics, 29(13), pp. i335-i343, 2013.](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btt220)




