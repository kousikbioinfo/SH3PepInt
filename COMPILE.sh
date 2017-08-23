#!/bin/bash

#########################################
# function for error handling
#########################################

function error_exit
{
    echo ${1:-"Unknown Error"} 1>&2
    exit 1
}

perl --help > /dev/null || error_exit "No Perl installation available. Please install first!";

mkdir -p $INSTALL_PREFIX/bin;

#########################################
echo " compile svmsgdnspdk tools";
#########################################

cd svmsgdnspdk_dir;
make;
cd ..;

#########################################
echo " compilation and installation finished!";
#########################################

echo
echo " You can run now ./SH3PepInt.sh sample.fasta models/SRC-84-145.model";
echo

