#!/bin/bash

# command name: mkblastdb.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 03 AGO 2023
# objetive: Create blast dbs
# Syntax: ./mkblastdb.sh <parameters>
# Link: https://www.ncbi.nlm.nih.gov/books/NBK569841/

# Syntax: ./mkblastdb.sh <path/dataset> <dbtype: nucl or prot> <dbname>

# Define variables
DATASET=$1
DBTYPE=$2
DBNAME=$3
DBDIR=data/BLASTDB

[ ! -d $DBDIR/$DBNAME ] || mkdir $DBDIR/$DBNAME

# Makeblastdb
makeblastdb -in $DATASET -parse_seqids -dbtype $DBTYPE -out $DBDIR/$DBNAME

