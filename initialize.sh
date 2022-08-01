#!/usr/bin/sh

# first clone repo and cd into it: 
# git clone https://github.com/alperyilmaz/rnaseq-staj.git
# cd rnaseq-staj
# then run this script with ./initialize.sh
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
chmod 755 Mambaforge-Linux-x86_64.sh
./Mambaforge-Linux-x86_64.sh -b
source .bashrc
mamba install -c bioconda -c conda-forge  snakemake-minimal
mamba install pandas
rm Mambaforge-Linux-x86_64.sh
#snakemake -c16 --use-conda
