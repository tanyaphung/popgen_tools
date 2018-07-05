# popgen_tools
This repository contains scripts that are often used in population genetics analyses (generating the site frequency spectrum or calculating genetic diversity)

## Getting started

* To get started, clone the directory (via a command-line):

```
git clone https://github.com/tnphung/popgen_tools.git
```

## Prepare the input files

## Usage

* This script generates both the folded site frequency spectrum and calculates genetic diversity in nonoverlapping windows from the VCF file

* For usage, do:

```
python popgen_tools.py -h                                       
usage: popgen_tools.py [-h] --vcf_file VCF_FILE --target_bed TARGET_BED
                       --names_list NAMES_LIST [--sfs_out SFS_OUT]
                       [--pi_out PI_OUT] [--no_sfs] [--no_pi]

This script generates a site frequency spectrum and calculates genetic
diversity, pi. Pi is defined as the average number of DNA differences between
all pairs of sequence).

optional arguments:
  -h, --help            show this help message and exit
  --vcf_file VCF_FILE   REQUIRED. Input the path to a VCF file. Either gzipped
                        or not gzipped file works.
  --target_bed TARGET_BED
                        REQUIRED. Input the path to the BED file that
                        specifies the coordinates to generate SFS or to
                        calculate pi. For example, this file specificies the
                        coordinates for the putatively neutral regions
  --names_list NAMES_LIST
                        REQUIRED. Input the path to the file that lists the
                        individuals from the VCF file that you want to
                        calculate genetic diversity or to generate the SFS.
  --sfs_out SFS_OUT     Input the path for the output file for the site
                        frequency spectrum. If this parameter is not
                        specified, an output file called sfs.out will be
                        outputted in the current directory.
  --pi_out PI_OUT       Input the path for the output file for pi. If this
                        parameter is not specified, an output file called
                        pi.out will be outputted in the current directory.
  --no_sfs              Turn on this flag if you do not want to generate a
                        site frequency spectrum.
  --no_pi               Turn on this flag if you do not want to calculate
                        genetic diversity.

```
