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
usage: popgen_tools.py [-h] --vcf_file VCF_FILE --num_chr NUM_CHR
                       [--window_bed WINDOW_BED] [--sfs_out SFS_OUT]
                       [--pi_out PI_OUT] [--no_sfs] [--no_pi]

This script generates a site frequency spectrum and calculates genetic
diversity, pi in nonoverlapping windows (average number of DNA differences
between all pairs of sequence).

optional arguments:
  -h, --help            show this help message and exit
  --vcf_file VCF_FILE   REQUIRED. Input the path to a VCF file
  --num_chr NUM_CHR     REQUIRED. Input the number of chromosomes in your VCF
                        file. The number of alleles is equal to the number of
                        individuals multiplying by 2.
  --window_bed WINDOW_BED
                        Input the path to the BED file that specifies the
                        coordinates for each nonoverlapping window. This
                        argument is required for the analysis of genetic
                        diversity.
  --sfs_out SFS_OUT     Input the path for the output file for the site
                        frequency spectrum. If this parameter is not
                        specified, an output file called sfs.out will be
                        output in the current directory.
  --pi_out PI_OUT       Input the path for the output file for pi. If this
                        parameter is not specified, an output file called
                        pi.out will be output in the current directory.
  --no_sfs              Turn on this flag if you do not want to generate a
                        site frequency spectrum.
  --no_pi               Turn on this flag if you do not want to calculate
                        genetic diversity.

```

## Examples
1. To generate the folded site frequency spectrum and calculate genetic diversity:

```
python popgen_tools.py --vcf_file ~/chr21_10YRI.vcf --num_chr 20 --window_bed ~/chr21_100kb_nonoverlapping_windows.txt
```

2. To generate the folded site frequency spectrum only:

```
python popgen_tools.py --vcf_file ~/chr21_10YRI.vcf --num_chr 20 --window_bed ~/chr21_100kb_nonoverlapping_windows.txt --no_pi
```

3. To calculate genetic diversity only:

```
python popgen_tools.py --vcf_file ~/chr21_10YRI.vcf --num_chr 20 --window_bed ~/chr21_100kb_nonoverlapping_windows.txt --no_sfs
```

4. To plot the folded site frequency spectrum:

```
Rscript plot_sfs.R sfs.out sfs_plot.png
```

5. To plot genetic diversity:

```
Rscript plot_genetic_diversity.R pi.out pi_plot.png
```
