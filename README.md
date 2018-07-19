# popgen_tools
This repository contains scripts that are often used in population genetics analyses (generating the site frequency spectrum or calculating genetic diversity). 

* Note that these scripts are still under development to add more features. 

## Getting started

* To get started, clone the directory (via a command-line):

```
git clone https://github.com/tnphung/popgen_tools.git
```

## Usage and Examples

* This script generates a site frequency spectrum (SFS) and calculates genetic diversity (pi). Users have the option turning on or off parts of the pipeline (i.e. generating the SFS only or calculating pi only). 
* The example command-line is written using the example files provided in the folder `example_input_files`.

1. Calculate genetic diversity only
* Here, genetic diversity is defined as the average number of differences in DNA across all pairs of sequences in a sample

  a. Calculate genetic diversity within regions of the genome that is specified by a BED file
    - We are typically interested in calculating genetic diversity in certain regions of the genome, such as in the neutral regions of the genome. 
    - The command is:
    ```
    python popgen_tools.py --vcf_file <path/to/VCF> --target_bed <path/to/BED> --pi_out <path/to/output/for/pi> --total_SNPs <path/to/output/for/total_SNPs/in/regions> --no_sfs
    ```
    - Basically, one needs to supply a VCF file, a BED file specifing regions of the genome where you want to calculate genetic diversity, and the file names for output.
    - Spefifically, for the files provided in `example_input_files`:
    ```
    python popgen_tools.py --vcf_file example_input_files/example_vcf.vcf.gz --target_bed example_input_files/example_neutral_regions.bed --pi_out example_output_files/pi.out --total_SNPs example_output_files/total_SNPs.out --no_sfs
    ```
    - The above command will calculate genetic diversity in each region of the genome whose start and end coordinates are specified by the BED file for all of the individuals in your VCF file. If you want to calculate genetic diversity for only a subset of individuals in your sample, the command is:
    
    ```
    python popgen_tools.py --vcf_file <path/to/VCF> --target_bed <path/to/BED> --pi_out <path/to/output/for/pi> --total_SNPs <path/to/output/for/total_SNPs/in/regions> --names_list <path/to/a/list/of/individuals> --no_sfs
    ```
    - For example:
    ```
    python popgen_tools.py --vcf_file example_input_files/example_vcf.vcf.gz --target_bed example_input_files/example_neutral_regions.bed --pi_out example_output_files/pi_subset.out --total_SNPs example_output_files/total_SNPs_subset.out --no_sfs --names_list example_input_files/individuals_list.txt
    ```
     - The file `individuals_list.txt` lists the names of the individuals whose genetic diversity you want to calculate. List each individual per line. 
     - The output `pi.out` file has 3 columns. The first column is the start coordinate, the second column is the end coordinate, and the third column is pi. 
     - Typically we are interested in pi/site. To calculate pi/site, we can do this in R. I will incorporate this in the Python script soon. 
     ```
     module load R
     R
     data = read.table("pi.out")
     head(data)
        V1       V2 V3
        1  4748512  4748534  0
        2  4609189  4609191  0
        3  4702756  4702793  0
        4  4764986  4764995  0
        5 15711545 15711549  0
        6  4917404  4917405  0
     pi_per_site = sum(data[,3])/sum(data[,2]-data[,1])
     pi_per_site
     [1] 0.000642915
     ```
     
  b. Calculate genetic diversity in non-overlapping windows
    - Often time, we are interested in calculating genetic diversity in each non-overlapping windows where we divide the genome into windows of 50kb, 100kb, or 1Mb (or any size of choice). 
    - In each window, we would like to calculate genetic diversity. 
    - Right now, the scripts are written such that in each non-overlapping window, we want to calculate genetic diversity within regions that are specified by another BED file. I will add options to remove this soon!
    - For this, the VCF file should be first filtered to contain only the SNPs that are found in the BED file. GATK can be used to do this. TODO: do this in the Python script, instead of using GATK
    - The command below calculate genetic diversity for each window using regions that are specified by the BED file. It returns a file where the first column is the start coordinate of the window, the second column is the end coordinate of the window, the third column is the total number of callable sites within that window, and the final column is pi in that window. The reason why we are interested in the total number of callable sites is because we are ultimately interested in calculating pi/site, which is equal to pi divided by the total number of callable sites. 
    ```
    python popgen_tools.py --vcf_file <path/to/VCF> --target_bed <path/to/BED> --pi_out <path/to/output/for/pi>  --window_bed <path/to/BED/file/specifying/windows> --window --no_sfs
    ```
    - Specifically:
    
