# popgen_tools
popgen_tools is a program written in Python used to calculate summary statistics commonly used in population genetics. The main input of popgen_tools is a VCF (variant calling format) file. Currently, popgen_tools supports genetic diversity calculation and generation of a site frequency spectrum.

Updated: 04/08/2020
## Generate a site frequency spectrum
```
python popgen_tools.py --vcf_file {input} --sfs_all --sfs_all_out {output} --ploidy {insert haploid or diploid}
```
- The above command generates the SFS using all of the sites within the VCF file. This command will not do any filtering of variants. Therefore, if you want to generate a site frequeny spectrum for a subset of your sites in the VCF file, you would need to do the subset prior to running this command. 

## Calculate genetic diversity
```
python popgen_tools.py --vcf_file {input} --pi --ploidy {insert haploid or diploid} --pi_all
```
