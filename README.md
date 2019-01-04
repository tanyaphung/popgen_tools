# popgen_tools
popgen_tools is a program written in Python used to calculate summary statistics commonly used in population genetics. The main input of popgen_tools is a VCF (variant calling format) file. Currently, popgen_tools supports genetic diversity calculation.  

## Getting started
* To get started, clone the directory (via a command-line):
```
git clone https://github.com/tnphung/popgen_tools.git
```

## Dependencies
* popgen_tools is tested with Python3.7 and requires the following packages:
  - pandas (version 0.23.4)

## Calculating genetic diversity
* Genetic diversity is defined as pi, average number of differences between pairs of sequences in the sample (Tajima 1983). 

### 1. Calculate genetic diversity using all of the variants in the VCF file
* Scenario: When you have a VCF file and you want to know what pi for all of the individuals or a subset of individuals in the VCF file.
* Usage:
```
python
```
