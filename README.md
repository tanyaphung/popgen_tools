# popgen_tools
popgen_tools is a program written in Python used to calculate summary statistics commonly used in population genetics. Currently, popgen_tools supports genetic diversity calculation.  

* Note that these scripts are still under development to add more features. 

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
