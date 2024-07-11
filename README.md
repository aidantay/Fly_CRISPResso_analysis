# Fly CRISPResso Analysis

This repository contains various scripts for analysing the allele frequency tables from [CRISPResso](https://github.com/pinellolab/CRISPResso2)

## Installation ## 

Create the conda environment: `conda create env -f environment.yml`.

## Usage ##

**Convert allele frequency table to SAM format**
```
python analyse_crispresso_allele_freq.py freqtosam \
    Alleles_frequency_table.txt
```

**Extract reference sequence from allele frequency table**
```
python analyse_crispresso_allele_freq.py extractref \
    Alleles_frequency_table.txt
```

**Align gRNA to reference sequence** 
```
python analyse_crispresso_allele_freq.py grnatosam \
    Alleles_frequency_table.txt \
    GCCACAATTGTCGATCGTCA CACTGCTGGCCATCTGGAAG
```

**Find and count stop codons**
```
python analyse_crispresso_allele_freq.py countstopcodons \
    Alleles_frequency_table.txt \
    GCCACAATTGTCGATCGTCA CACTGCTGGCCATCTGGAAG
```

**Count mutations**
```
python analyse_crispresso_allele_freq.py countmutation \
    Alleles_frequency_table.txt \
    GCCACAATTGTCGATCGTCA CACTGCTGGCCATCTGGAAG
```

**Count frameshift mutations**
```
python analyse_crispresso_allele_freq.py countframeshift \
    Alleles_frequency_table.txt \
    GAAGTGGCATTGAAGAGCCAACCAGGATGTCAGTCCTTCTCACCACCTCAGTATTATTCGCAGCTGATAACAGAG \
    EVALKSQPGCQSFSPPQYYSQLITEMEILGWD
```

## Contribution guidelines ##

The source codes are licensed under GPL less public licence. Users can contribute by making comments on the issues tracker, the wiki or direct contact via e-mail (see below).

## Contact ##

* Aidan Tay: a.tay@unswalumni.com
* Michael Clark: michael.clark@mq.edu.au
