# Trend analysis of PubMed abstracts for microbial ecology and ecology informatics

This analysis aims to elucidate the trends of research focus and disciplines intesections. Trends will be found using specific keyword search in the PubMed corpus. When a keyword is found in an article the article ID (PMID) and year of publication are saved. Then we search the co-occurrences of keywords in the same article to infer when these keywords - research topics - started to combine.

## PubMed corpus

The data used for this analysis are the abstracts and article information of the PubMed database. PubMed corpus is freely available to download through * ftp * and the E-utilities API. For more information see [PubMed Download](https://www.nlm.nih.gov/databases/download/pubmed_medline_documentation.html).

The output of the ftp method are multiple xml files which then were transformed to tab separated files. All these files are zipped with the ``` gzip command```. Each .tsv.gz file contains 6 fields with the following structure:

 | PMID DOI |  Authors | Journal.volume:pages  | year  | title | Abstract  |
 | --- | --- | --- | --- | --- | --- |

## Keywords

A .txt file with the following keywords is read from the bash script

```
microbiome
text mining
biogeochemical cycles
biogeochemistry
metabolism
metabolome
data mining
ecosystem
ecosystem function
ecosystem service
flux balance analysis
biodiversity
fba
flux
metabolic process
microbial ecology
environment
omics
ecology
```

## Prerequisites

In bash the the search is done by the ``` zgrep ``` command and results are written in a file using ``` awk ```.

For results file transformation, statistics and plotting we used R version 3.5.2 and the ``` tidyverse 1.3.0 ``` packages ``` readr, dplyr, ggplot2, tidyr```.

## Running

### Directory

### Search

### Plots



