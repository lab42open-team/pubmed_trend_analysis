# DIG: trenD | Gap analysis of scientific literature

![logo of DIG](DIG_trends_gaps_analysis.png)

## Trend and Gap analysis of PubMed abstracts

This analysis aims to elucidate the trends of research focus and disciplines intersections. Trends will be found using specific keyword search in the PubMed corpus. When a keyword is found in an article the article ID (PMID) and year of publication are saved. Then we search the co-occurrences of keywords in the same article to infer when these keywords - research topics - started to combine.

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
biogeochemical cycle
biogeochemistry
metabolism
metabolome
data mining
ecosystem
ecosystem function
ecosystem service
flux balance analysis
biodiversity
flux
metabolic process
microbial ecology
environment
omics
ecology
```
But these keywords need regular expressions to retrive also plural forms suffixes and complex words

## Prerequisites

In bash the the search is done in awk patterm ```match``` command and results are written in a file after being checked for duplicates and missing information.

For results file transformation, statistics and plotting we used R version 3.5.2 and the ``` tidyverse 1.3.0 ``` packages ``` readr, dplyr, ggplot2, tidyr```.

## Running

The analysis is divided in two scripts, one shell script for keyword mining and one R script for data transformations and plotting. Plots are saved in ```.png``` format in a subdirectory named ``` plots ```.

The user is necessary to pass 3 arguments in the executable bash script. These arguments are passed through flags

* -k |  expects a the path of the txt file with the keywords
* -d |  expects the absolute path to the PubMed data
* -p |  expects a string with prefix of all the generated files (txt and plots)


```
./scripts/dig_analysis.sh -k keywords.txt -d ~/dig_analysis/dataset -p "ecology_trends"

```

When indexing is finished the script calls the r script to generate the plots.

### Search

The heart of the search is the ```match``` function after all abstracts and titles where transformed to lower letters.

### Plots

Three types of plots are created to provide insight to trends of the specific keywords.

#### Keywords frequency

In how many abstracts in PubMed each keyword is mentioned?

![Keyword mentions in PubMed](demo/species_2020-04-14_03-10_pubmed_keyword_frequency.png)

How many abstracts are published each year that contain each keyword?

![Keywords per year in PubMed](demo/species_2020-04-14_03-10_pubmed_keyword_per_year.png)

#### Heatmap

Which keywords are mentioned in the same abstracts and in how many?

This question is awnsered with heatmaps, which illustrate the co-mention frequencies of all keyword pairs. Frequencies of mentions are symmetric so only the upper or lower triangle contain all information.

In order to find the co-mention frequencies we used a simple matrix multiplication. We start from an edgelist containing 2 columns, one with article ID's and one with keywords. We transform this edgelist to a n*m matrix and fill it with 0's and 1's whether a keyword is absent from an article or present, respectively. Then we multiplied this matrix with it's transposed matrix resulting in keyword co-mention frequencies. This method is also called one mode projection of a bipartite graph in graph theory.


![Heatmap](demo/species_2020-04-14_03-17_pubmed_keyword_heatmap.png)

Because the difference of keywords distance can be several orders of multitude we plotted the log2 values to reduce the divergence.

![Heatmap in log scale](demo/species_2020-04-14_03-17_log_pubmed_keyword_heatmap.png)


#### Jaccard similarity

We also wanted to quantify the co-mention of keywords in terms of similarity. We choose the Jaccard similarity because is intuitive and easy to calculate; it is defined as the intersection over the union of two sets. Hence it is located in the [0,1] space.

![Heatmap with Jaccard similarity](demo/species_2020-04-14_03-17_pubmed_jaccard_heatmap.png)

#### Random expectation 

In this heatmap we compare the co-mention frequency with random expectation. 

All frequencies are measured with the total unique abstracts in PubMed.

## Two Keyword co-occurences

Sometimes it is insightful to view keyword co-occurences per year to view the trend. To do this we created a different script. The user must specify the data generated from the basic script and then two keywords. The script automatically checks if the both keywords exist in the data. If they exist the script calculates the abstracts per year that both keywords are mentioned. Then it generates a plot with each keyword and then their co-occurrance.

```
Rscript trends_2_keywords_plot.r data_file "keyword1" "keyword2" "file_prefix"

```

