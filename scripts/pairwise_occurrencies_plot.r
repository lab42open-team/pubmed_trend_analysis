#!/usr/bin/Rscript

# Author:   Savvas Paragkamian (s.paragkamian@hcmr.gr)
#           Institute of Marine Biology Biotechnology and Aquaculture (IMBBC)
#           Hellenic Centre for Marine Research (HCMR)
#
# Created:  05/06/2020
# License:  GNU GPLv3 license
# Developer : Savvas Paragkamian
# This script is the analysis of co-occurrence of keywords. This requires manual input.


suppressPackageStartupMessages({
    library(tidyverse) # version tidyverse 1.3.0, used packages from tidyverse are readr_1.3.1, ggplot2_3.3.0, dplyr_0.8.5
    library(Matrix) # version Matrix 1.2-15
    library(scales)
})


# load the user arguments
#args <- commandArgs(trailingOnly=TRUE)

args <- c("../data/beach_litter_2021-04-13_01-30_dig_analysis.tsv", "marine litter", "macroplastic","test")
keyword1 <- args[2]
keyword2 <- args[3]
keywords_1_and_2 <- paste0(keyword1,"_AND_",keyword2)
user_prefix <- args[4]

## file loading
trends_pubmed <- read_delim(args[1], delim="\t", col_names=F,col_types = cols())

colnames(trends_pubmed) <- c("PMID","year","keyword")
#colnames(trends_pubmed) <- c("line","PMID","keyword","file","year")

all_keywords <- unique(trends_pubmed$keyword)

# check if the provided exist in the data file

check_keywords <- c(keyword1,keyword2) %in% all_keywords

if (check_keywords[1]==F){
  print(paste0(keyword1, " is not in the dataset"))

} 

if (check_keywords[2]==F){
  print(paste0(keyword2, " is not in the dataset"))
}

if (F %in% check_keywords){

  stop("missing keywords from data")
}

keywords_id_per_year <- trends_pubmed %>% filter(keyword %in% c(keyword1,keyword2)) %>% distinct(PMID, keyword,year) %>% mutate(values=1) %>% pivot_wider(names_from=keyword,values_from=values)

keywords_AND <- keywords_id_per_year %>% mutate(oneandtwo=if_else(.[[3]]==1 & .[[4]]==1,1,0)) %>% pivot_longer(-c(PMID,year), names_to="keyword",values_to="presence") %>% filter(presence==1)

###################### Plotting ###############################


## trends per year
keywords_per_year <- keywords_AND %>% group_by(year, keyword) %>% summarize(counts=n()) %>% ungroup() %>% arrange(year) %>% group_by(keyword) %>% mutate(cumulative_counts=cumsum(counts))
# the article ID is a line in the pubmed files so it is the foundation of our analysis. We run the distinct function to eliminate possible duplicated lines.

legend_breaks <- unique(keywords_per_year$keyword)
legend_labels <- gsub("oneandtwo",keywords_1_and_2,legend_breaks)

pubmed_keyword_per_year <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=counts, color=keyword), show.legend=T)+
#    scale_x_continuous(breaks=seq(1940,2020,10),limits=c(min(keywords_per_year$year),2020))+
    scale_y_continuous(breaks=pretty_breaks(n=5))+
    scale_color_manual(breaks=legend_breaks,labels=legend_labels, values=c("dodgerblue","darkorange","firebrick1"))+
    ggtitle("Occurrences of keywords per year")+
    theme_bw()+
    theme(legend.position=c(0.1,0.9))
    
ggsave(paste0("../plots/",user_prefix,"_", format(Sys.time(), "%Y-%m-%d_%H-%M"),"_pubmed_keyword_coocurence_per_year.png"), plot = pubmed_keyword_per_year, device = "png", dpi = 150)

# cummulative records of keywords in abstracts over the years
pubmed_keyword_per_year_cumulative <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=cumulative_counts, color=keyword), show.legend=T)+
#    scale_x_continuous(breaks=seq(1940,2020,10),limits=c(min(keywords_per_year$year),2020))+
    scale_y_continuous(breaks=pretty_breaks(n=5))+
    scale_color_manual(breaks=legend_breaks,labels=legend_labels, values=c("dodgerblue","darkorange","firebrick1"))+
    ggtitle("Cumulative occurrences of keywords per year")+
    theme_bw()+
    theme(legend.position=c(0.1,0.9))
ggsave(paste0("../plots/", user_prefix,"_",format(Sys.time(), "%Y-%m-%d_%H-%M"),"_pubmed_keyword_per_year_cumulative.png"), plot = pubmed_keyword_per_year_cumulative, device = "png", dpi = 150)



print("done")
