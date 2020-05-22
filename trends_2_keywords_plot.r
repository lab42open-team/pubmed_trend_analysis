#!/usr/bin/Rscript

library(tidyverse) # version tidyverse 1.3.0, used packages from tidyverse are readr_1.3.1, ggplot2_3.3.0, dplyr_0.8.5
library(Matrix) # version Matrix 1.2-15

### This script creates plot to compare coocurrences between two keywords.
### The user feeds the script with 4 arguments,
### 1. data file name
### 2. keyword 1
### 3. keyword 2 
### 4. plot prefix


# load the user arguments
args <- commandArgs(trailingOnly=TRUE)
keyword1 <- args[2]
keyword2 <- args[3]
keywords_1_and_2 <- paste0(keyword1,"_AND_",keyword2)
user_prefix <- args[4]

## file loading
trends_pubmed <- read_delim(args[1], delim="\t", col_names=F)

colnames(trends_pubmed) <- c("line","PMID","keyword","file","year")

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
    scale_x_continuous(breaks=seq(1940,2020,10),limits=c(1940,2020))+
    scale_color_manual(breaks=legend_breaks,labels=legend_labels, values=c("dodgerblue","darkorange","firebrick1"))+
    ggtitle("Occurrences of keywords per year")+
    theme_bw()
    
ggsave(paste0("plots/",user_prefix,"_", format(Sys.time(), "%Y-%m-%d_%H-%M"),"_pubmed_keyword_coocurence_per_year.png"), plot = pubmed_keyword_per_year, device = "png", dpi = 150)



print("done")
