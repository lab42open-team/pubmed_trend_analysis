#!/usr/bin/Rscript

library(tidyverse)

trends_pubmed <- read_delim("review_2020-04-14_11-50_trends_pubmed.tsv",delim="\t",col_names=F)
colnames(trends_pubmed) <- c("line","PMID","keyword","file","year")

trends_counts <- trends_pubmed %>% distinct(PMID,keyword) %>% group_by(keyword) %>% summarise(counts=n())

## bar plot of keyword frequencies
pubmed_keyword_frequency <- ggplot()+
    geom_col(data=trends_counts,aes(x=keyword,y=counts),width=0.8)+
    geom_text(data=trends_counts,aes(x=keyword,y=counts,label=counts), vjust=-0.4, color="grey70", size=2)+
    scale_y_continuous(name="Number of papers",limits=c(0,max(trends_counts$counts)),n.breaks=10)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid.major.x = element_blank() ,panel.grid.minor=element_blank())
    
ggsave("bar_pubmed_keyword_frequency.png", plot = pubmed_keyword_frequency, device = "png", dpi = 150)


## trends per year
keywords_per_year <- trends_pubmed %>% distinct(PMID, keyword,year) %>% group_by(year, keyword) %>% summarize(counts=n()) %>% ungroup() %>% arrange(year) %>% group_by(keyword) %>% mutate(cumulative_counts=cumsum(counts))
# the article ID is a line in the pubmed files so it is the foundation of our analysis. We run the distinct function to eliminate possible duplicated lines.

pubmed_keyword_per_year <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=counts, color=keyword), show.legend=T)+
    scale_x_continuous(breaks=seq(1960,2020,10),limits=c(1960,2020))+
    ggtitle("Occurrences of keywords per year")+
    scale_colour_viridis_d(option = "inferno")+
    theme_bw()
    
ggsave("year_color_pubmed_keyword_per_year.png", plot = pubmed_keyword_per_year, device = "png", dpi = 150)

# cummulative records of keywords in abstracts over the years
pubmed_keyword_per_year_cumulative <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=cumulative_counts, color=keyword), show.legend=T)+
    scale_x_continuous(breaks=seq(1960,2020,10),limits=c(1960,2020))+
    ggtitle("Cumulative occurrences of keywords per year")+
    theme_bw()
   
ggsave("c_pubmed_keyword_per_year_cumulative.png", plot = pubmed_keyword_per_year_cumulative, device = "png", dpi = 150)
