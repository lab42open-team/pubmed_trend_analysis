#!/usr/bin/Rscript

library(tidyverse)

trends_pubmed <- read_delim("2020-04-01_trends_pubmed.tsv", delim="\t", col_names=F)

colnames(trends_pubmed) <- c("line","PMID","keyword","file","year")

pubmed_keyword_frequency <- ggplot()+
    geom_bar(data=trends_pubmed,aes(keyword))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
ggsave("plots/pubmed_keyword_frequency.pdf", plot = pubmed_keyword_frequency, device = "pdf", dpi = 150)


keywords_per_year <- trends_pubmed %>% distinct(PMID, keyword,year) %>% group_by(year, keyword) %>% summarize(counts=n()) %>% ungroup() %>% arrange(year) %>% mutate(cumulative_counts=cumsum(counts))
# the article ID is a line in the pubmed files so it is the foundation of our analysis. We run the distinct function to eliminate possible duplicated lines.

pubmed_keyword_per_year <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=counts, color=keyword), show.legend=T)+
    scale_x_continuous(breaks=seq(1940,2020,20),limits=c(1940,2020))+
    ggtitle("Occurrences of keywords per year")+
    theme_bw()
    
ggsave("plots/pubmed_keyword_per_year.pdf", plot = pubmed_keyword_per_year, device = "pdf", dpi = 150)

pubmed_keyword_per_year_cumulative <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=cumulative_counts, color=keyword), show.legend=T)+
    ggtitle("Cumulative occurrences of keywords per year")+
    theme_bw()
    
ggsave("plots/pubmed_keyword_per_year_cumulative.pdf", plot = pubmed_keyword_per_year_cumulative, device = "pdf", dpi = 150)

papers_multiple <- trends_pubmed %>% dplyr::distinct(PMID, keyword) %>% group_by(PMID) %>% summarise(n=n()) %>% ungroup() %>% group_by(n) %>% summarise(counts=n())

write_delim(papers_multiple, "papers_multiple.txt", delim="\t")




keywords_coocurence <- trends_pubmed %>% dplyr::distinct(PMID, keyword) %>% group_by(PMID) %>% summarise(keywords_c=toString(keyword)) %>% group_by(keywords_c) %>% summarise(counts=n()) %>% filter(grepl(",",keywords_c)) %>% mutate(b=gsub(".*\\, ","",keywords_c),a=gsub("\\,.*","",keywords_c))

write_delim(keywords_coocurence,"keywords_coocurence.txt",delim="\t")

pubmed_keyword_coocurrence_heatmap <- ggplot(keywords_coocurence,aes(x=a,y=b,fill=counts))+
    geom_tile()+
    scale_fill_gradientn(colours = terrain.colors(3))+
    scale_x_discrete(limits=unique(rev(keywords_coocurence$a)))+
    ggtitle("Heatmap of co-occurrence of keywords")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/pubmed_keyword_coocurrence_heatmap.pdf", plot = pubmed_keyword_coocurrence_heatmap, device = "pdf", dpi = 150)

