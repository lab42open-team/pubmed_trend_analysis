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
    scale_x_continuous(breaks=seq(1940,2020,10),limits=c(1940,2020))+
    ggtitle("Occurrences of keywords per year")+
    theme_bw()
    
ggsave("plots/pubmed_keyword_per_year.pdf", plot = pubmed_keyword_per_year, device = "pdf", dpi = 150)

pubmed_keyword_per_year_cumulative <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=cumulative_counts, color=keyword), show.legend=T)+
    ggtitle("Cumulative occurrences of keywords per year")+
    theme_bw()
    
ggsave("plots/pubmed_keyword_per_year_cumulative.pdf", plot = pubmed_keyword_per_year_cumulative, device = "pdf", dpi = 150)

# Heatmap

library(Matrix)

papers_keywords_network <- trends_pubmed %>% group_by(PMID, keyword) %>% distinct(PMID, keyword) %>% ungroup()

papers_keywords_matrix <- spMatrix(nrow=length(unique(papers_keywords_network$PMID)),ncol=length(unique(papers_keywords_network$keyword)),i=as.numeric(factor(papers_keywords_network$PMID)),j=as.numeric(factor(papers_keywords_network$keyword)),x = rep(1, length(as.numeric(papers_keywords_network$PMID))))

row.names(papers_keywords_matrix) <- levels(factor(papers_keywords_network$PMID))
colnames(papers_keywords_matrix) <- levels(factor(papers_keywords_network$keyword))

keywords_heatmap <- tcrossprod(t(papers_keywords_matrix))

keywords_heatmap_long <- as.data.frame(as.matrix(keywords_heatmap)) %>% rownames_to_column() %>% pivot_longer(-rowname,names_to="colname",values_to="count" ) # %>% filter(count>0)

write_delim(keywords_heatmap_long, "keywords_heatmap.txt", delim="\t")

keywords_heatmap_long$colname <- factor(keywords_heatmap_long$colname, level=unique(keywords_heatmap_long$colname))
keywords_heatmap_long$rowname <- factor(keywords_heatmap_long$rowname, level=unique(keywords_heatmap_long$rowname))

pubmed_keyword_coocurrence_heatmap <- ggplot(keywords_heatmap_long,aes(x=colname,y=rowname,fill=count))+
    #geom_tile()+
    geom_point(aes(colour = count, size =count))  +    ## geom_point for circle illusion
    #scale_color_gradient(low = "yellow",high = "red")+       ## color of the corresponding aes
    scale_size(range = c(1, 10))+
    #scale_color_brewer(palette = "PuOr")+
    ggtitle("Heatmap of co-occurrence of keywords")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/pubmed_keyword_coocurrence_heatmap.pdf", plot = pubmed_keyword_coocurrence_heatmap, device = "pdf", dpi = 150)

