#!/usr/bin/Rscript

# In this script we take the data from trends extraction from pubmed and transform them to create plots.
# The data file contains 5 fields, the line of the original file, the id of the article - PMID, the keyword that was matched, the file that contains the keyword and the year of the publication.
# Three types of plots are created, bar plot, yearly line plot and a heatmap of co-occurrances
# 
library(tidyverse) # version tidyverse 1.3.0, used packages from tidyverse are readr_1.3.1, ggplot2_3.3.0, dplyr_0.8.5
library(Matrix) # version Matrix 1.2-15

## file loading
args = commandArgs(trailingOnly=TRUE)

trends_pubmed <- read_delim(args[1], delim="\t", col_names=F)

colnames(trends_pubmed) <- c("line","PMID","keyword","file","year")

## bar plot of keyword frequencies
pubmed_keyword_frequency <- ggplot()+
    geom_bar(data=trends_pubmed,aes(keyword))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
ggsave("plots/pubmed_keyword_frequency.png", plot = pubmed_keyword_frequency, device = "png", dpi = 150)

## trends per year
keywords_per_year <- trends_pubmed %>% distinct(PMID, keyword,year) %>% group_by(year, keyword) %>% summarize(counts=n()) %>% ungroup() %>% arrange(year) %>% mutate(cumulative_counts=cumsum(counts))
# the article ID is a line in the pubmed files so it is the foundation of our analysis. We run the distinct function to eliminate possible duplicated lines.

pubmed_keyword_per_year <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=counts, color=keyword), show.legend=T)+
    scale_x_continuous(breaks=seq(1940,2020,10),limits=c(1940,2020))+
    ggtitle("Occurrences of keywords per year")+
    theme_bw()
    
ggsave("plots/pubmed_keyword_per_year.png", plot = pubmed_keyword_per_year, device = "png", dpi = 150)

# cummulative records of keywords in abstracts over the years
pubmed_keyword_per_year_cumulative <- ggplot()+
    geom_line(data=keywords_per_year, aes(x=year,y=cumulative_counts, color=keyword), show.legend=T)+
    scale_x_continuous(breaks=seq(1940,2020,10),limits=c(1940,2020))+
    ggtitle("Cumulative occurrences of keywords per year")+
    theme_bw()
   
ggsave("plots/pubmed_keyword_per_year_cumulative.png", plot = pubmed_keyword_per_year_cumulative, device = "png", dpi = 150)

##################### Heatmap #######################

# create the edglist of keywords and PMID's
papers_keywords_network <- trends_pubmed %>% group_by(PMID, keyword) %>% distinct(PMID, keyword) %>% ungroup()

# create a matrix class spMatrix (handles better sparse matrices) to do inverse table multiplication
papers_keywords_matrix <- spMatrix(nrow=length(unique(papers_keywords_network$PMID)),ncol=length(unique(papers_keywords_network$keyword)),i=as.numeric(factor(papers_keywords_network$PMID)),j=as.numeric(factor(papers_keywords_network$keyword)),x = rep(1, length(as.numeric(papers_keywords_network$PMID))))

row.names(papers_keywords_matrix) <- levels(factor(papers_keywords_network$PMID))
colnames(papers_keywords_matrix) <- levels(factor(papers_keywords_network$keyword))

# with the inverse cross product we do the projection of the edgelist to keywords in order to calculate how many times keyword pairs appear together in abstracts.
keywords_heatmap <- tcrossprod(t(papers_keywords_matrix))

# because the matrix is summetric we keep the triangle
keywords_heatmap[upper.tri(keywords_heatmap)] <- NA

# transform to long format for plotting and remove zero's and NA's and assign -1 to loops (self occurrence)
keywords_heatmap_long <- as.data.frame(as.matrix(keywords_heatmap)) %>% rownames_to_column() %>% pivot_longer(-rowname,names_to="colname",values_to="count" ) %>% mutate(count=if_else(colname==rowname,-1,count)) %>% filter(count!=0) %>% na.omit()

#factor values so they appear in alphabetical order
keywords_heatmap_long$colname <- factor(keywords_heatmap_long$colname, level=unique(keywords_heatmap_long$colname))

keywords_heatmap_long$rowname <- factor(keywords_heatmap_long$rowname, level=unique(keywords_heatmap_long$rowname))

# elements of the plot
## only the loops in order to color the heatmap diagonal differently
diagonal <- keywords_heatmap_long %>% filter(count==-1)

g <- guide_legend("no of abstracts") # legend title and combination of different colors and shapes into one legend
breaks=c(1,50,500,1000,4000,8000,12000) # the breaks of the legend and point size

# running the plot
pubmed_keyword_coocurrence_heatmap <- ggplot()+
  geom_tile(data=keywords_heatmap_long,aes(x=colname, y=rowname,fill=count),alpha=0, show.legend = F)+
  geom_point(data=keywords_heatmap_long,aes(x=colname, y=rowname,colour = count, size =count))  +
  scale_size_binned(name="co-occurrence",range = c(0.5, 10),breaks=breaks)+
  scale_colour_gradientn(colours =c("gray80","steelblue1","yellowgreen","yellow","goldenrod1","orange"),breaks=breaks) +
  geom_point(data=diagonal,aes(x=colname, y=rowname),colour="lightyellow4",size=1,show.legend = F)+
  scale_x_discrete(position = "top")+
  guides(colour = g, size = g)+
  ggtitle("Heatmap of co-occurrence of keywords")+
  xlab("") +
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0),legend.position = c(.85, .25))

ggsave("plots/pubmed_keyword_heatmap.png", plot = pubmed_keyword_coocurrence_heatmap, device = "png", dpi = 150)


