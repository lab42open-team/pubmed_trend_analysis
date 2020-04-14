#!/usr/bin/Rscript

library(tidyverse)

keywords_heatmap_long <- read_delim("keywords_heatmap_long.tsv",delim="\t",col_names=T)

# we defined here the diagonal because the raw values don't include them. In addition we need the diagonal seperate from the raw data because we will paint it differently
keywords <- unique(c(keywords_heatmap_long$colname,keywords_heatmap_long$rowname))
diagonal <- tibble(rowname=keywords,colname=keywords,count=-1,jaccard=0) 

## summaries to dynamically set the break points and limits of the plot

summary <- broom::tidy(summary(keywords_heatmap_long$count))

## For the heatmap we need breaks to define the specific points of the legend and limits to ensure that all values will be included in the plot. To create breaks and limits and make them scalable we used the base R functions summary and quantile.
## also quantiles because the raw counts are far apart

quantile <- as.vector(quantile(keywords_heatmap_long$count,probs = c(50,90, 95,98)/100)) # big probabilities because of order of magnitude difference of values 

breaks <- c(floor(summary$minimum),round(quantile[1]),round(quantile[2]),round(quantile[3]),round(quantile[4]),ceiling(summary$maximum)) 

limits=c(min(breaks),max(breaks))

g <- guide_legend("no of abstracts") # legend title and combination of different colors and shapes into one legend

# running the plot
pubmed_keyword_coocurrence_heatmap <- ggplot()+
  geom_tile(data=keywords_heatmap_long,aes(x=colname, y=rowname,fill=count),alpha=0, show.legend = F)+
  geom_point(data=keywords_heatmap_long,aes(x=colname, y=rowname,colour = count, size=count))  +
  scale_size(name="co-occurrence",range = c(0.5, 10),breaks=breaks,limits=limits)+
  scale_colour_gradientn(colours =c("steelblue1","yellowgreen","yellow","goldenrod1","orange"),breaks=breaks,limits=limits) +
  geom_point(data=diagonal,aes(x=colname, y=rowname),colour="lightyellow4",size=1,show.legend = F)+
  scale_x_discrete(position = "top")+
  guides(colour = g, size = g)+
  ggtitle("Heatmap of co-occurrence of keywords")+
  xlab("") +
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0),legend.position = c(.85, .25))

ggsave("count_pubmed_keyword_heatmap.png", plot = pubmed_keyword_coocurrence_heatmap, device = "png", dpi = 150)


# running the log2 plot
### log2 transformation was the best way to gather together counts. Log10 was too aggresive and sqrt too soft.
keywords_heatmap_long$log2 <- log2(keywords_heatmap_long$count)

## summaries to dynamically set the break points and limits of the plot

summary_log <- broom::tidy(summary(keywords_heatmap_long$log2))

breaks_log <- c(floor(summary_log$minimum),round(summary_log$q1),round(summary_log$median),round(summary_log$q3),ceiling(summary_log$maximum))

limits_log=c(min(breaks_log),max(breaks_log))

g_log <- guide_legend("log2(no of abstracts)") # legend title and combination of different colors and shapes into one legend


pubmed_keyword_coocurrence_heatmap <- ggplot()+
  geom_tile(data=keywords_heatmap_long,aes(x=colname, y=rowname,fill=log2),alpha=0, show.legend = F)+
  geom_point(data=keywords_heatmap_long,aes(x=colname, y=rowname,colour = log2, size=log2))  +
  scale_size(name="co-occurrence",range = c(0.5, 10),breaks=breaks_log,limits=limits_log)+
  scale_colour_gradientn(colours =c("steelblue1","yellowgreen","yellow","goldenrod1","orange"),breaks=breaks_log,limits=limits_log) +
  geom_point(data=diagonal,aes(x=colname, y=rowname),colour="lightyellow4",size=1,show.legend = F)+
  scale_x_discrete(position = "top")+
  guides(colour = g_log, size = g_log)+
  ggtitle("Heatmap of co-occurrence of keywords (log2)")+
  xlab("") +
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0),legend.position = c(.85, .25))

ggsave("log_pubmed_keyword_heatmap.png", plot = pubmed_keyword_coocurrence_heatmap, device = "png", dpi = 150)

################################ Jaccard Index ######################################
keywords_n_papers <- read_delim("keywords_n_papers.tsv",delim="\t", col_names=T)

## Jaccard index is the intersection over the union. So we join for each node - keyword the total occurrencies. The join is double because we have two columns of keywords and this way is easier for the calculations

keywords_heatmap_jaccard <- keywords_heatmap_long %>% left_join(keywords_n_papers,c("rowname"="keyword")) %>% left_join(keywords_n_papers,c("colname"="keyword")) %>% mutate(jaccard_index=count/(n_papers.x+n_papers.y-count))

colors_j <- c("steelblue1","yellow","goldenrod1","orange")
limits<- c(0, round(max(keywords_heatmap_jaccard$jaccard_index[keywords_heatmap_jaccard$jaccard_index < 1])*1.07,4) )
# the limits are the maximum value multiplied by >1 so it is bigger than the maximum value

pubmed_jaccard_heatmap <- ggplot()+ 
    geom_tile(data=keywords_heatmap_jaccard,aes(x=colname, y=rowname,fill=jaccard_index),alpha=1, show.legend = T,colour = "white")+
    scale_fill_gradientn(colours =colors_j,limits=limits)+
    geom_tile(data=diagonal,aes(x=colname, y=rowname,fill=count),show.legend = F)+
    scale_x_discrete(position = "top")+
    ggtitle("Heatmap of Jaccard similarity between keywords")+
    xlab("") +
    ylab("")+
    guides(fill = guide_legend(title="Jaccard similarity"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),axis.text.x = element_text(angle = 90, hjust = 0),legend.position = c(.85, .25))

ggsave("pubmed_jaccard_heatmap.png", plot = pubmed_jaccard_heatmap, device = "png", dpi = 150)
