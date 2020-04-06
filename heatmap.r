#
#
library(tidyverse)

keywords_heatmap_long <- read_delim("keywords_heatmap.txt",delim="\t", col_names=T)


#factor values so they appear in alphabetical order
keywords_heatmap_long$colname <- factor(keywords_heatmap_long$colname, level=unique(keywords_heatmap_long$colname))

keywords_heatmap_long$rowname <- factor(keywords_heatmap_long$rowname, level=unique(keywords_heatmap_long$rowname))

# remove zero's and loops (self occurrence)
keywords_heatmap_long <- keywords_heatmap_long %>% filter(!rowname==colname, count>0)

# elements of the plot
g <- guide_legend("no of abstracts") # legend title and combination of different colors and shapes into one legend
breaks=c(50,500,1000,4000,8000,12000) # the breaks of the legend and point size

#
pubmed_keyword_coocurrence_heatmap <- ggplot()+
  geom_tile(data=keywords_heatmap_long,aes(x=colname, y=rowname,fill=count),alpha=0, show.legend = F)+
  geom_point(data=keywords_heatmap_long,aes(x=colname, y=rowname,colour = count, size =count))  +
  scale_size_binned(name="co-occurrence",range = c(1, 10),breaks=breaks)+
  scale_colour_gradientn(colours =c("gray80","steelblue1","yellowgreen","yellow","orange"),breaks=breaks) +
  guides(colour = g, size = g)+
  ggtitle("Heatmap of co-occurrence of keywords")+
  xlab("") +
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/pubmed_keyword_coocurrence_heatmap.pdf", plot = pubmed_keyword_coocurrence_heatmap, device = "pdf", dpi = 150)
