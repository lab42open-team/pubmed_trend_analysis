#!/usr/bin/Rscript

# In this script we take the data from trends extraction from pubmed and transform them to create plots.
# The data file contains 5 fields, the line of the original file, the id of the article - PMID, the keyword that was matched, the file that contains the keyword and the year of the publication.
# Three types of plots are created, bar plot, yearly line plot and a heatmap of co-occurrances
#
# Author:   Savvas Paragkamian (s.paragkamian@hcmr.gr)
#           Institute of Marine Biology Biotechnology and Aquaculture (IMBBC)
#           Hellenic Centre for Marine Research (HCMR)
#
# Created:  05/06/2020
# License:  GNU GPLv3 license
# Developer : Savvas Paragkamian
#
# Change the n_papers_pubmed when pubmed is refreshed in order to estimate the random chance of keyword occurrance.
#
# Run
# Rscript trends_plots.r ../demo/total_data.tsv "test_ccmri" ../demo/keywords.txt
#
# Packages

suppressPackageStartupMessages({
    library(tidyverse) # version tidyverse 1.3.0, used packages from tidyverse are readr_1.3.1, ggplot2_3.3.0, dplyr_0.8.5
    library(Matrix) # version Matrix 1.2-15
    library(igraph)
    library(ggraph)
    library(tidygraph)
})

## file loading. Plotting takes two arguments,
## 1) file name
## 2) prefix for the plots names
## 3) user keywords file
#args <- c("../demo/total_data.tsv", "test_ccmri", "../demo/keywords.txt")

n_papers_pubmed <- 36304541 # number of unique papers in PubMed. Must be updated when PubMed will be refreshed!!!!!!!

args <- commandArgs(trailingOnly=TRUE)
user_prefix <- args[2]

trends_pubmed <- read_delim(args[1], delim="\t", col_names=F,col_types = cols())

colnames(trends_pubmed) <- c("year","PMID","synonym")

trends_categories <- read_delim(args[3], delim="\t", col_names=F,col_types = cols()) |> arrange(X3)
colnames(trends_categories) <- c("synonym","keyword","category")

trends_categories_only <- trends_categories |> distinct(keyword,category)
## filter only the keywords that are listed in the trends_categories and then
## join them to keep the general categories. Also remove the the synonyms to
## keep only the unique number of PMIDs per keyword.

trends_pubmed <- trends_pubmed |>
    filter(synonym %in% trends_categories$synonym) |> 
    dplyr::left_join(trends_categories, by=c("synonym"="synonym")) |> 
    dplyr::distinct(year,PMID,keyword,category)

## bar plot of keyword frequencies

trends_counts <- trends_pubmed |>
    distinct(PMID,keyword,category) |>
    group_by(keyword,category) |>
    summarise(counts=n(), .groups="keep")

trends_counts$keyword <- fct_reorder(trends_counts$keyword,trends_counts$category)

pubmed_keyword_frequency <- ggplot()+
    geom_col(data=trends_counts,
             aes(x=keyword,y=counts, fill=category),
             width=0.8)+
    geom_text(data=trends_counts,
              aes(x=keyword,y=counts,label=counts),
              vjust=-0.4,
              color="grey70",
              size=2)+
    scale_y_continuous(name="Number of abstracts",
                       limits=c(0,max(trends_counts$counts)),
                       n.breaks=10)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor=element_blank())

ggsave(paste0("../plots/",user_prefix,"_", format(Sys.time(), "%Y-%m-%d_%H-%M"),"_pubmed_keyword_frequency.png"),
       plot = pubmed_keyword_frequency,
       device = "png",
       dpi = 150)


## trends per year
keywords_per_year <- trends_pubmed %>% 
    distinct(PMID, keyword,category,year) %>% 
    group_by(year, keyword,category) %>% 
    summarize(counts=n()) %>% 
    ungroup() %>% 
    arrange(year) %>% 
    group_by(keyword,category) %>% 
    mutate(cumulative_counts=cumsum(counts)) %>% 
    ungroup() %>% 
    mutate(keyword=fct_reorder(keyword,category, .desc=TRUE)) %>% 
    mutate(count_bin=cut(counts, breaks=c(0,10, 50, 100, 500, 2000, max(counts,na.rm=T)),
                         labels=c("1-10","10-50", "50-100", "100-500","500-2000","2000<")))

keywords_per_year <- trends_pubmed |>
    distinct(PMID, keyword,category, year) |>
    group_by(year, keyword, category) |>
    summarize(counts=n(),.groups="keep") |>
    ungroup() |>
    arrange(year) |>
    group_by(keyword, category) |>
    mutate(cumulative_counts=cumsum(counts)) |>
    ungroup() |> 
#    mutate(keyword=fct_reorder(keyword, category, .desc=TRUE)) |>
    mutate(count_bin=cut(counts, breaks=c(0,10, 50, 100, 500, 2000, max(counts,na.rm=T)),
                         labels=c("1-10","10-50", "50-100", "100-500","500-2000","2000<")))


# change the order to descreasing to appear with alphabetical order
keywords_per_year$keyword <- factor(keywords_per_year$keyword,
                                    levels=unique(keywords_per_year$keyword[order(keywords_per_year$category,
                                                                                  keywords_per_year$keyword,
                                                                                  decreasing=T)]))

# the article ID is a line in the pubmed files so it is the foundation of our analysis. We run the distinct function to eliminate possible duplicated lines.

pubmed_keyword_per_year <- ggplot()+
    geom_line(data=keywords_per_year,
              aes(x=year,y=counts, color=keyword),
              show.legend=T)+
    scale_x_continuous(n.breaks=6,
                       limits=c(min(keywords_per_year$year,na.rm=T),max(keywords_per_year$year, na.rm=T)))+
    scale_y_continuous(n.breaks=5)+
    ggtitle("Occurrences of keywords per year")+
    theme_bw()
    
ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(),
                     "%Y-%m-%d_%H-%M"),
              "_pubmed_keyword_per_year.png"),
       plot = pubmed_keyword_per_year,
       device = "png",
       dpi = 150)

# cummulative records of keywords in abstracts over the years
pubmed_keyword_per_year_cumulative <- ggplot()+
    geom_line(data=keywords_per_year,
              aes(x=year,y=cumulative_counts, color=keyword),
              show.legend=T)+
    scale_x_continuous(n.breaks=6,
                       limits=c(min(keywords_per_year$year,na.rm=T),max(keywords_per_year$year, na.rm=T)))+
    ggtitle("Cumulative occurrences of keywords per year")+
    theme_bw()
   
ggsave(paste0("../plots/",
              user_prefix,"_",format(Sys.time(), "%Y-%m-%d_%H-%M"),
              "_pubmed_keyword_per_year_cumulative.png"),
       plot = pubmed_keyword_per_year_cumulative,
       device = "png",
       dpi = 150)

################################## Heatmaps #############################

############################## timeline heatmap #########################

pubmed_keyword_per_year_heatmap <- ggplot()+
    geom_tile(data=keywords_per_year,
              aes(x=year,y=keyword,
                  fill=count_bin,
                  height=1,
                  width=1),
              color="white",
              show.legend=T)+
    scale_x_continuous(n.breaks=10)+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35))+
    scale_fill_manual(values=c("#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#4a1486"))+
    ylab("")+
    xlab("")+
    coord_equal()+
    theme_bw()+
    guides(fill=guide_legend(title="# of abstracts"))+
    theme(panel.border=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor=element_blank(),
          legend.position="bottom",
          legend.direction="horizontal",
          legend.key.height=unit(0.2,"cm"),
          legend.key.width = unit(0.8,"cm"),
          axis.text.y = element_text(size = 9),
          plot.margin=margin(0,0,0,0,"cm"))
   
ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(), "%Y%m%d%H%M"),
              "_key_time_heatmap.png"),
       plot=pubmed_keyword_per_year_heatmap,
       width = 45,
       height = 15,
       units='cm',
       device = "png",
       dpi = 300)

##### with facet of categories
key_time_heatmap_facet <- ggplot(data=keywords_per_year)+
    geom_tile(aes(x=year,y=keyword, fill=count_bin),size=0.2,color="white", show.legend=T)+
    scale_x_continuous(expand=c(0,0),limits=c(1960,2022),breaks=seq(1960,2030,5))+
#    scale_y_discrete(expand=c(0,0),limits = rev(levels(keywords_per_year$keyword)))+
    scale_fill_manual(values=c("#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#4a1486"))+
    ylab("")+
    xlab("")+
    theme_bw()+
    guides(fill=guide_legend(title="# of abstracts"))+
    theme(plot.background=element_blank(),
          panel.border=element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor=element_blank(), 
          legend.position="right",legend.direction="vertical",
          legend.key.height=unit(0.8,"cm"),
          legend.key.width = unit(0.2,"cm"),
          plot.margin=margin(0.1,0.2,0.2,0.2,"cm"),
          axis.text=element_text(size=10),
          strip.text.y = element_text(angle =0,size=9),
          strip.background = element_rect(colour = "white", fill = "white") ,
          aspect.ratio = 1) +
    facet_grid(rows = vars(category),space = "free",scales = "free_y",labeller = labeller(category = label_wrap_gen(9)))
   
ggsave(paste0("../plots/", user_prefix,"_",format(Sys.time(), "%Y%m%d%H%M"),"_key_time_heatmap_facet.png"), 
       plot =key_time_heatmap_facet ,height=12, width = 30, units='cm',device = "png", dpi = 300)



# create the edglist of keywords and PMID's
papers_keywords_network <- trends_pubmed |>
    group_by(PMID, keyword) |> 
    distinct(PMID, keyword) |> 
    ungroup()

keywords_n_papers <- papers_keywords_network |>
    group_by(keyword) |>
    summarise(n_papers=n()) |> 
    mutate(freq=n_papers/n_papers_pubmed) 


# create a matrix class spMatrix (handles better sparse matrices) to do inverse table multiplication
papers_keywords_matrix <- spMatrix(nrow=length(unique(papers_keywords_network$PMID)),
                                   ncol=length(unique(papers_keywords_network$keyword)),
                                   i=as.numeric(factor(papers_keywords_network$PMID)),
                                   j=as.numeric(factor(papers_keywords_network$keyword)),
                                   x = rep(1, length(as.numeric(papers_keywords_network$PMID))))

row.names(papers_keywords_matrix) <- levels(factor(papers_keywords_network$PMID))
colnames(papers_keywords_matrix) <- levels(factor(papers_keywords_network$keyword))

# with the inverse cross product we do the projection of the edgelist to keywords in order to calculate how many times keyword pairs appear together in abstracts.
keywords_heatmap <- tcrossprod(t(papers_keywords_matrix))

# becaue the matrix is summetric we keep the triangle
keywords_heatmap[upper.tri(keywords_heatmap)] <- 0

keywords_heatmap <- as.data.frame(as.matrix(keywords_heatmap))
#write_delim(keywords_heatmap,"keywords_heatmap.tsv",delim="\t")

# transform to long format for plotting and remove zero's and NA's and assign -1 to loops (self occurrence)
keywords_heatmap_long <- as.data.frame(as.matrix(keywords_heatmap)) |> 
    rownames_to_column() |> 
    pivot_longer(-rowname,names_to="colname",values_to="count" ) |> 
    filter(count!=0,colname!=rowname) |> 
    na.omit()

######################################### Network analysis ###########################################

coword_graph <- graph_from_adjacency_matrix(as.matrix(keywords_heatmap),weighted=T, mode="undirected")

coword_graph <- simplify(coword_graph,remove.loops=T)

### Centralities calculation

V(coword_graph)$Degree <- degree(coword_graph)
V(coword_graph)$betweenness <- betweenness(coword_graph, weights=NA)
V(coword_graph)$betweenness_w <- betweenness(coword_graph,weights=E(coword_graph)$weight)
V(coword_graph)$closeness <- closeness(coword_graph)
V(coword_graph)$transitivity <- transitivity(coword_graph,type="local", weights=E(coword_graph)$weight)
V(coword_graph)$page_rank <- page_rank(coword_graph)$vector

E(coword_graph)$link_width <- as_tibble(round(log(E(coword_graph)$weight),2)) |> 
    mutate(value=if_else(value==0,0.5,value)) |>
    mutate(value=value/max(value)) |>
    pull(value)

layout_circle <- layout_in_circle(coword_graph)
layout_fr <- layout_with_fr(coword_graph)
layout_nice <- layout_nicely(coword_graph)


coword_graph_tidy <- as_tbl_graph(coword_graph)

p <- ggraph(coword_graph_tidy,layout = 'linear', circular = TRUE) +
        geom_edge_link(aes(edge_color=weight))+
        scale_edge_color_continuous(low="gray90", high="gray20")+
#        scale_edge_colour_continuous(low = "white", high = "black" na.value = "grey50")+
#        geom_node_text(aes(label = name), , repel = TRUE)+
        geom_node_point(aes(size=Degree), colour="darkgoldenrod")+
        coord_fixed()+
        theme_graph()

p <- p + geom_node_text(aes(label = name),
                        nudge_x = p$data$x * .15,
                        nudge_y = p$data$y * .15)

ggsave(paste0("../plots/", user_prefix,"_", format(Sys.time(), "%Y-%m-%d_%H-%M"),"_pubmed_keyword_network.png"),
       plot = p,
       width = 25,
       height = 25,
       units='cm' ,
       device = "png",
       dpi = 300)

#png(file=paste0("../plots/", user_prefix,"_", format(Sys.time(), "%Y-%m-%d_%H-%M"),"_pubmed_keyword_network_circular.png"), width = 25, height = 25, units='cm',res=300)
#plot(coword_graph, vertex.label=V(coword_graph)$name, vertex.size=V(coword_graph)$degree, vertex.label.cex=.9,vertex.color= "lightblue", vertex.frame.color="white", edge.size=E(coword_graph)$weights, layout=layout_circle)
#dev.off()

######################################### Heatmap plotting ###########################################

# elements of the plot

# we defined here the diagonal because the raw values don't include them. In addition we need the diagonal seperate from the raw data because we will paint it differently
keywords <- unique(c(keywords_heatmap_long$colname,keywords_heatmap_long$rowname))
diagonal <- tibble(rowname=keywords,colname=keywords,count=-1,jaccard=0) 

## summaries to dynamically set the break points and limits of the plot

summary <- summary(keywords_heatmap_long$count)

## For the heatmap we need breaks to define the specific points of the legend and limits to ensure that all values will be included in the plot. To create breaks and limits and make them scalable we used the base R functions summary and quantile.
## also quantiles because the raw counts are far apart

quantile <- as.vector(quantile(keywords_heatmap_long$count,probs = c(50,90, 95,98)/100)) # big probabilities because of order of magnitude difference of values 

breaks <- c(floor(min(summary)),
            round(quantile[1]),
            round(quantile[2]),
            round(quantile[3]),
            round(quantile[4]),
            ceiling(max(summary))) 

limits=c(min(breaks),max(breaks))

g <- guide_legend("no of abstracts") # legend title and combination of different colors and shapes into one legend

# running the plot
pubmed_keyword_coocurrence_heatmap <- ggplot()+
  geom_tile(data=keywords_heatmap_long,
            aes(x=colname, 
                y=rowname,
                fill=count),
            alpha=0,
            show.legend = F)+
  geom_point(data=keywords_heatmap_long,
             aes(x=colname,
                 y=rowname,
                 colour = count,
                 size=count))  +
  scale_size(name="co-occurrence",
             range = c(0.5, 10),
             breaks=breaks,limits=limits)+
  scale_colour_gradientn(colours =c("steelblue1","yellowgreen","yellow","goldenrod1","orange"),
                         breaks=breaks,
                         limits=limits) +
  geom_point(data=diagonal,
             aes(x=colname, y=rowname),
             colour="lightyellow4",
             size=1,
             show.legend = F)+
  scale_x_discrete(position = "top")+
  guides(colour = g, size = g)+
  ggtitle("Heatmap of co-occurrence of keywords")+
  xlab("") +
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0),legend.position = c(.85, .25))

ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(),
                     "%Y-%m-%d_%H-%M"),
              "_pubmed_keyword_heatmap.png"),
       plot = pubmed_keyword_coocurrence_heatmap,
       device = "png",
       dpi = 150)

#write_delim(keywords_heatmap_long,"heatmap_data.txt", delim="\t")

####################################### running the log2 heatmap ###################################################
### log2 transformation was the best way to gather together counts. Log10 was too aggresive and sqrt too soft.
keywords_heatmap_long$log2 <- log2(keywords_heatmap_long$count)

## summaries to dynamically set the break points and limits of the plot

summary_log <- summary(keywords_heatmap_long$log2)

breaks_log <- c(floor(min(summary)),
                round(summary_log[2]),
                round(summary_log[3]),
                round(summary_log[4]),
                ceiling(max(summary_log)))

limits_log=c(min(breaks_log),max(breaks_log))

g_log <- guide_legend("log2(no of abstracts)") # legend title and combination of different colors and shapes into one legend


pubmed_keyword_coocurrence_heatmap <- ggplot()+
  geom_tile(data=keywords_heatmap_long,
            aes(x=colname,
                y=rowname,
                fill=log2),
            alpha=0,
            show.legend = F)+
  geom_point(data=keywords_heatmap_long,
             aes(x=colname,
                 y=rowname,
                 colour = log2,
                 size=log2))  +
  scale_size(name="co-occurrence",
             range = c(0.5, 10),
             breaks=breaks_log,
             limits=limits_log)+
  scale_colour_gradientn(colours =c("steelblue1","yellowgreen","yellow","goldenrod1","orange"),
                         breaks=breaks_log,
                         limits=limits_log) +
  geom_point(data=diagonal,
             aes(x=colname,
                 y=rowname),
             colour="lightyellow4",
             size=1,
             show.legend = F)+
  scale_x_discrete(position = "top")+
  guides(colour = g_log, size = g_log)+
  ggtitle("Heatmap of co-occurrence of keywords (log2)")+
  xlab("") +
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        legend.position = c(.85, .25))

ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(),
                     "%Y-%m-%d_%H-%M"),
              "_log_pubmed_keyword_heatmap.png"),
       plot = pubmed_keyword_coocurrence_heatmap,
       device = "png",
       dpi = 150)

################################ Jaccard Index ######################################

## Jaccard index is the intersection over the union. So we join for each node - keyword the total occurrencies. The join is double because we have two columns of keywords and this way is easier for the calculations

keywords_heatmap_jaccard <- keywords_heatmap_long |>
    left_join(keywords_n_papers,c("rowname"="keyword")) |> 
    left_join(keywords_n_papers,c("colname"="keyword")) |> 
    mutate(jaccard_index=count/(n_papers.x+n_papers.y-count),
           random_expectation=(count/n_papers_pubmed)/(freq.x*freq.y))

colors_j <- c("steelblue1","yellow","goldenrod1","orange")
limits<- c(0,
           round(max(keywords_heatmap_jaccard$jaccard_index[keywords_heatmap_jaccard$jaccard_index < 1])*1.07,4) )
# the limits are the maximum value multiplied by >1 so it is bigger than the maximum value

pubmed_jaccard_heatmap <- ggplot()+ 
    geom_tile(data=keywords_heatmap_jaccard,
              aes(x=colname,
                  y=rowname,
                  fill=jaccard_index),
              alpha=1,
              show.legend = T,
              colour = "white")+
    scale_fill_gradientn(colours =colors_j,limits=limits)+
    geom_tile(data=diagonal,
              aes(x=colname,
                  y=rowname,
                  fill=count),
              show.legend = F)+
    scale_x_discrete(position = "top")+
    ggtitle("Heatmap of Jaccard similarity between keywords")+
    xlab("") +
    ylab("")+
    guides(fill = guide_legend(title="Jaccard similarity"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0),
          legend.position = c(.85, .25))

ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(), "%Y-%m-%d_%H-%M"),
              "_pubmed_jaccard_heatmap.png"),
       plot = pubmed_jaccard_heatmap,
       device = "png",
       dpi = 150)

################################ Actual vs Random  ######################################

# todo : impove colors, below 1 must be shown. rescale .

pubmed_random_heatmap <- ggplot()+ 
    geom_tile(data=keywords_heatmap_jaccard,
              aes(x=colname,
                  y=rowname,
                  fill=random_expectation),
              alpha=1,
              show.legend = T,
              colour = "white")+
    geom_tile(data=diagonal,
              aes(x=colname,
                  y=rowname,
                  fill=count),
              show.legend = F)+
    scale_fill_gradient2(low = "steelblue1",
                         mid = "white",
                         high = "orange",
                         midpoint = 1,
                         na.value = "grey50")+
    scale_x_discrete(position = "top")+
    ggtitle("Heatmap of Actual against random expectation co-occurrence between keywords")+
    xlab("") +
    ylab("")+
    guides(fill = guide_legend(title="Times above random"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0),
          legend.position = c(.8, .25))

ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(), "%Y-%m-%d_%H-%M"),
              "_pubmed_random_heatmap.png"),
       plot = pubmed_random_heatmap,
       device = "png",
       dpi = 150)
