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
    library(pheatmap)
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

ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(), "%Y-%m-%d_%H-%M"),
              "_pubmed_keyword_frequency.png"),
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
#    coord_equal()+
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

key_time_heatmap_facet <- pubmed_keyword_per_year_heatmap + 
    facet_grid(rows = vars(category),
               space = "free",
               scales = "free_y",
               labeller = labeller(category = label_wrap_gen(9)))
   
ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(), "%Y%m%d%H%M"),
              "_key_time_heatmap_facet.png"), 
       plot=key_time_heatmap_facet,
       height=12,
       width = 30,
       units='cm',
       device = "png",
       dpi = 300)

#################################### Co-occerrence of keywords #####################################

# create the edglist of keywords and PMID's
papers_keywords_network <- trends_pubmed |>
    group_by(PMID, keyword) |> 
    distinct(PMID, keyword) |> 
    ungroup() |>
    left_join(trends_categories_only, by=c("keyword"="keyword"))

keywords_n_papers <- papers_keywords_network |>
    group_by(keyword) |>
    summarise(n_papers=n()) |> 
    mutate(freq=n_papers/n_papers_pubmed) 

# very important for the correct order of the keywords based on categories.
papers_keywords_network$keyword <- factor(papers_keywords_network$keyword,
                                          levels=unique(papers_keywords_network$keyword[order(papers_keywords_network$category,
                                                                                              papers_keywords_network$keyword)]))

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

colnames(keywords_heatmap_long) <- c("from","to","count")

keywords_heatmap_long$count_bin <- cut(keywords_heatmap_long$count, 
                                       breaks=c(0,5,50, 100, 500, 800, 1000),
                                       labels=c("1-5","5-50", "50-100", "100-500","500-800","800<"))

keywords_heatmap_long$count_n <- cut(keywords_heatmap_long$count, breaks=10)

# assign the order levels of the count_bin
keywords_heatmap_long$count_bin <- factor(as.character(keywords_heatmap_long$count_bin),
                                          levels=rev(levels(keywords_heatmap_long$count_bin)))

## add the categories for the keywords
keywords_heatmap_long <- keywords_heatmap_long %>% 
    left_join(trends_categories_only, by=c("from"="keyword")) %>% 
    left_join(trends_categories_only, by=c("to"="keyword"))
 

######################################### Network analysis ###########################################

colnames(trends_counts) <- c("name","category","abstracts" )

keywords_heatmap_long_net <- keywords_heatmap_long |>
    filter(count!=0, from!=to)

coword_graph <- graph_from_data_frame(keywords_heatmap_long_net, directed=FALSE, vertices=trends_counts)

V(coword_graph)$color <- V(coword_graph)$category
### Centralities calculation

V(coword_graph)$Degree <- degree(coword_graph)
V(coword_graph)$betweenness <- betweenness(coword_graph, weights=NA)
V(coword_graph)$betweenness_w <- betweenness(coword_graph,weights=E(coword_graph)$weight)
V(coword_graph)$closeness <- closeness(coword_graph)
V(coword_graph)$transitivity <- transitivity(coword_graph,type="local", weights=E(coword_graph)$weight)
V(coword_graph)$page_rank <- page_rank(coword_graph)$vector

E(coword_graph)$link_width <- as_tibble(round(log(E(coword_graph)$count),2)) %>% mutate(value=if_else(value==0,0.5,value)) %>% mutate(value=value/max(value)) %>% pull(value)

layout_circle <- layout_in_circle(coword_graph)
layout_fr <- layout_with_fr(coword_graph)
layout_nice <- layout_nicely(coword_graph)

coword_graph_tidy <- as_tbl_graph(coword_graph)

#p <- ggraph(coword_graph_tidy,layout = 'linear', circular = TRUE) +
#p <- p + geom_node_text(aes(label = name),repel = TRUE)#, nudge_x = p$data$x * .15, nudge_y = p$data$y * 25)
#p <- ggraph(coword_graph_tidy,layout = "centrality",cent = graph.strength(coword_graph_tidy)) +

p <- ggraph(coword_graph_tidy,layout = 'stress') +
        geom_edge_link(aes(edge_color=count,
                           edge_width = count))+
        geom_node_point(aes(size=Degree,
                            shape=category,
                            color=category))+
        scale_edge_color_continuous(low="#bdbdbd",
                                    high="#636363")+
        geom_node_text(aes(color=category,
                           label = name),
                       fontface = "bold",
                       nudge_y = 0.1,
                       check_overlap = TRUE,
                       show.legend=FALSE)+
#        scale_color_manual(values=c("gut microbiome"="#e66101","flux"="#5e3c99","proteomics"="#bababa","medical"="#fdb863"))+
        scale_size(range = c(0.5,6),
                   breaks=c(5,15,25))+
        scale_edge_width(range = c(0.1,2))+
        guides(edge_color= guide_legend("# of co-occurrence\nin abstracts",order = 3),
               edge_width=guide_legend("# of co-occurrence\nin abstracts",order = 3), 
               shape=guide_legend("Compound",order = 1),
               color=guide_legend("Compound",order = 1), 
               size=guide_legend("Keywords co-occurrences\n(degree)",
                                 order = 2,
                                 override.aes = list(color="gray50")))+
        theme_graph()+
        coord_cartesian(clip = "off")+
        theme(legend.justification = "top",
              legend.box="horizontal",
              legend.direction= "vertical",
              legend.margin=margin(unit(0.4,"cm")),
              legend.position = "bottom",
              legend.spacing.x=unit(1.5,"cm"),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 14),
              legend.key.size = unit(0.8,"cm"))

ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(),"%Y%m%d%H%M"),"_net.png"),
       plot = p,
       width = 25,
       height = 25,
       units='cm',
       device = "png",
       dpi = 300)


######################################### Heatmap plotting ###########################################

keywords_heatmap_wide <- as.data.frame(as.matrix(keywords_heatmap)) 

diag(keywords_heatmap_wide) <- 0

png(paste0("../plots/",
           user_prefix,"_",
           format(Sys.time(),"%Y-%m-%d_%H-%M"),
           "_heatmap.png"),
    res=300,unit="cm",height=30,width=30)

pheatmap(keywords_heatmap_wide,
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","gray","skyblue","gold", "palevioletred3"))(30))

dev.off() 

# elements of the plot

keywords <- trends_categories |>
    filter(keyword %in% unique(c(keywords_heatmap_long$from,
                                 keywords_heatmap_long$to))) |>
    mutate(from=factor(keyword,
                       levels=as.character(unique(keyword)))) |>
    mutate(to=factor(keyword,
                     levels=as.character(unique(unique(keyword)))),
           count=0,
           count_bin="0",
           jaccard=0) |>
    dplyr::select(from,to,count,count_bin)

diagonal <- tibble(from=factor(keywords$from,
                               levels = unique(keywords_heatmap_long$from[order(keywords_heatmap_long$category.x,keywords_heatmap_long$from)])),
                   to=factor(keywords$to,
                                   levels = unique(keywords_heatmap_long$to[order(keywords_heatmap_long$category.y,
                                                                                  keywords_heatmap_long$to)])),
                   count=-1,
                   jaccard=0) 

# we defined here the diagonal because the raw values don't include them.
# In addition we need the diagonal seperate from the raw data because we will paint it differently

## summaries to dynamically set the break points and limits of the plot

summary <- summary(keywords_heatmap_long$count)

## For the heatmap we need breaks to define the specific points of the legend
## and limits to ensure that all values will be included in the plot.
## To create breaks and limits and make them scalable we used the base R functions summary and quantile.
## also quantiles because the raw counts are far apart

quantile <- as.vector(quantile(keywords_heatmap_long$count,probs = c(50,90, 95,98)/100)) # big probabilities because of order of magnitude difference of values 

breaks <- c(floor(min(summary)),
            round(quantile[1]),
            round(quantile[2]),
            round(quantile[3]),
            round(quantile[4]),
            ceiling(max(summary))) 

limits=c(min(breaks),max(breaks))

# add the order based on the categories so they appear in that order
keywords_heatmap_long$from <- factor(keywords_heatmap_long$from, 
                                     levels = unique(keywords_heatmap_long$from[order(keywords_heatmap_long$category.x,
                                                                                      keywords_heatmap_long$from)]))
keywords_heatmap_long$to <- factor(keywords_heatmap_long$to, 
                                   levels = unique(keywords_heatmap_long$to[order(keywords_heatmap_long$category.y,
                                                                                  keywords_heatmap_long$to)]))


# running the plot

g <- guide_legend("no of abstracts") # legend title and combination of different colors and shapes into one legend
keywords_heatmap_long <- keywords_heatmap_long %>% filter(!is.na(count_bin))
pubmed_keyword_coocurrence_heatmap <- ggplot()+
  geom_tile(data=keywords_heatmap_long,
            aes(x=from,
                y=to,
                fill=count_bin),
            alpha=1,
            show.legend = T)+
  geom_tile(data=keywords_heatmap_long,
            aes(x=from, y=from),
            fill="white",
            alpha=1,
            show.legend = F)+
  geom_point(data=keywords_heatmap_long,aes(x=from, y=from, color="gray50"),alpha=1, show.legend = F)+
#  scale_fill_manual(values=c("#d73027","#fc8d59","#fee090","#4575b4","#91bfdb","#e0f3f8","white")) + #,breaks=breaks,limits=limits) +
  scale_color_manual(values=c("gray80"))+
  scale_x_discrete(position = "top")+
  scale_y_discrete(limits = rev(levels(keywords_heatmap_long$to)))+
  guides(fill = guide_legend("# of abstracts"), color=FALSE)+
  xlab("") +
  ylab("")+
  theme_bw()+
  theme(plot.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        text = element_text(size=17), 
        axis.text.x = element_text(angle = 90, hjust = 0),
        legend.position = c(.90, .83))

  ggsave(paste0("../plots/",
              user_prefix,"_",
              format(Sys.time(),
                     "%Y-%m-%d_%H-%M"),
              "_heatmap2.png"),
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
            aes(x=from,
                y=to,
                fill=log2),
            alpha=0,
            show.legend = F)+
  geom_point(data=keywords_heatmap_long,
             aes(x=from,
                 y=to,
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
             aes(x=from,
                 y=to),
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
    left_join(keywords_n_papers,c("from"="keyword")) |> 
    left_join(keywords_n_papers,c("to"="keyword")) |> 
    mutate(jaccard_index=count/(n_papers.x+n_papers.y-count),
           random_expectation=(count/n_papers_pubmed)/(freq.x*freq.y))

colors_j <- c("steelblue1","yellow","goldenrod1","orange")
limits<- c(0,
           round(max(keywords_heatmap_jaccard$jaccard_index[keywords_heatmap_jaccard$jaccard_index < 1])*1.07,4) )
# the limits are the maximum value multiplied by >1 so it is bigger than the maximum value

pubmed_jaccard_heatmap <- ggplot()+ 
    geom_tile(data=keywords_heatmap_jaccard,
              aes(x=from,
                  y=to,
                  fill=jaccard_index),
              alpha=1,
              show.legend = T,
              colour = "white")+
    scale_fill_gradientn(colours =colors_j,limits=limits)+
    geom_tile(data=diagonal,
              aes(x=from,
                  y=to,
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
              aes(x=from,
                  y=to,
                  fill=random_expectation),
              alpha=1,
              show.legend = T,
              colour = "white")+
    geom_tile(data=diagonal,
              aes(x=from,
                  y=to,
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
