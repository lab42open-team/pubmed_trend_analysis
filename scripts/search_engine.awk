#! /usr/bin/gawk -f
#  gunzip -c dataset/*.tsv.gz | ./scripts/search_engine.awk keywords.txt - > results.tsv
#
#
BEGIN {
    FS="\t"
    }

# Treat keywords
(ARGIND==1){
    $key=tolower($1)
    keywords[$key]++
}

# Treat corpus
(ARGIND==2){

#    if ( $4 ~ /[[:alnum:]]/ && $5 ~ /[[:alnum:]]/ && $6 ~ /[[:alnum:]]/) {

        line=tolower($5$6)

        sub("\\|.*$","",$1)
        
        for (k in keywords){
            regex=k #"[:blank:]"
            if (match(line, regex)){
                corpus[$1]= $4 "\t" line "\t" k
            }
        }
#    }
}

END{ 
#    print length(corpus)
#    print length(keywords)
    for (i in corpus){
            print i "\t" corpus[i]
    }
}

