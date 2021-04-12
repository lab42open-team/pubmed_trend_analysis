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

    $line=tolower($0)
    sub("\\|.*$","",$1)

    for (k in keywords){
        if ($line ~ k){
            corpus[$1]= $4 "\t" k
        }
    }
}

END{ 
#    print length(corpus)
#    print length(keywords)
    for (i in corpus){
            print i "\t" corpus[i]
    }
}

