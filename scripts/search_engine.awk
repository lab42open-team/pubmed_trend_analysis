#! /usr/bin/gawk -f
BEGIN {
    FS="\t"
    }

# Treat corpus
(ARGIND==1){

    corpus[$1]=$4 "\t" tolower($5) "\t" tolower($6)
}

# Treat keywords
(ARGIND==2){
    
    $key=tolower($1)
    keywords[$key]++

}

END{ 
    print length(corpus)
    print length(keywords)
    
    for (i in corpus){
        for (k in keywords){
            if (corpus[i] ~ k){
            print i "\t" corpus[i] "\t" k
            }

        }
    }
}

