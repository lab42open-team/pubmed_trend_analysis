#! /usr/bin/gawk -f
# Author Savvas Paragkamian
#
# Run :  gunzip -c dataset/*.tsv.gz | ./scripts/search_engine.awk keywords.txt - > results.tsv
# Note the - , it takes the input from the pipe |
#
BEGIN {
    FS="\t"
    }

# Treat keywords
(ARGIND==1){
    keywords[$1]=tolower($1)
}

# Treat corpus. The input comes from the - from the execution of the script
(ARGIND==2){

    if ( $4 ~ /[[:alnum:]]/ && $5 ~ /[[:alnum:]]/ && $6 ~ /[[:alnum:]]/) {

        line=tolower($5$6)

        sub("\\|.*$","",$1)
        
        for (k in keywords){
            #regex=k "$"  #"[:blank:]"
            if (match(line,"[^[:alpha:]]" keywords[k] "[^[:alpha:]]")){
                print $1 "\t" $4 "\t" line "\t" keywords[k]
            }
        }
    }
}

