#!/bin/sh
# this script is about finding matches of keywords in the pubmed tsv files

DATE=$(date '+%Y-%m-%d') 
touch ${DATE}_trends_pubmed.tsv

files=$(find /data/databases/pubmed/ -name \*.tsv.gz) # the files that we will search for patterns in the directory that they are stored

while read keyword; do # there is a text file containing the keywords, each line has a keyword
  for file in $files
  do
    if zgrep -qinw "$keyword" $file; then # here the program checks if there is a match of the keyword in the specific file, if there is a match then it continuous to print the output
   
      zgrep -inw "$keyword" $file | awk -v var="$keyword" -v file="$file" 'BEGIN{FS="\t"; OFS="\t"} {sub(":","\t"); print $1,$2,var,file,$5}' >> ${DATE}_trends_pubmed.tsv # grep in zipped files, -i for case insensitive,-w for searching the whole pattern and -n to return the line of the match. 
    
    else
      echo "no matches $keyword in $file"
  fi
  
  done
done < keywords.txt

echo "done!!"
