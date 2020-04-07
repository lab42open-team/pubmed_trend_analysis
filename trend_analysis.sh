#!/bin/bash
# this script is about finding matches of keywords in the pubmed tsv files

user_keywords=$1
pubmed_path=$2

DATE=$(date '+%Y-%m-%d')
time_start=`date +%s`

touch ${DATE}_trends_pubmed.tsv

files=$(find $pubmed_path -name \*.tsv.gz) # the files that we will search for patterns in the directory that they are stored
files_number=`find $pubmed_path -name \*.tsv.gz| wc -l`

mapfile -t keywords_file < $user_keywords # load keywords to array
keywords_number=${#keywords_file[@]}

total_repeats=$(($keywords_number*$files_number))

i=0

echo "Processing $files_number files with $keywords_number keywords"
echo "Percentage% =" `expr $((($i/$total_repeats)*100))`

while read keyword; do # there is a text file containing the keywords, each line has a keyword
  for file in $files
  do
    #if zgrep -qinw "$keyword" $file; then # here the program checks if there is a match of the keyword in the specific file, if there is a match then it continuous to print the output
   
      zgrep -inw "$keyword" $file | awk -v var="$keyword" -v file="$file" 'BEGIN{FS="\t"; OFS="\t"} {sub(":","\t"); print $1,$2,var,file,$5}' >> ${DATE}_trends_pubmed.tsv # grep in zipped files, -i for case insensitive,-w for searching the whole pattern and -n to return the line of the match. 
      ((i+=1))
      
    #else
      #echo "no matches $keyword in $file"  
  done
  
  echo $keyword
  echo "Percentage% =" `expr $((($i/$total_repeats)*100))`
done < keywords.txt

time_end=`date +%s`
time_exec=`expr $(( $time_end - $time_start ))`

echo "PubMed Trends file created! Execution time was $time_exec seconds"
echo "Plotting results..."

Rscript "trends_plots.r"

echo "Finished!"
