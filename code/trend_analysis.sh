#!/bin/bash
# this script is about finding matches of keywords in the pubmed tsv files

########################### user arguments#####################################
# There are 3 user arguments. -k is for the keywords, -d is for pubmed destination and -p is for the user defined prefix of all files created from the analysis (e.g text files and plots). -d has a "default" option to use the inhouse pubmed downloads.

while [ "$1" != "" ]; do
  case $1 in
    -k | --keywords )
      shift # after the match of the -k | --keywords it moves the counter to the next number to take the import. Otherwise the value of $1 would have been "-k" and not the user argument
      user_keywords=$1
      ;;

    -d | --data )
      shift
      if [ "$1" == "default" ];then # the default option of Pubmed location
        pubmed_path="/data/databases/pubmed/"
      else
        pubmed_path=$1
      fi
      ;;
    -p | --prefix)
      shift
      user_prefix=$1
      ;;

  esac
  shift

done

echo 'The arguments passed are: '
echo 'keywords file: ' $user_keywords
echo 'PubMed folder: '$pubmed_path
echo 'Prefix: ' $user_prefix

###############################################################################

################################## Initiation #################################
time_start=`date +%s`
DATE=$(date +"%Y-%m-%d_%H-%M")

output="../data/${user_prefix}_${DATE}_trends_pubmed.tsv"
touch $output

files=$(find $pubmed_path -name \*.tsv.gz) # the files that we will search for patterns in the directory that they are stored
files_number=`find $pubmed_path -name \*.tsv.gz| wc -l`

mapfile -t keywords_file < $user_keywords # load keywords to array
keywords_number=${#keywords_file[@]}

total_repeats=$(($keywords_number*$files_number))

i=0

echo "Processing $files_number files with $keywords_number keywords"
echo "Percentage% =" `expr $((($i/$total_repeats)*100))`

############################## GREP SEARCH #####################################

while IFS= read -r keyword; do # there is a text file containing the keywords, each line has a keyword
  

      regex_keyword="$(echo $keyword | gawk '{gsub(/ /,"[- ]",$0)}{key="\\b"$0"[[:alpha:]]?\\b"; print key}')"
      
      echo $keyword
      echo $regex_keyword
  
  for file in $files; do
      
      
      zgrep -inE "$regex_keyword" $file | awk -v var="$keyword" -v file="$file" 'BEGIN{FS="\t"; OFS="\t"} {sub(":","\t"); print $1,$2,var,file,$5}' >> $output # grep in zipped files, -i for case insensitive,-w for searching the whole pattern and -n to return the line of the match. 
      ((i+=1))
      
      done

  
  awk -v i=$i -v total_repeats=$total_repeats 'BEGIN { print "Percentage = " (i/total_repeats) }'

done < $user_keywords
############################## END OF SEARCH ####################################

time_end=`date +%s`
time_exec=`expr $(( $time_end - $time_start ))`

echo "PubMed Trends file created! Execution time was $time_exec seconds"
echo "Plotting results..."

########################### CALLING R SCRIPT #####################################
Rscript "trends_plots.r" "$output" "$user_prefix"

echo "Finished!"
