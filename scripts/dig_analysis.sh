#!/bin/bash
# this script is about finding matches of keywords in the pubmed tsv files
#
# Author:   Savvas Paragkamian (s.paragkamian@hcmr.gr)
#           Institute of Marine Biology Biotechnology and Aquaculture (IMBBC)
#           Hellenic Centre for Marine Research (HCMR)
#
# Created:  05/06/2020
# License:  GNU GPLv3 license
# Developer : Savvas Paragkamian
#
#
########################### user arguments#####################################
# There are 3 user arguments. -k is for the keywords, -d is for pubmed destination and -p is for the user defined prefix of all files created from the analysis (e.g text files and plots). -d has a "default" option to use the inhouse pubmed downloads.

## Usage of the script
usage="Use the parameter -k for the keywords file, -d for the corpus directory (expects the path to the PubMed data) and -p a string with prefix of all the generated files (txt and plots). \nExample: ./scripts/dig_analysis.sh -k keywords.txt -d default -p "ecology_trends" \n"

## Set the working directory

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # determine the path of the script
cd $DIR # set the scripts directory as working directory

## User input parameters of keywords, corpus directore and prefix of the analysis.
while getopts "k:d:p:" option
do
   case "$option" in
        k)   user_keywords="../${OPTARG}";;
        d)   pubmed_path="${OPTARG}"
            if [ "$pubmed_path" == "default" ];then # the default option of Pubmed location
                pubmed_path="/data/databases/pubmed/"
            else
                pubmed_path="${OPTARG}"
            fi;;
        p) user_prefix="${OPTARG}";;
        ?|:)   echo -e "$usage" ; exit 1;;
        *)   echo -e "option ${OPTARG} unknown. \n$usage" ; exit 1 ;;
   esac
done

## Detect if no options were passed
if [ $OPTIND -eq 1 ]; 
    then echo -e "No options were passed. \n$usage "; exit 1; 
fi

## Detect if a single option is missing
if [ -z "${user_keywords}" ]; then
    echo -e "Option -k empty. $usage"; exit 1;
fi

if [ -z "${pubmed_path}" ]; then
    echo -e "Option -d empty. $usage"; exit 1;
fi
if [ -z "${user_prefix}" ]; then
    echo -e "Option -p empty. $usage"; exit 1;
fi
# Initiation of Trend Analysis
echo 'The arguments passed are: '
echo 'keywords file: ' $user_keywords
echo 'PubMed folder: '$pubmed_path
echo 'Prefix: ' $user_prefix

###############################################################################

################################## Initiation #################################

time_start=`date +%s`
DATE=$(date +"%Y-%m-%d_%H-%M")
output="../data/${user_prefix}_${DATE}_dig_analysis.tsv"
touch $output

echo 'Data will be stored in ' $output

files=$(find $pubmed_path -name \*.tsv.gz) # the files that we will search for patterns in the directory that they are stored
files_number=`find $pubmed_path -name \*.tsv.gz| wc -l`

mapfile -t keywords_file < $user_keywords # load keywords to array
keywords_number=${#keywords_file[@]}

total_repeats=$(($keywords_number*$files_number))

#echo "Processing $files_number files with $keywords_number keywords"
echo "Percentage% =" `expr $((($i/$total_repeats)*100))`

############################## PATTERN SEARCH #####################################

gunzip -c $files | ./search_engine.awk $user_keywords - > $output

################################################################################

#i=0
#while IFS= read -r keyword; do # there is a text file containing the keywords, each line has a keyword
#  
#
#      regex_keyword="$(echo $keyword | gawk '{gsub(/ /,"[- ]",$0)}{key="\\b"$0"[[:alpha:]]?\\b"; print key}')"
#      
#      echo $keyword
#      echo $regex_keyword
#  
#  for file in $files; do
#      
#      
#      zgrep -inE "$regex_keyword" $file | awk -v var="$keyword" -v file="$file" 'BEGIN{FS="\t"; OFS="\t"} {sub(":","\t"); print $1,$2,var,file,$5}' >> $output # grep in zipped files, -i for case insensitive,-w for searching the whole pattern and -n to return the line of the match. 
#      ((i+=1))
#      
#      done
#
#  
#  awk -v i=$i -v total_repeats=$total_repeats 'BEGIN { print "Percentage = " (i/total_repeats) }'
#
#done < $user_keywords
############################### END OF SEARCH ####################################

time_end=`date +%s`
time_exec=`expr $(( $time_end - $time_start ))`

echo "D|G analysis file created! Execution time was $time_exec seconds"
echo "Plotting results..."

########################### CALLING R SCRIPT #####################################
Rscript "trends_plots.r" "$output" "$user_prefix"

echo "Finished!"
