#!/bin/bash

while [ "$1" != "" ]; do
  case $1 in
    -k | --keywords )
      shift
      user_keywords=$1
      ;;

    -d | --data )
      shift
      if [ "$1" == "default" ];then
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

