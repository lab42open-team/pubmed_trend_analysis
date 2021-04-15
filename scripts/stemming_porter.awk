#! /usr/bin/gawk -f

#  This is the Porter stemming algorithm, coded up in awk by Gregory Grefenstette
#  July 5, 2012
#  It follows the  algorithm presented in
#
#  Porter, 1980, An algorithm for suffix stripping, Program, Vol. 14,
#  no. 3, pp 130-137,
#
# and more precisely the code for the ANSI C version found at
#
#  http://www.tartarus.org/~martin/PorterStemmer
#
# This endioding of the algorthm can be used free of charge for any purpose




# TRUE if last two characters are a double consonant

function doublec(s)
{ if(substr(s,length(s),1) == substr(s,length(s)-1,1) && (substr(s,length(s),1)  !~ /[aeiou]/)) return 1;
                                                         else return 0 }

# antyhing other than one of a, e, i, o, u, for the case of "y" it checks if
# there is a preceding vowel, in which case it is considered a consonant

function cons(str,i)
{ if(substr(str,i,1) ~ /[aeiou]/) return 0 ;
  if(i==1) return 1;
  if(substr(str,i,1) == "y") { if(substr(str,i-1,1)  ~ /[aeiou]/) return 1; else return 0 }
  return 1
}

# cvc(i) is TRUE if last three characters of str are  consonant - vowel - consonant
# and also if the second c is not w,x or y. this is used when trying to
#   restore an e at the end of a short word. e.g.
#
#  cav(e), lov(e), hop(e), crim(e), but
#  snow, box, tray.

function cvc(str) {
  if(length(str) <= 2) return 0;
  if( str ~ /[wxy]$/ ) return 0;
  if ( cons(str,length(str)-2) && !cons(str,length(str)-1) &&  cons(str,length(str)) ) return 1;
  return 0
}


# m() measures the number of consonant sequences between k0 and j. if c is
# a consonant sequence and v a vowel sequence, and <..> indicates optional
#
#      <c><v>       gives 0
#      <c>vc<v>     gives 1
#      <c>vcvc<v>   gives 2
#      <c>vcvcvc<v> gives 3
#      ....
# this version returns "2" as the maximum value

function m(str) {
  # skip initial consonants
  mreturns=0;
  mindex=1;
  while((mindex <= length-str) && cons(str,mindex)) mindex++ ;
  while (1)
    { while(1)
        { if (mindex > length(str)) return mreturns;
          if( cons(str,mindex) ) break;
          mindex++
        }
      mindex++;
      mreturns++;
      if(mreturns > 2) return mreturns;
      while(1)
        { if (mindex > length(str)) return mreturns;
          if( ! cons(str,mindex) ) break;
          mindex++
        }
    }
}




# step1ab() gets rid of plurals and -ed or -ing. e.g.
#
#      caresses  ->  caress
#      ponies    ->  poni
#      ties      ->  ti
#      caress    ->  caress
#      cats      ->  cat
#
#      feed      ->  feed
#      agreed    ->  agree
#      disabled  ->  disable
#
#      matting   ->  mat
#      mating    ->  mate
#      meeting   ->  meet
#      milling   ->  mill
#      messing   ->  mess
#
#      meetings  ->  meet


function step1ab(str) {

if(str ~ /sses$/ || str ~ /ies$/ ) str=substr(str,1,length(str)-2) ;
else if(str ~ /ss$/) ;
else if (str ~ /s$/) str=substr(str,1,length(str)-1) ;

if(str ~ /eed$/) { if(m(substr(str,1,length(str)-3))>0 ) str=substr(str,1,length(str)-1) ; }
else {trunc=0;
    if (str ~ /[aeiouy].*ed$/) trunc=2;
    if (str ~ /[aeiouy].*ing$/) trunc=3 ;
    if(trunc>0) { str=substr(str,1,length(str)-trunc) ;
    if(str ~ /(at|bl|iz)$/) str=str"e" ;
    else
     if (doublec(str)==1 && (str !~ /[lsz]$/)) {  str=substr(str,1,length(str)-1); }
     else if( m(str)==1 && cvc(str)) str=str"e";
    }}
    return str }

# step1c() turns terminal y to i when there is another vowel in the stem.

function step1c(str) {
if(str ~/[aeiouy].*y$/) str=substr(str,1,length(str)-1)"i" ;
 return str }

# step2() maps double suffices to single ones. so -ization ( = -ize plus
#   -ation) maps to -ize etc. note that the string before the suffix must give
#           m() > 0. */

function step2(str) {
  if( str ~ /[aeiouy][^aeiouy].*ational$/ ) str=substr(str,1,length(str)-5)"e" ;
  else if ( str ~ /[aeiouy][^aeiou].*tional$/ ) str=substr(str,1,length(str)-2) ;
  else if ( str ~ /[aeiouy][^aeiou].*[ae]nci$/ ) str=substr(str,1,length(str)-1)"e" ;
  else if ( str ~ /[aeiouy][^aeiou].*izer$/ ) str=substr(str,1,length(str)-1) ;
  else if ( str ~ /[aeiouy][^aeiou].*bli$/ ) str=substr(str,1,length(str)-1)"e" ;
  else if ( str ~ /[aeiouy][^aeiou].*alli$/ ) str=substr(str,1,length(str)-2);
  else if ( str ~ /[aeiouy][^aeiou].*entli$/ ) str=substr(str,1,length(str)-2);
  else if ( str ~ /[aeiouy][^aeiou].*eli$/ ) str=substr(str,1,length(str)-2);
  else if ( str ~ /[aeiouy][^aeiou].*ousli$/ ) str=substr(str,1,length(str)-2);
  else if ( str ~ /[aeiouy][^aeiou].*ization$/ ) str=substr(str,1,length(str)-5)"e" ;
  else if ( str ~ /[aeiouy][^aeiou].*ation$/ ) str=substr(str,1,length(str)-3)"e" ;
  else if ( str ~ /[aeiouy][^aeiou].*ator$/ ) str=substr(str,1,length(str)-2)"e" ;
  else if ( str ~ /[aeiouy][^aeiou].*alism$/ ) str=substr(str,1,length(str)-3);
  else if ( str ~ /[aeiouy][^aeiou].*iveness$/ ) str=substr(str,1,length(str)-4) ;
  else if ( str ~ /[aeiouy][^aeiou].*fulness$/ ) str=substr(str,1,length(str)-4) ;
  else if ( str ~ /[aeiouy][^aeiou].*ousness$/ ) str=substr(str,1,length(str)-4) ;
  else if ( str ~ /[aeiouy][^aeiou].*aliti$/ ) str=substr(str,1,length(str)-3) ;
  else if ( str ~ /[aeiouy][^aeiou].*iviti$/ ) str=substr(str,1,length(str)-3)"e" ;
  else if ( str ~ /[aeiouy][^aeiou].*biliti$/ ) str=substr(str,1,length(str)-5)"le" ;
  else if ( str ~ /[aeiouy][^aeiou].*logi$/ ) str=substr(str,1,length(str)-1) ;
 return str }

# step3() deals with -ic-, -full, -ness etc. similar strategy to step2.
function step3(str) {
  if( str ~ /[aeiouy][^aeiouy].*icate$/ ) str=substr(str,1,length(str)-3) ;
  else if ( str ~ /[aeiouy][^aeiou].*ative$/ ) str=substr(str,1,length(str)-5);
  else if ( str ~ /[aeiouy][^aeiou].*alize$/ ) str=substr(str,1,length(str)-3) ;
  else if ( str ~ /[aeiouy][^aeiou].*iciti$/ ) str=substr(str,1,length(str)-3) ;
  else if ( str ~ /[aeiouy][^aeiou].*ical$/ ) str=substr(str,1,length(str)-2) ;
  else if ( str ~ /[aeiouy][^aeiou].*ful$/ ) str=substr(str,1,length(str)-3);
  else if ( str ~ /[aeiouy][^aeiou].*ness$/ ) str=substr(str,1,length(str)-4);
return str }


# step4() takes off -ant, -ence etc., in context <c>vcvc<v>.
function step4(str) {
  if( str ~ /al$/ ) { if ( m(substr(str,1,length(str)-2)) > 1 )  str=substr(str,1,length(str)-2) }
  else if ( str ~ /[ae]nce$/ ) { if ( m(substr(str,1,length(str)-4)) > 1) str=substr(str,1,length(str)-4) }
  else if ( str ~ /(er|ic)$/ ) { if (  m(substr(str,1,length(str)-2)) > 1 ) str=substr(str,1,length(str)-2) }
  else if ( str ~ /[ai]ble$/ ) { if (  m(substr(str,1,length(str)-4)) > 1 ) str=substr(str,1,length(str)-4) }
  else if ( str ~ /ant$/ ) { if (  m(substr(str,1,length(str)-3)) > 1 )  str=substr(str,1,length(str)-3) }
  else if ( str ~ /ement$/) { if(  m(substr(str,1,length(str)-5)) > 1 ) str=substr(str,1,length(str)-5) }
  else if ( str ~ /ment$/) { if (  m(substr(str,1,length(str)-4)) > 1 ) str=substr(str,1,length(str)-4) }
  else if ( str ~ /ent$/) { if (  m(substr(str,1,length(str)-3)) > 1 ) str=substr(str,1,length(str)-3) }
  else if ( str ~ /[st]ion$/) { if (  m(substr(str,1,length(str)-3)) > 1 )  str=substr(str,1,length(str)-3) }
  else if ( str ~ /ou$/) { if (  m(substr(str,1,length(str)-2)) > 1 )  str=substr(str,1,length(str)-2) }
  else if ( str ~ /(ism|ate|iti|ous|ive|ize)$/) { if (  m(substr(str,1,length(str)-3)) > 1 )  str=substr(str,1,length(str)-3)
}
return str}

# step5() removes a final -e if m() > 1, and changes -ll to -l if
#  m() > 1.
function step5(str) {
  if ( str ~ /e$/ && ( m(str)>1 || (m(str)==1  && !cvc(substr(str,1,length(str)-1)))))  str=substr(str,1,length(str)-1) ;
  if( str ~ /ll$/ && m(str)>1 ) str=substr(str,1,length(str)-1) ;
 return str }


function stem(str)
{ str=tolower(str);
  if(length(str)<=2) return str;
  str=step1ab(str);
  str=step1c(str);
  str=step2(str);
  str=step3(str);
  str=step4(str);
  str=step5(str);
return str
}


# main
{  printf("%s",stem($1));
   for(i=2;i<=NF;i++) printf("%s%s",FS,stem($i));
   print ""
}
