#!/bin/bash
j=1
l=0
input_file=$1
output_file=$2
num_lines=0
#Poison reset
poison=0
#Lines reset due to possible error message
num_lines=0
#Simple reset
Simple=0
#Geoms reset
Geoms=0
#CSV reset
CSV=0
#Check for folders
Folders=`ls -d */`
#Functions
error_popup () {
  echo "ERROR: $1 file $2 does not exist." ; echo "Write \". TOP.bat -h\" for help."
}
help_us_obi1 () {
  echo "======================================================================  HELP FILE  ======================================================================"
  echo "This script is used to transform results of AUTOOPT to a more human-readable format. Made by Alexandr Zaykov."
  echo "It needs two files to run for the first time: simple input file and AUTOOPT output file named OptResultsReduced"
  echo "To run the script write \". TOP.bat simple_input_file AUTOOPT_output_file\". Optionally add a number of lines (molecular pairs) you want to investigate."
  echo "The script will then create separate folder for each structure and copy and modify the input as \"Simple.inp\"."
  echo "It will then ask you if you want to run Simple in these folders. If so, it will create a logfile named \"logSimple\"."
  echo "Also, the script can extract geometries into molden readable file - \"Geometries.xyz\". Each structure will have one frame."
  echo "The last thing this script can do is to extract the data from Simple calculation into .csv formatted file \"Table_Data.csv\" readable by table processors."
  echo
  echo "TL; DR:"
  echo "INPUT: Original Simple input file, AUTOOPT output named \"OptResultsReduced\""
  echo "START: \$ . TOP.bat simple_input_file AUTOOPT_output_file"
  echo "OUTPUT: Table_Data.csv, Geometries.xyz, Structure_xy.xyz, Simple log files in each of the folders"
  echo "========================================================================================================================================================="
}
poison () {
  poison=1
}
num_pairs () {
  echo "The script will use $1 lines. $2"
  num_lines=$1
}
yn_switch () {
  echo
  read -n1 -p "$1" doit
  case $doit in
      y|Y) printf -v "$2" '%s' '1' ;;  ## Assigns a variable with variable name, avoids eval and other evils
      n|N) printf -v "$2" '%s' '0' ;;
      *) ;;
  esac
}
run_simple () {
  [[ -d $1 ]] && cd $1 &> /dev/null && let Counter=$2*100 && let Percent=Counter/$3 && echo -ne "$Percent %\r" && Simple Simple.inp > logSimple ; cd - &> /dev/null
}
check_file () {
  [[ $1 == Structure* ]] && { cd $j &> /dev/null ; [ ! -f $2 ] && printf -v "$3" '%s' '0' ; cd - &> /dev/null ; }
}
show_percent () {
  let Counter=$1*100
  let Percent=Counter/$2
  echo -ne "$Percent %\r"
}
#Check if the files even exist, if it does not prints a message and performs harakiri, also check for help request

[ "$input_file" = "-h" ] && help_us_obi1 && poison || { [ ! -f "$input_file" ] && error_popup "Simple input file" $input_file && poison ; [ ! -f "$output_file" ] && error_popup "AUTOOPT output file" $output_file && poison ; }

#Check optional lines input and precheck if harakiri is needed, ends at the bottom of the file - very long if, I am sorry.

if [ ! $poison -eq 1 ] ; then

#Read number of pairs chosen by the user
  [[ -z ${3} ]] && num_pairs 20 "(default)" || num_pairs $3 "(user set)"

  #Run the body of the liner
  k=1
  while [[ $k -le $num_lines ]] ; do
    j=$(($k + 1))
    if [ -d "Structure$k" ] ; then
      echo -ne "Folder(s) with the name Structure already exists. Continuing...\r"
    else
      echo -ne "Last made dir is Structure$k!                                      \r"
      mkdir Structure$k
      cd Structure$k &> /dev/null
        cp ../$input_file Simple.inp
        sed -i '/Const :/d' Simple.inp
        sed -i '/T */d' Simple.inp
        sed -i '/R */d' Simple.inp
        echo 'Const :' >> Simple.inp 
        awk -v var=$j 'NR==var {print "T Z " $3  >> "Simple.inp"}' ../$output_file
        awk -v var=$j 'NR==var {print "T Y " $4  >> "Simple.inp"}' ../$output_file
        awk -v var=$j 'NR==var {print "T X " $5 >> "Simple.inp"}' ../$output_file
        awk -v var=$j 'NR==var {print "R Z " $6 >> "Simple.inp"}' ../$output_file
        awk -v var=$j 'NR==var {print "R Y " $7 >> "Simple.inp"}' ../$output_file
        awk -v var=$j 'NR==var {print "R X " $8 >> "Simple.inp"}' ../$output_file
        echo 'dE(CT) : 1.0' >> Simple.inp
        echo 'Reduce : 1.0' >> Simple.inp
      cd - &> /dev/null
    fi
    k=$(($k + 1))
  done
  yn_switch "Do you want to run Simple now? [y,n]" Simple
  yn_switch "Do you want geometries to be extracted into a separate file? [y,n]" Geoms
  yn_switch "Do you want the results to be extracted into a .CSV file? [y,n]" CSV
  echo

#Check if Simple input files exist
ifiles=1
[ $Simple -eq 1 ] && for j in $Folders ; do
                       check_file $j "Simple.inp" ifiles
                     done  
#Simpler part of the script, pun inteded on all levels
if [ $ifiles -eq 1  ] ; then
  if [ $Simple -eq 1 ] ; then
    x=1
    Count=0

#Counter
    for j in $Folders ; do
      [[ $j == Structure* ]] && Count=$(($Count + 1)) 
    done 
    k=1
    while [[ $k -le $num_lines ]] ; do 
      run_simple Structure$k k Count
      k=$(($k + 1))
    done
    echo Simpled!
  fi
#Check if Simple input files exist, but first check if it is even necessary to check
ofiles=1
[ $Geoms -eq 1 ] || [ $CSV -eq 1 ] && for j in $Folders ; do
                                        check_file $j "logSimple" ofiles
                                      done 
#Geometry part of the script
if [ $ofiles -eq 1 ] ; then
  if [ $Geoms -eq 1 ] ; then
#resets
    num_atoms=0
    twice_atoms=0
    atoms=1
    x=1
#Geometries part
    rm -r Geoms &> /dev/null
    mkdir Geoms
    touch Geoms/Geometries.xyz
#Start gathering structures
    cd Structure1 &> /dev/null
      num_atoms=$(sed -n '/Geometry A/,/Geometry B/p' logSimple | wc -l)  #number of atoms counter
      num_atoms=$(($num_atoms - 2))
      twice_atoms=$(($num_atoms + $num_atoms))
      echo The molecule has $num_atoms atoms.
    cd - &> /dev/null
#body
    k=1
    while [[ $k -le $num_lines ]] ; do
      if [[ -d "Structure$k" ]] ; then        
        j=Structure$k
        show_percent k Count
        cd $j &> /dev/null
#Find string and print it to the file
#Geometries
          echo $twice_atoms >> ../Geoms/Geometries.xyz
          echo $k: >> ../Geoms/Geometries.xyz
          awk '/Geometry A:/ {for(i=1; i<='$num_atoms'; i++) {getline; print}}' logSimple >> ../Geoms/Geometries.xyz
          awk '/Transformed geometry of B/ {for(i=1; i<='$num_atoms'; i++) {getline; print}}' logSimple >> ../Geoms/Geometries.xyz
          touch ../Geoms/Structure_$k.xyz
          echo $twice_atoms >> ../Geoms/Structure_$k.xyz
          echo Structure $k: >> ../Geoms/Structure_$k.xyz
          awk '/Geometry A:/ {for(i=1; i<='$num_atoms'; i++) {getline; print}}' logSimple >> ../Geoms/Structure_$k.xyz
          awk '/Transformed geometry of B/ {for(i=1; i<='$num_atoms'; i++) {getline; print}}' logSimple >> ../Geoms/Structure_$k.xyz
        cd - &> /dev/null
#Necessary deletions
        while [ $atoms -lt 120 ]; do
          sed -i "s/   $atoms //g" ./Geoms/Structure_$k.xyz
          sed -i "s/  $atoms //g" ./Geoms/Structure_$k.xyz
          let atoms=atoms+1
        done
        atoms=1
      fi
      k=$(($k + 1))
    done
#deletes and mods
    echo -ne "Deleting unnecessary characters\r"
    atoms=1
    while [ $atoms -le 120 ]; do
      sed -i "s/   $atoms //g" ./Geoms/Geometries.xyz
      sed -i "s/  $atoms //g" ./Geoms/Geometries.xyz
      let atoms=atoms+1
    done
    atoms=1
    echo Geometries written in Geoms folder!  
  fi
#Geoms ending, CSV starting
  if [ $CSV -eq 1 ] ; then
    rm Table_Data.csv &> /dev/null
    touch Table_Data.csv
#Headers
    echo "Structure, T A (pert.), T B(pert.), T (S+), T (S-), k(Marcus eq.), 2Vab, dE (TT - S+), dE (TT - S-), dE (Davydov spl.), dE (Biexciton binding energy), dE (Overall endoergicity)" >> Table_Data.csv
#Body
    k=1
    while [[ $k -le $num_lines ]] ; do
      if [[ -d "Structure$k" ]] ; then        
        show_percent k Count
        cd Structure$k &> /dev/null
          line=$(awk '/CSV:/ {getline; print}' logSimple)    #awk a line
          echo "$k, $line" >> ../Table_Data.csv
        cd - &> /dev/null
        k=$(($k + 1))
      fi
    done
    echo Results written in Table_Data.csv!
fi


#FATAL errors
else
  error_popup "Simple output" "logSimple"
  echo SOLUTION: Run Simple first within this script.
fi
else
  error_popup "Simple input" "Simple.inp"
  echo SOLUTION: The shortest route to fix is to delete the Structure folders and start this script again.
fi
#END OF POISON
fi

