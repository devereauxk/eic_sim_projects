#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

WORKING_DIR=./eAu_10_100_qhat0_nlo

if [ -z "$1" ]
then
        echo "No job number set."
        echo "Please run as ./run_kdebug.sh jobnumber"
        echo "Exiting..."
        exit 1
fi

INPUT=$1

echo "start running in directory $PWD"

echo "Running Job Number $1"
while read line; do
  # if line of file starts with "@kdebug " then prints whole line to $fout
  first_word=`echo $line | awk '{print $1;}'`
  if [[ $first_word == "@kdebug" ]]
  then
    echo $line >> kdebug_bin.txt
  fi
done < $WORKING_DIR/logs/eAu_${INPUT}.log
echo "Completed!!!"

echo ""

#Move output files and cleanup
echo "Cleaning Up..."
mv -v kdebug_bin.txt ${WORKING_DIR}/logs/kdebug_bin_${INPUT}.txt
echo "DONE!!!"
