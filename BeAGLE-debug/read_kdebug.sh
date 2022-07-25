fin_dir=eAu_10_100_qhat0_nlo/logs
fout=kdebug_bin.txt
#fin_dir=test_dir
#fout=test_dir/kdebug_bin_test.txt
rm $fout

n_files=`ls $fin_dir/*.log | wc -l`
echo "reading through all $n_files $fin_dir/*.log files for @kdebug statements ... "

n=1
for fin in $fin_dir/*.log; do

  #n_lines=`wc -l $fin | awk '{ print $1 }'`
  #n=0

  while read line; do
    # if line of file starts with "@kdebug " then prints whole line to $fout
    first_word=`echo $line | awk '{print $1;}'`
    if [[ $first_word == "@kdebug" ]]
    then
      echo $line >> $fout
    fi

    # counts through lines of file and prints progress to screen
    #if [[ $(($n%10)) == 0 ]]
    #then
    #  echo "$n/$n_lines"
    #fi
    #n=$((n+1))

  done < $fin

  # counts through number of files read and print progress to screen
  echo "$n/$n_files"
  n=$((n+1))

done
