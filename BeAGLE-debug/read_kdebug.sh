fin='eAu_10_100_qhat0_nlo/logs/eAu_0.log'

echo "reading through $fin for @kdebug statements ... "

while read line; do
  first_word=`echo $line | awk '{print $1;}'`
  if [[ $first_word == "@kdebug" ]]
  then
    echo $line
  fi
done < $fin
