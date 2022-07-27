fin_dir=/eic/u/kdevereaux/eAu_10_100_qhat0_nlo/logs
fout=kdebug_bin.txt
rm $fout

for fin in $fin_dir/kdebug_bin*; do

  cat $fin >> $fout

done
