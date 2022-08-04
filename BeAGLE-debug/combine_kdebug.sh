fin_dir=/eic/u/kdevereaux/work/BeAGLE-debug/eAu_10_100_qhat0_nlo_run2/logs
fout=/eic/u/kdevereaux/work/BeAGLE-debug/eAu_10_100_qhat0_nlo_run2/kdebug_bin.txt
rm $fout

for fin in $fin_dir/kdebug_bin*; do

  cat $fin >> $fout

done
