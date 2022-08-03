fin_dir=/eic/u/kdevereaux/work/BeAGLE-debug/eAu_10_100_tauforOff_qhat0_nlo/logs
fout=/eic/u/kdevereaux/work/BeAGLE-debug/eAu_10_100_tauforOff_qhat0_nlo/kdebug_bin.txt
rm $fout

for fin in $fin_dir/kdebug_bin*; do

  cat $fin >> $fout

done
