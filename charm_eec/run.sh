
DIR=./analysis/ep_18_275_D0decay_pow1_inclusive

#rm $DIR/merged.root
#hadd $DIR/merged.root $DIR/*.root
root -l -q "plot_charm_eec_hists.C(\"$DIR/merged.root\", \"$DIR/\")"

DIR=./analysis/ep_18_275_D0decay_pow1_D0injet

#rm $DIR/merged.root
#hadd $DIR/merged.root $DIR/*.root
root -l -q "plot_charm_eec_hists.C(\"$DIR/merged.root\", \"$DIR/\")"