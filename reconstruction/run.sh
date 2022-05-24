EP_DIR=./BeAGLE_v101/ep_10_100_baseline_parp2
EA_DIR=./BeAGLE_v102/eAu_10_100

root -l -q "ana_hadron_gen.C(\"'$EA_DIR'/outForPythiaMode/merged.root\", \"'$EA_DIR'/outForPythiaMode/ana_merged.root\")"

root -l -q "plot_chadron_gen.C(\"'$EP_DIR'/outfiles/ana_merged.root\", 0, 1, \"'$EA_DIR'/outForPythiaMode/ana_merged.root\", 1, 1)"

cp topng.sh $EA_DIR/figs
bash $EA_DIR/figs/topng.sh
rm $EA_DIR/figs/topng.sh
