#!/usr/bin/bash 

if [ -z "$1" ]
then
	echo "No job number set."
        echo "Please run as ./run_minbias.sh jobnumber"
	echo "Exiting..."
	exit 1
fi

chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

INPUT=$(( 100 + $1 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR
echo start running in directory $DIR

echo "-----------------------------------"
echo "Running PYTHIA Simulation for ep Collider!!!"
echo "-----------------------------------"
echo "$2 GeV e + $3 GeV p"
echo "Performing Job $1"
echo "Generate $4 Events with PythiaeRHIC"
echo "..."
echo ""

E_electron=$2
E_proton=$3
NEVT=$4

echo "Running PYTHIA..."
cp /gpfs/mnt/gpfs02/eic/wfan/event_gen/pythiaeRHIC/${E_electron}_${E_proton}.inp input_minbias.inp
sed -i '4s/.*/'$NEVT',10           ! Number of events/' input_minbias.inp

pythiaeRHIC < input_minbias.inp

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
ln -s /gpfs/mnt/gpfs02/eic/wfan/HF_vtx_study/make_tree.C .
root -l -b -q "make_tree.C(\"ep_minbias.out\")"
echo "Done!!!"
echo ""

echo "Cleaning up..."
#mv -v ep_minbias.out /gpfs/mnt/gpfs02/eic/wfan/event_gen/pythiaeRHIC/${E_electron}_${E_proton}/ep_minbias_${INPUT}.out
mv -v ep_minbias.root /gpfs/mnt/gpfs02/eic/wfan/event_gen/pythiaeRHIC/${E_electron}_${E_proton}/ep_minbias_${INPUT}.root
echo "Done!!!"

popd


