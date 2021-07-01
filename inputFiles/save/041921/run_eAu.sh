#!/usr/bin/bash 

if [ -z "$1" ]
then
	echo "No job number set."
        echo "Please run as ./run_eAu.sh jobnumber"
	echo "Exiting..."
	exit 1
fi

echo "-----------------------------------"
echo "Running BeAGLE Simulation for eAu Collider!!!"
echo "-----------------------------------"
echo "Performing Job $1"
echo "..."
echo ""

VAR1=$1
PythiaCard="s/eAu.txt/eAu_${VAR1}.txt/g" 
BeAGLECard="s/S3ALL002/InpAu_${VAR1}/g"

sed "${PythiaCard}" S3ALL002 > InpAu_$1
sed "${BeAGLECard}" ./inputFiles/eAu.inp > ./inputFiles/eAu_$1.inp

#VAR2=$((5 * ${VAR1}))
#echo "Sleeping for ${VAR2} Seconds"
#sleep ${VAR2}

SEED=`od -vAn -N2 -tu2 < /dev/urandom`
echo "The Random SEED is ${SEED// /}"
echo ""
sed -i "s/1234567/${SEED// /}/g" InpAu_$1


echo "Running BeAGLE..."
$BEAGLESYS/BeAGLE < inputFiles/eAu_$1.inp > logs/eAu_$1.log 

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"eAu_$1.txt\")"
echo "Done!!!"
echo ""

echo "Cleaning up..."
rm -vf ./outForPythiaMode/eAu_$1.txt
rm -vf InpAu_$1
rm -vf inputFiles/eAu_$1.inp
echo "Done!!!"
echo ""



