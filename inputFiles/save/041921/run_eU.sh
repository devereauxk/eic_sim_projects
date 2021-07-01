#!/usr/bin/bash 

if [ -z "$1" ]
then
	echo "No job number set."
        echo "Please run as ./run_eU.sh jobnumber"
	echo "Exiting..."
	exit 1
fi

echo "-----------------------------------"
echo "Running BeAGLE Simulation for eU Collider!!!"
echo "-----------------------------------"
echo "Performing Job $1"
echo "..."
echo ""

VAR1=$1
PythiaCard="s/eAu.txt/eU_${VAR1}.txt/g" 
BeAGLECard="s/S3ALL002/InpU_${VAR1}/g"

sed "${PythiaCard}" S3ALL002 > InpU_$1
sed "${BeAGLECard}" ./inputFiles/eU.inp > ./inputFiles/eU_$1.inp

#VAR2=$((5 * ${VAR1}))
#echo "Sleeping for ${VAR2} Seconds"
#sleep ${VAR2}

SEED=`od -vAn -N2 -tu2 < /dev/urandom`
echo "The Random SEED is ${SEED// /}"
echo ""
sed -i "s/1234567/${SEED// /}/g" InpU_$1


echo "Running BeAGLE..."
$BEAGLESYS/BeAGLE < inputFiles/eU_$1.inp > logs/eU_$1.log 

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"eU_$1.txt\")"
echo "Done!!!"
echo ""

echo "Cleaning up..."
rm -vf ./outForPythiaMode/eU_$1.txt
rm -vf InpU_$1
rm -vf inputFiles/eU_$1.inp
echo "Done!!!"
echo ""



