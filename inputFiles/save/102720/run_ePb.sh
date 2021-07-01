#!/usr/bin/bash 

if [ -z "$1" ]
then
	echo "No job number set."
        echo "Please run as ./run_ePb.sh jobnumber"
	echo "Exiting..."
	exit 1
fi

echo "-----------------------------------"
echo "Running BeAGLE Simulation for ePb Collider!!!"
echo "-----------------------------------"
echo "Performing Job $1"
echo "..."
echo ""

VAR1=$1
PythiaCard="s/eAu.txt/ePb_${VAR1}.txt/g" 
BeAGLECard="s/eAS1noq/PyInp_${VAR1}/g"

sed "${PythiaCard}" eAS1noq > PyInp_$1
sed "${BeAGLECard}" ./inputFiles/ePb.inp > ./inputFiles/ePb_$1.inp

#VAR2=$((5 * ${VAR1}))
#echo "Sleeping for ${VAR2} Seconds"
#sleep ${VAR2}

SEED=`od -vAn -N2 -tu2 < /dev/urandom`
echo "The Random SEED is ${SEED// /}"
echo ""
sed -i "s/1234567/${SEED// /}/g" PyInp_$1


echo "Running BeAGLE..."
$BEAGLESYS/BeAGLE < inputFiles/ePb_$1.inp > logs/ePb_$1.log 

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q "make_tree.C(\"ePb_$1.txt\")"
echo "Done!!!"
echo ""

echo "Cleaning up..."
rm -vf ./outForPythiaMode/ePb_$1.txt
rm -vf PyInp_$1
rm -vf inputFiles/ePb_$1.inp
echo "Done!!!"
echo ""



