#! /bin/bash

# This script takes an input file path as an argument, and runs a python script to 
# process the input file and write an output ROOT file.
# The main use is to give this script to a slurm script.

# Take two command line arguments -- (1) input file path, (2) output dir prefix
if [ "$1" != "" ]; then
  INPUT_FILE=$1
  #echo "Input file: $INPUT_FILE"
else
  echo "Wrong command line arguments"
fi

if [ "$2" != "" ]; then
  JOB_ID=$2
  echo "Job ID: $JOB_ID"
else 
  echo "Wrong command line arguments"
fi

if [ "$3" != "" ]; then
  TASK_ID=$3
  echo "Task ID: $TASK_ID"
else
  echo "Wrong command line arguments"
fi

# Define output path from relevant sub-path of input file
OUTPUT_PREFIX="AnalysisResults/wenqing/$JOB_ID"
# Note: suffix depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f5-11)
#echo $OUTPUT_SUFFIX
OUTPUT_DIR="/rstorage/alice/$OUTPUT_PREFIX/$OUTPUT_SUFFIX"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# Load modules
module use /software/users/ploskon/heppy/modules
module load heppy/1.0
module use /software/users/wenqing/pyjetty/modules
module load pyjetty/1.0
module list
echo $PYJETTY_DIR

# Run python script via pipenv
cd /software/users/wenqing/pyjetty/pyjetty/alice_analysis
#pipenv run python process/user/wenqing/process_data_ENC.py -c config/ENC/pp/process_pp.yaml -f $INPUT_FILE -o $OUTPUT_DIR
python process/user/wenqing/process_data_ENC.py -c config/ENC/pp/process_pp.yaml -f $INPUT_FILE -o $OUTPUT_DIR

# Move stdout to appropriate folder
mv /rstorage/alice/AnalysisResults/wenqing/slurm-${JOB_ID}_${TASK_ID}.out /rstorage/alice/AnalysisResults/wenqing/${JOB_ID}/
