TARGET_DIR=$1

cd $TARGET_DIR

hadd -j merged_temp_1.root *_1*
hadd -j merged_temp_2.root *_2*
hadd -j merged_temp_3.root *_3*
hadd -j merged_temp_4.root *_4*
hadd -j merged_temp_5.root *_5*
hadd -j merged_temp_6.root *_6*
hadd -j merged_temp_7.root *_7*
hadd -j merged_temp_8.root *_8*
hadd -j merged_temp_9.root *_9*

hadd -j merged.root merged_temp_*

rm merged_temp_*
