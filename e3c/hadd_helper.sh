TARGET_DIR=$1

cd $TARGET_DIR

hadd -j merged_1.root *_1*
hadd -j merged_2.root *_2*
hadd -j merged_3.root *_3*
hadd -j merged_4.root *_4*
hadd -j merged_5.root *_5*
hadd -j merged_6.root *_6*
hadd -j merged_7.root *_7*
hadd -j merged_8.root *_8*
hadd -j merged_9.root *_9*

hadd -j merged.root merged_temp_*

rm merged_temp_*
