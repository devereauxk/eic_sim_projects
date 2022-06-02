# removes all .txt files in all outForPythiaMode subfolders of all folders in current dirrectory

for folder in ./*
do
	echo "${folder}"
	if [[ -d "${folder}" ]]
	then
		echo "is folder"
		rm outForPythiaMode/*.txt
	fi
done
