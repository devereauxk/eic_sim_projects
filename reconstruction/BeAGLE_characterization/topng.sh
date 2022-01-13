for folder in ./*
do
	echo "${folder}"
	if [[ -d "${folder}" ]]
	then
		echo "is folder"
		for pdfile in $folder/figs/*
		do
			pdftoppm "${pdfile}" "${pdfile%.*}" -png -x 500 -q
		done
	fi
done
