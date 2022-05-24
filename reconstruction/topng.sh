# converts all pdfs in containing folder into pngs

for pdffile in ./*
do
	echo "converting ${pdffile}"
	pdftoppm "${pdffile}" "${pdffile%.*}" -png -x 500 -q
done
