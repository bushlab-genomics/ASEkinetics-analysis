#!/bin/sh

file=$1

scol1=`head -1 $file  | xargs -n1 -d$"\t" | grep "File\sformat" -n | cut -d':' -f1 `
scol2=`head -1 $file  | xargs -n1 -d$"\t" | grep "Output\stype" -n | cut -d':' -f1 `
scol3=`head -1 $file  | xargs -n1 -d$"\t" | grep "File\sassembly" -n | cut -d':' -f1 `
scol4=`head -1 $file  | xargs -n1 -d$"\t" | grep "File\sanalysis\stitle" -n | cut -d':' -f1 `
col1=`head -1 $file  | xargs -n1 -d$"\t" | grep "Experiment\starget" -n | cut -d':' -f1 `
#col1=`head -1 $file  | xargs -n1 -d$"\t" | grep "Assay" -n | cut -d':' -f1 `
col3=`head -1 $file  | xargs -n1 -d$"\t" | grep "File\sdownload\sURL" -n | cut -d':' -f1 `
col2=`head -1 $file  | xargs -n1 -d$"\t" | grep "Biological\sreplicate(s)"  -n | cut -d':' -f1 `
col4=`head -1 $file  | xargs -n1 -d$"\t" | grep "Experiment\saccession"  -n | cut -d':' -f1 `
col5=`head -1 $file  | xargs -n1 -d$"\t" | grep "Biosample\sterm\sname"  -n | cut -d':' -f1 `

#generate list for bam files
cat $file |  awk -F"\t" -v svar1="$scol1" -v svar2="$scol2" -v svar3="$scol3" -v svar4="$scol4" \
			'($svar1=="bam") && ($svar2=="alignments") && ($svar3=="GRCh38") && ($svar4~/ENCODE4/) {print}'| \
	awk -F"\t" -v var1="$col1" -v var2="$col2" -v var3="$col3" -v var4="$col4"  -v var5="$col5" \
			'{print $var4,$var1,$var2,$var3,$var5}' | sed 's/-human//g' |
	awk '{print $0, $5"_"$1"-"$2"_rep"$3".bam"}'> ${file%.*}_bam.list

#generate list for bigwig files
cat $file |  awk -F"\t" -v svar1="$scol1" -v svar2="$scol2" -v svar3="$scol3" -v svar4="$scol4" \
			'($svar1=="bigWig") && ($svar2=="signal p-value") && ($svar3=="GRCh38") && ($svar4~/ENCODE4/) {print}'| \
	awk -F"\t" -v var1="$col1" -v var2="$col2" -v var3="$col3" -v var4="$col4" \
			'{print $var4,$var1,$var2,$var3}' | sed 's/-human//g' > ${file%.*}_bigwig.list


file_default=$2
cat $file_default |  awk -F"\t" -v svar1="$scol1" -v svar2="$scol2" -v svar3="$scol3" -v svar4="$scol4"  \
			'($svar1~/bed narrowPeak/)  && ($svar3=="GRCh38") && ($svar4~/ENCODE4/)  {print}'| \
	awk -F"\t" -v var1="$col1" -v var2="$col2" -v var3="$col3" -v var4="$col4" -v var5="$col5" \
			'{print $var4,$var1,$var3,$var5}' | sed 's/-human//g' | 
	awk '{print $0, $4"_"$1"-"$2".bed"}'> ${file_default%.*}.list



