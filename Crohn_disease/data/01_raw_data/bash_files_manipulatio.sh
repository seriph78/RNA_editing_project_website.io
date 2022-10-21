#!/bin/bash
yourfilenames=`ls ./*.gz`
for eachfile in $yourfilenames
do
    dir_name=$(echo $eachfile | cut -d "_" -f 2)
	if [[ ! -e $dir_name ]]; then
   	mkdir $dir_name
	elif [[ ! -d $dir_name ]]; then
   	echo "$dir_name already exists but is not a directory" 1>&2
	fi
    new_name=$(echo $eachfile | cut -d "_" -f 3)
	cp $eachfile $dir_name"/"$new_name
	if [ $new_name == "genes.tsv.gz" ]; then
		cp $eachfile $dir_name"/features.tsv.gz"
		rm $dir_name"/genes.tsv.gz"
		rm $dir_name"/feature.tsv.gz"
	else 
		cp $eachfile $dir_name"/"$new_name
	fi

done

