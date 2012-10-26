#!/bin/bash
# ##################################################################
# Script for rigid movement compensation on phallusya experiments
# Oscar Esteban - oesteban@die.upm.es
# v. 1.0 - 2011-07-07 00:15
####################################################################

start_t=$1
stop_t=$2
data_fld=${3%/}
result_fld=${4%/}
channel=$5


previous=""

procesed_files=( )
last=0

mkdir -p $result_fld

for file in $(find ${data_fld}/ -name *_ch0${channel}.vtk | sort -u)
do 
	tstep=$( echo $file | sed -e 's/.*_t\([0-9]\{3\}\)_.*/\1/')

	clean_t=$( echo $((10#$tstep)) )
#	clean_t=$(expr  ${tstep##+(0)} + 0 )

	if (( $clean_t > $stop_t  || $clean_t < $start_t )); then	
		continue
	fi
	
	if [ -n "$previous" ]
	then
		procesed_files=( ${procesed_files[@]} $file )
		command="autoreg $previous $file --no-output-moving -o ${result_fld}/${tstep}_"; 
		echo $command;
#		eval $command; 
		last=$clean_t
	fi
	previous=$file
done

# Compose final transforms file
tfms_cmd="mkdir -p ${result_fld}/tforms; mv ${result_fld}/*tforms.txt ${result_fld}/tforms/"
logs_cmd="mkdir -p ${result_fld}/logs; mv ${result_fld}/*log*.txt ${result_fld}/logs/"

echo $tfms_cmd
eval $tfms_cmd

echo $logs_cmd
eval $logs_cmd

if [ -e ${result_fld}/tforms.txt ]; then rm ${result_fld}/tforms.txt; fi

for i in $(find ${result_fld}/tforms -name "*tforms.txt" | sort -u)
do 
	eval "cat $i >> ${result_fld}/tforms.tfm"
done


for infile in ${procesed_files[@]}
do
	tstep=$( echo $infile | sed -e 's/.*_t\([0-9]\{3\}\)_.*/\1/' )
	clean_t=$( echo $((10#$tstep)) )
	caseid=$( expr $clean_t - $start_t - 1 )

	outfile=$(echo ${infile##$data_fld[/]} | sed -e 's/^\(.*\)\(_t[0-9]\{3\}.*\)$/\1r\2/' )
	command="batch_transform $infile ${result_fld}/$outfile -F ${result_fld}/tforms.tfm --case-id $caseid -o ${result_fld}/"
	echo $command
	eval $command
done 
