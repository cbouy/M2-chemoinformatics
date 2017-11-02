#!/bin/bash

# Reads the Tools/Predict output from jatoonSOM for Kohonen or CPG maps
# Outputs the number of classes correctly predicted
# Max number of classes for CPG: 5
# If you have more classes, just add it to the 'arr[number]="letter";' part of the code

usage() {
echo "\
Usage: $0 csvFile mapType
mapType:	K for Kohonen
        	C for CPG"
}

if [ "$#" -ne 2 ]
then
	usage
else
	if [ "$2" = "K" ]
	then
		awk '
		BEGIN {
			count=0
		}{
			if ($1 ~ $4) {count+=1}
		} END {
			print count "/" NR " instances ("count/NR*100"%) were correctly predicted"
		}' $1
	
	elif [ "$2" = "C" ]
	then
		awk '
		BEGIN {
			count = 0;
			arr[1]= "A";
			arr[2]= "B";
			arr[3]= "C";
			arr[4]= "D";
			arr[5]= "E";
		}{
			max=0;
			field=0;
			for (i=4;i<=NF;i++){
				if ($i>max) {max=$i;field=i}
			};
			col=field-3;
			if ($1 ~ arr[col]) {count+=1}
		} END {
		print count "/" NR " instances ("count/NR*100"%) were correctly predicted"}' $1
	else
		usage
	fi
fi
