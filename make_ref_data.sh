#!/bin/bash
echo "#define INSTALLDIR \""$PWD/src"\"" > src/installdir.h

# just add the path and does not check the file exists or not
addPath()
{
    echo "#define $1 \""$2"\"" >> src/installdir.h
}

checkfile()
{
    filename=$1
    echo "take the paths from $filename and add them into installdir.h"
    # ignore the comment lines started with # 
    cat $filename |grep -v "^#"| while read -r line
    do
	#ignore empty lines 
	if [[ -n $line ]] 
	then
	    
	    name=`echo $line |cut -d "=" -f1`
	    path=`echo $line |cut -d "=" -f2`
        #echo "$line"
	    addPath $name $path
	fi
    done
}
#check all the fils in reference.txt
checkfile "reference-data.txt"
