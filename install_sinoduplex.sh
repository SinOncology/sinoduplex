#!/bin/bash
cd thirdparty/samtools/samtools-1.3.1/
make -j 16
cd ../../../

mkdir -p build
cd build && rm -rf *
cmake ..
cmake ..
cd ../
bash ./make_ref_data.sh
cd build
make -j  16
make install
cd ../

echo "generate version.txt file"
head -1 ./reference-data.txt > version.txt
branch=`git rev-parse --abbrev-ref HEAD`
shaID=`git log|head -1|awk '{print $2}'`

echo "#Sinotools_branch: ${branch}" >>version.txt
echo "#Sinotools_shaID: ${shaID}" >>version.txt
cat version.txt

