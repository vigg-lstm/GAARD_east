#!/bin/bash

# Could run it by sample set, but there is no reason to expect Moshi PM to be different to Moshi Delta
#pops=(
#	Moshi.arabiensis.Delta
#	Moshi.arabiensis.PM
#	Muleba.arabiensis.Delta
#)
pops=(
	Moshi
	Muleba
)

for chrom in 2L 2R 3L 3R X
do
	for pop in ${pops[@]}
	do
		python H12_tanzania.py $pop $chrom 2000 &
	done
	wait
done
