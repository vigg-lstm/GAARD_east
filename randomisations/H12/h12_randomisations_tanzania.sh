#!/bin/bash

pops=(
	Moshi.arabiensis.Delta
	Moshi.arabiensis.PM
	Muleba.arabiensis.Delta
)

for chrom in 2L 2R 3L 3R X
do
	for pop in ${pops[@]}
	do
		python H12.py $pop $chrom 200 &
	done
	wait
done
