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
		python H12_tanzania.py $pop $chrom 2000 &
	done
	wait
done
