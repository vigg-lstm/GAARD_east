#!/bin/bash

pops=(
	Moshi.arabiensis.Delta
	Moshi.arabiensis.PM
	Muleba.arabiensis.Delta
)
numpops=${#pops[@]}

for chrom in 2L 2R 3L 3R X
do
	for i in $(seq 0 $(($numpops-1)))
	do
		for j in $(seq $(($i+1)) $(($numpops-1)))
		do
			python H1x_tanzania.py ${pops[$i]} ${pops[$j]} $chrom 2000 &
		done
	done
	wait
done
