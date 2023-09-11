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
		python PBS.py $pop $chrom 500 &
		# We wait a few minutes between runs so that they're not all downloading at once
		sleep 20m
	done
	wait
done
