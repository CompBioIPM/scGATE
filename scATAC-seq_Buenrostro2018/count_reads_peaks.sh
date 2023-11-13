#!/bin/bash

mkdir -p count_reads_peaks_output
bampath=../../input/sc-bams_nodup/
dirlist=(`ls $bampath*.bam`)

for ((i=1; i<=2034; i++)); do
    echo ./count_reads_peaks_output/$(basename ${dirlist[i-1]}).peaks.txt
    bedtools coverage -a ../../input/combined.sorted.merged.bed -b ${dirlist[i-1]} > ./count_reads_peaks_output/$(basename ${dirlist[i-1]}).peaks.txt
done