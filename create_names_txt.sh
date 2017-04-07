#!/bin/bash  

BAM=$1
i=1
for cur_bam in $BAM/*.bam
do
  # use bam file name as sample name
  bam_file_name=$(basename "${cur_bam%.*}")
  # remove whitespaces from name
  SM1="$(echo -e "${bam_file_name}" | tr -d '[[:space:]]')"

  # extract sample name from bam file read group info field
  SM2=$(samtools view -H $cur_bam | grep "^@RG" | tail -n1 | sed "s/.*SM:\\([^	]*\\).*/\\1/" | tr -d '[:space:]')

  printf "$SM1	$SM2\\n" >> names.txt
  i=$((i+1))
done
