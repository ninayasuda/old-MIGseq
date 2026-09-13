#!/bin/bash

#ディレクトリの指定
work_dir=/Users/masumi/work
read_dir=/Users/masumi/migseq/Acropora/raw

cd "$work_dir"

echo -e "Sample\tRaw_R1\tRaw_R2\tFiltered_R1\tFiltered_R2\tFinal_R1_FASTA\tFinal_R2_FASTA\tFinal_cat_FASTA" > read_count_summary.tsv

cd "$read_dir"
ls *_R1_001.fastq > ${work_dir}/filelist1.txt

cd "$work_dir"

while read line
do
    SAMPLE=${line%%_R1_001.fastq}
    CAT_SAMPLE=${SAMPLE%%_S*}

    R1_raw=${read_dir}/${SAMPLE}_R1_001.fastq
    R2_raw=${read_dir}/${SAMPLE}_R2_001.fastq

    R1_filter=${work_dir}/filter/${SAMPLE}_R1_001.fastq
    R2_filter=${work_dir}/filter/${SAMPLE}_R2_001.fastq

    R1_fasta=${work_dir}/form/${SAMPLE}_R1_001.fasta
    R2_fasta=${work_dir}/form/${SAMPLE}_R2_001.fasta

    cat_fasta=${work_dir}/cat/${CAT_SAMPLE}.fasta

    raw_r1=$(seqkit stats -T "$R1_raw" | awk 'NR==2{print $4}')
    raw_r2=$(seqkit stats -T "$R2_raw" | awk 'NR==2{print $4}')

    filt_r1=$(seqkit stats -T "$R1_filter" | awk 'NR==2{print $4}')
    filt_r2=$(seqkit stats -T "$R2_filter" | awk 'NR==2{print $4}')

    final_r1=$(grep -c "^>" "$R1_fasta")
    final_r2=$(grep -c "^>" "$R2_fasta")
    final_cat=$(grep -c "^>" "$cat_fasta")

    echo -e "${SAMPLE}\t${raw_r1}\t${raw_r2}\t${filt_r1}\t${filt_r2}\t${final_r1}\t${final_r2}\t${final_cat}"

done < filelist1.txt >> read_count_summary.tsv

#このコードはchatGPTを使いながら作成しました。