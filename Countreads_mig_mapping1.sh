#!/bin/bash
#作業ディレクトリの指定
work_dir=/Users/masumi/work
#リードのあるディレクトリの指定
read_dir=/Users/masumi/migseq/Acropora/raw

#ファイル解凍
#リードからリストを作成
cd $read_dir
gunzip *

#出力ディレクトリの作成
cd $work_dir
mkdir filter yuasa1 yuasa2 remove form trim cat

#リードの処理
#forward側の処理
cd $read_dir
ls *1_001.fastq > ${work_dir}/filelist1.txt
cd $work_dir
while read line
do
    FASTA=`echo "$line" | sed s/fastq/fasta/g`
    fastq_quality_filter -v -Q 33 -q 30 -p 40 -i ${read_dir}/$line -o filter/$line
    cutadapt -e 0.05 -b GTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ./filter/$line > yuasa1/$FASTA #Read1seq.fasta
    fastq_to_fasta -i yuasa1/$FASTA -o yuasa2/$FASTA
    python /usr/local/bin/short_read_remover.py yuasa2/$FASTA 80 > remove/$FASTA
    fasta_formatter -i remove/$FASTA -o form/$FASTA
done < filelist1.txt 

#reverse側の処理
cd $read_dir
ls *2_001.fastq > ${work_dir}/filelist2.txt
cd $work_dir
while read line
do
    FASTA=`echo "$line" | sed s/fastq/fasta/g`
    fastx_trimmer -Q 33 -f 15 -i ${read_dir}/$line -o trim/$line
    fastq_quality_filter -v -Q 33 -q 30 -p 40 -i trim/$line -o filter/$line
    cutadapt -e 0.05 -b CAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAC ./filter/$line > yuasa1/$FASTA #Read2seq.fasta
    fastq_to_fasta -i yuasa1/$FASTA -o yuasa2/$FASTA
    python /usr/local/bin/short_read_remover.py yuasa2/$FASTA 80 > remove/$FASTA
    fasta_formatter -i remove/$FASTA -o form/$FASTA
done < filelist2.txt

#リード数確認用リスト作成
while read line
do
    FASTA=`echo "$line" | sed s/fastq/fasta/g`
    FASTA2=`echo "$FASTA" | sed s/_R1_/_R2_/g`
    SAMPLE=`echo "$line" | cut -d'_' -f1`
    cat ./form/$FASTA ./form/$FASTA2 > cat/${SAMPLE}.fasta
done < filelist1.txt

#リード数確認
cd $read_dir
wc -l *R1_001.fastq > ${work_dir}/wc-raw
cd $work_dir
cd filter
wc -l * > ../wc-filter
cd $work_dir
cd cat
wc -l * > ../wc-cat
cd $work_dir

