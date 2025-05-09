#!/bin/bash

#作業ディレクトリの指定
work_dir=/Users/heliopora/work
#リードのあるディレクトリの指定
read_dir=/Users/heliopora/migseq/Gizeru/raw

#ファイル解凍
#リードからリストを作成
cd $read_dir
gunzip *
ls *1_001.fastq > ${work_dir}/read_name_list.tMp
cd $work_dir
python /Users/heliopora/local/bin/list_maker_1st.py read_name_list.tMp > read_name_list.txt

#出力ディレクトリの作成
mkdir i ii iii trimed

#リードの処理
while read line
do
	#forward側
	#low qualityをトリム
	fastq_quality_filter -v -Q 33 -q 30 -p 40 -i ${read_dir}/${line}_L001_R1_001.fastq -o ./ii/${line}_L001_R1_001.fastq
	#index配列のトリム
	cutadapt -e 0.05 -b GTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ./ii/${line}_L001_R1_001.fastq > ./iii/${line}_L001_R1_001.fastq
	#fastaqから短い配列のトリム
	python /Users/heliopora/local/bin/short_read_remover.py ./iii/${line}_L001_R1_001.fastq 80 > ./trimed/${line}.1.fastq &
	#reverse側
	#頭15bpのトリム
　　　　fastx_trimmer -Q 33 -f 0  -i ${read_dir}/${line}_L001_R2_001.fastq -o ./i/${line}_L001_R2_001.fastq
	#元のスクリプト：fastx_trimmer -Q 33 -f 15 -i ${read_dir}/${line}_L001_R2_001.fastq -o ./i/${line}_L001_R2_001.fastq
	#low qualityをトリム
	fastq_quality_filter -v -Q 33 -q 30 -p 40 -i ./i/${line}_L001_R2_001.fastq -o ./ii/${line}_L001_R2_001.fastq
	#index配列のトリム
	cutadapt -e 0.05 -b CAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAC ./ii/${line}_L001_R2_001.fastq > ./iii/${line}_L001_R2_001.fastq
	#fastaqから短い配列のトリム
	python /Users/heliopora/local/bin/short_read_remover.py ./iii/${line}_L001_R2_001.fastq 80 > ./trimed/${line}.2.fastq
	wait
done < read_name_list.txt

#中間ファイルおよびディレクトリの削除
rm *.tMp
rm -r i ii iii

#!/bin/bash

#リファレンスゲノム
genome=/Users/heliopora/migseq/refGenome/Seashell/Estearnsii_consensus500.fa
#使用するthread数
t=4

#BWAのindex作成
bwa index -p $genome $genome

#出力ディレクトリの作成
mkdir repair mapping stacks kakunin

while read line
do
	#pair情報のないリードを除去
	repair.sh -Xmx50g in1=./trimed/${line}.1.fastq in2=./trimed/${line}.2.fastq out1=./repair/${line}.1.fastq out2=./repair/${line}.2.fastq
	#マッピング
	bwa mem -t $t $genome ./repair/${line}.1.fastq ./repair/${line}.2.fastq > ./mapping/${line}.sam
	#ファイルサイズ確認用
	sed '/^@/d' ./mapping/${line}.sam > ./kakunin/${line}.sam
done < read_name_list.txt

#ファイルサイズがあるサンプルのみのリスト作成
ls -l ./kakunin/*.sam | awk '{if($5>0)print $9}' > kakunin_list.txt
awk -F "./kakunin/" '{print $2}' kakunin_list.txt |awk -F ".sam" '{print $1}' > nakamiari_list.txt

#ファイルサイズ確認用のディレクトリ削除
rm -r ./kakunin/

#ファイルサイズあるサンプルのみsamをbamに変換
while read line
do
	#samからbamに変換
	samtools view -@ $t -Sb ./mapping/${line}.sam > ./mapping/${line}.raw.bam
	#samの消去
	rm ./mapping/${line}.sam &
	#sort
	samtools sort -m 5G -@ $t -o ./mapping/${line}.bam ./mapping/${line}.raw.bam
	#raw.bamの消去
	rm ./mapping/${line}.raw.bam &
	#index
	samtools index ./mapping/${line}.bam
	wait
done < nakamiari_list.txt
wait

#samの消去(ファイルサイズ大きいため)
#rm ./mapping/*.sam
#ファイル名の整形
#rm ./mapping/*.raw.bam

#made_by_HideakiYUASA&HirokiTANINAKA_November13th2019

