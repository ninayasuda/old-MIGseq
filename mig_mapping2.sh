#!/bin/bash

#作業ディレクトリの指定
work_dir=/Users/heliopora/work

#ファイルサイズがあるサンプルのみのpopulation-information
popmap_file=/Users/heliopora/migseq/Gizeru/no_utsumi_popmap.txt

#使用するthread数
t=4

#Stacks2.2解析
gstacks -I ./mapping/ -M $popmap_file -O ./stacks/ -t $t
populations -P ./stacks/ -M $popmap_file -r 0.7 --max_obs_het 0.7 --min_maf 0.01 --vcf --structure -t $t 

#中間ディレクトリを削除
rm -r repair &
rm -r mapping &
wait

#made_by_HideakiYUASA&HirokiTANINAKA_November13th2019



