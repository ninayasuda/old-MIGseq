MIG-seq Processing & SNP Calling Pipeline
A two-stage shell workflow for quality control, read mapping, and variant detection with Stacks

English README (scroll further down for Japanese)
1. Quick Summary
This pipeline automates a typical MIG-seq workflow:

Read QC & trimming (mig_mapping1.sh, first half)

Read-pair repair, reference mapping, BAM processing (mig_mapping1.sh, second half)

Stacks 2.2 SNP calling & export (mig_mapping2.sh)

The scripts assume Illumina R1/R2 FASTQ files produced on a MiSeq and a draft reference genome.

2. Software Requirements
Tool	Tested version	Purpose
FASTX-Toolkit (fastq_quality_filter, fastx_trimmer)	0.0.14	read QC & fixed-length trimming
cutadapt	≥ 1.13	adapter removal
BWA-MEM	≥ 0.7.17	read mapping
samtools	≥ 1.10	BAM conversion/sorting/indexing
BBMap repair.sh	≥ 38	restore read pairing
Stacks	2.2	gstacks, populations
Python 3	≥ 3.8	helper scripts list_maker_1st.py, short_read_remover.py

Add all tools to $PATH or specify absolute paths inside the scripts.

3. Directory Layout & Variables
Edit the variables at the top of each script:

Variable	Meaning	Example
work_dir	main working directory	/Users/heliopora/work
read_dir	raw FASTQ location	/Users/heliopora/migseq/Gizeru/raw
genome	reference FASTA	/Users/heliopora/migseq/refGenome/Seashell/Estearnsii_consensus500.fa
popmap_file	Stacks population map	/Users/heliopora/migseq/Gizeru/no_utsumi_popmap.txt
t	threads used by BWA/Stacks	4

Intermediate folders (i/ ii/ iii/ trimed/, repair/, mapping/, etc.) are created automatically and deleted at the end to save disk space.

4. Usage
bash
コピーする
編集する
# 1. Make scripts executable
chmod +x mig_mapping1.sh mig_mapping2.sh

# 2. Run QC + Mapping
bash mig_mapping1.sh      # creates trimmed FASTQ, BAM, BAM index

# 3. Run SNP calling
bash mig_mapping2.sh      # produces VCF & STRUCTURE files in ./stacks
Tip: run each stage inside a screen/tmux session or HPC job scheduler for long datasets.

5. Output Overview
Path	Contents
trimed/	quality-filtered & adapter-free FASTQ pairs
repair/	read pairs after BBMap repair.sh
mapping/	sorted/indexed BAM files (*.bam, *.bai)
stacks/	Stacks catalog, *.vcf, *.structure, logs

6. Custom Helper Scripts
Script	Role
list_maker_1st.py	converts a list of R1 filenames into sample basenames
short_read_remover.py	discards reads shorter than a user-defined threshold (default 80 bp)

Place them in the directory indicated in each script (/Users/heliopora/local/bin/).

7. Tuning & Troubleshooting
Scenario	What to change
Different adapters	edit the two cutadapt -b sequences in R1 and R2 sections
Lenient quality filter	lower -q 30 -p 40 thresholds in fastq_quality_filter
Long reference build time	comment out bwa index after the first successful run
Stacks population thresholds	change -r, --max_obs_het, --min_maf options in populations

Check stderr logs; most issues come from wrong file paths or missing tools.

8. Attribution
Original scripts by Hideaki Yuasa & Hiroki Taninaka (Nov 13 2019). Minor edits and this README prepared by ChatGPT for Nina Yasuda’s project.

日本語 README
1. 概要
本パイプラインは MIG-seq データの解析を 2 段階で自動化します。

リード QC・トリミング（mig_mapping1.sh 前半）

リード修復・リファレンスマッピング・BAM 整形（同スクリプト後半）

Stacks 2.2 による SNP コール・書き出し（mig_mapping2.sh）

MiSeq 由来のペアエンド FASTQ とリファレンスゲノムを前提にしています。

2. 必要ソフトウェア
ソフト	動作確認版	用途
FASTX-Toolkit	0.0.14	fastq_quality_filter, fastx_trimmer
cutadapt	1.13 以上	アダプタ除去
BWA-MEM	0.7.17 以上	マッピング
samtools	1.10 以上	BAM 変換・ソート・インデックス
BBMap repair.sh	38 系	ペアリード修復
Stacks	2.2	gstacks, populations
Python 3	3.8 以上	補助スクリプト実行

各ツールは $PATH に通すか、スクリプト内に絶対パスで指定してください。

3. ディレクトリと変数
スクリプト冒頭の変数を自分の環境に合わせて変更します。

変数	役割	例
work_dir	作業用トップディレクトリ	/Users/heliopora/work
read_dir	生 FASTQ 置き場	/Users/heliopora/migseq/Gizeru/raw
genome	参照ゲノム FASTA	/Users/heliopora/migseq/refGenome/Seashell/Estearnsii_consensus500.fa
popmap_file	Stacks 用 popmap	/Users/heliopora/migseq/Gizeru/no_utsumi_popmap.txt
t	使用スレッド数	4

中間フォルダは自動生成後に削除され、容量を節約します。

4. 使い方
bash
コピーする
編集する
# 実行権付与
chmod +x mig_mapping1.sh mig_mapping2.sh

# ① QC とマッピング
bash mig_mapping1.sh

# ② SNP コール
bash mig_mapping2.sh
データ量が多い場合は screen/tmux かジョブスケジューラでの実行を推奨します。

5. 主な出力
パス	内容
trimed/	QC 済み FASTQ
repair/	repair.sh 後のペアリード
mapping/	整形済み BAM (*.bam, *.bai)
stacks/	Stacks カタログ、VCF、STRUCTURE 形式ファイル

6. 補助 Python スクリプト
スクリプト	役割
list_maker_1st.py	R1 ファイル名の一覧からサンプル名を抽出
short_read_remover.py	長さ閾値（標準 80 bp）未満のリードを除外

指定パス（/Users/heliopora/local/bin/）に配置してください。

7. よくある調整ポイント
目的	変更箇所
別アダプタ配列	cutadapt -b の配列を編集
QC の厳しさ緩和	fastq_quality_filter の -q -p 値を下げる
インデックス作成を省略	1 回目完了後に bwa index 行をコメントアウト
Stacks パラメータ調整	populations の -r --max_obs_het --min_maf 変更

エラー時はパス・ツールの存在をまず確認してください。

