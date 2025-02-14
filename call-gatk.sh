#!/bin/bash

# スクリプトの説明
# このスクリプトはGATK解析を実行するためのものです。

# 引数チェック
if [ $# -ne 2 ]; then
    echo "Usage: $0 <sample> <index_dir>"
    exit 1
fi

s=$1
index_dir=$2

# 必要なインデックスファイルの確認
index_files=("oryCun2.fa" "oryCun2.fa.fai" "oryCun2.dict")
for file in "${index_files[@]}"; do
    if [ ! -f "${index_dir}/${file}" ]; then
        echo "Error: Required index file ${file} not found in ${index_dir}."
        exit 1
    fi
done

# モジュールのロード
modules=("fastqc" "bwa" "samtools" "picard" "gatk")
for m in "${modules[@]}"; do
    module load $m
    if [ $? -ne 0 ]; then
        echo "Error: Failed to load $m"
        exit 1
    fi
done

# 必要な入力ファイルの確認
input_files=("${s}_1.fastq" "${s}_2.fastq")
for file in "${input_files[@]}"; do
    if [ ! -s "$file" ]; then
        echo "Error: Input file $file is missing or empty."
        exit 1
    fi
done

# BWAによるマッピング
bwa mem -t 16 -R "@RG\tID:${s}\tSM:Dutch" ${index_dir}/oryCun2.fa ${s}_1.fastq ${s}_2.fastq > ${s}.sam
if [ $? -ne 0 ]; then echo "Error: BWA mapping failed"; exit 1; fi

# SAMをBAMに変換
samtools view -bS ${s}.sam > ${s}.bam
if [ $? -ne 0 ]; then echo "Error: SAM to BAM conversion failed"; exit 1; fi

# BAMのソート
samtools sort -T ${s}.sort -o ${s}_sort.bam ${s}.bam
if [ $? -ne 0 ]; then echo "Error: BAM sorting failed"; exit 1; fi

# BAMのインデックス作成
samtools index ${s}_sort.bam
if [ $? -ne 0 ]; then echo "Error: BAM indexing failed"; exit 1; fi

# 重複のマーク
picard MarkDuplicates -I ${s}_sort.bam -M ${s}_marked_dup_metrics -O ${s}_md.bam --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --SORTING_COLLECTION_SIZE_RATIO 0.1
if [ $? -ne 0 ]; then echo "Error: MarkDuplicates failed"; exit 1; fi

# 重複の除去
samtools view -bh -F 1036 ${s}_md.bam > ${s}_rmdup.bam
if [ $? -ne 0 ]; then echo "Error: Removing duplicates failed"; exit 1; fi

# BAMのインデックス作成
samtools index ${s}_rmdup.bam
if [ $? -ne 0 ]; then echo "Error: BAM indexing failed"; exit 1; fi

# GATK HaplotypeCaller
gatk --java-options "-Xmx32g" HaplotypeCaller -R ${index_dir}/oryCun2.fa -I ${s}_rmdup.bam -O ${s}.vcf
if [ $? -ne 0 ]; then echo "Error: HaplotypeCaller failed"; exit 1; fi

# バリアントの選択（SNPとINDEL）
gatk --java-options "-Xmx32g" SelectVariants -R ${index_dir}/oryCun2.fa -V ${s}.vcf -select-type SNP -O ${s}_snp.vcf
if [ $? -ne 0 ]; then echo "Error: SelectVariants (SNP) failed"; exit 1; fi

gatk --java-options "-Xmx32g" SelectVariants -R ${index_dir}/oryCun2.fa -V ${s}.vcf -select-type INDEL -O ${s}_indel.vcf
if [ $? -ne 0 ]; then echo "Error: SelectVariants (INDEL) failed"; exit 1; fi

# バリアントフィルタリング（SNP）
gatk --java-options "-Xmx32g" VariantFiltration -R ${index_dir}/oryCun2.fa -V ${s}_snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" -O ${s}_snp_filtered.vcf
if [ $? -ne 0 ]; then echo "Error: VariantFiltration (SNP) failed"; exit 1; fi

# バリアントフィルタリング（INDEL）
gatk --java-options "-Xmx32g" VariantFiltration -R ${index_dir}/oryCun2.fa -V ${s}_indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "INDEL_FILTER" -O ${s}_indel_filtered.vcf
if [ $? -ne 0 ]; then echo "Error: VariantFiltration (INDEL) failed"; exit 1; fi

echo "Pipeline completed successfully! Sample: ${s}"
# source ~/autoreminder.sh を、~/autoreminder.sh があれば実行する
if [ -f ~/autoreminder.sh ]; then
    source ~/autoreminder.sh
fi
