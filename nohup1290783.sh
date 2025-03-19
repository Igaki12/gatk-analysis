#!/bin/bash

# スクリプトの説明
# このスクリプトはGATK解析を実行するためのものです。
# Dutch rabbit SRR1290783 10sample分含まれているらしい


s="SRR1290783"
index_dir="."

# 必要なインデックスファイルの確認
index_files=("oryCun2.fa" "oryCun2.fa.fai" "oryCun2.dict")
for file in "${index_files[@]}"; do
    if [ ! -f "${index_dir}/${file}" ]; then
        echo "Error: Required index file ${file} not found in ${index_dir}."
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Required index file ${file} not found in ${index_dir}."
        fi
        exit 1
    fi
done

# モジュールのロード
modules=("fastqc" "bwa" "samtools" "picard" "gatk")
for m in "${modules[@]}"; do
    module load $m
    if [ $? -ne 0 ]; then
        echo "Error: Failed to load $m"
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Failed to load $m"
        fi
        exit 1
    fi
done

# 必要な入力ファイルの確認
input_files=("${s}_1.fastq" "${s}_2.fastq")
for file in "${input_files[@]}"; do
    if [ ! -s "$file" ]; then
        echo "Error: Input file $file is missing or empty."
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Input file $file is missing or empty."
        fi
        exit 1
    fi
done

# BWAによるマッピング
bwa mem -t 16 -R "@RG\tID:${s}\tSM:Dutch" ${index_dir}/oryCun2.fa ${s}_1.fastq ${s}_2.fastq > ${s}.sam
if [ $? -ne 0 ]; then echo "Error: BWA mapping failed";if [ -f ~/automessage.sh ]; then source ~/automessage.sh "Error: BWA mapping failed";
exit 1; fi

# SAMをBAMに変換
samtools view -bS ${s}.sam > ${s}.bam
if [ $? -ne 0 ]; then 
    echo "Error: SAM to BAM conversion failed"
    if [ -f ~/automessage.sh ]; then 
        source ~/automessage.sh "Error: SAM to BAM conversion failed"
    fi
    exit 1
fi

# BAMのソート
samtools sort -T ${s}.sort -o ${s}_sort.bam ${s}.bam
if [ $? -ne 0 ]; then echo "Error: BAM sorting failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: BAM sorting failed"
fi
exit 1
fi

# BAMのインデックス作成
samtools index ${s}_sort.bam
if [ $? -ne 0 ]; then echo "Error: BAM indexing failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: BAM indexing failed"
fi
exit 1
fi

# 重複のマーク
picard MarkDuplicates -I ${s}_sort.bam -M ${s}_marked_dup_metrics -O ${s}_md.bam --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --SORTING_COLLECTION_SIZE_RATIO 0.1
if [ $? -ne 0 ]; then echo "Error: MarkDuplicates failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: MarkDuplicates failed"
fi
exit 1
fi

# 重複の除去
samtools view -bh -F 1036 ${s}_md.bam > ${s}_rmdup.bam
if [ $? -ne 0 ]; then 
    echo "Error: Removing duplicates failed"
    if [ -f ~/automessage.sh ]; then 
        source ~/automessage.sh "Error: Removing duplicates failed"
    fi
    exit 1
fi

# BAMのインデックス作成
samtools index ${s}_rmdup.bam
if [ $? -ne 0 ]; then echo "Error: BAM indexing failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: BAM indexing failed"
fi
exit 1
fi

# GATK HaplotypeCaller
gatk --java-options "-Xmx32g" HaplotypeCaller -R ${index_dir}/oryCun2.fa -I ${s}_rmdup.bam -O ${s}.vcf
if [ $? -ne 0 ]; then echo "Error: HaplotypeCaller failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: HaplotypeCaller failed"
fi
exit 1
fi

# バリアントの選択（SNPとINDEL）
gatk --java-options "-Xmx32g" SelectVariants -R ${index_dir}/oryCun2.fa -V ${s}.vcf -select-type SNP -O ${s}_snp.vcf
if [ $? -ne 0 ]; then echo "Error: SelectVariants (SNP) failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: SelectVariants (SNP) failed"
fi
exit 1
fi

gatk --java-options "-Xmx32g" SelectVariants -R ${index_dir}/oryCun2.fa -V ${s}.vcf -select-type INDEL -O ${s}_indel.vcf
if [ $? -ne 0 ]; then echo "Error: SelectVariants (INDEL) failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: SelectVariants (INDEL) failed"
fi
exit 1
fi

# バリアントフィルタリング（SNP）
gatk --java-options "-Xmx32g" VariantFiltration -R ${index_dir}/oryCun2.fa -V ${s}_snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" -O ${s}_snp_filtered.vcf
if [ $? -ne 0 ]; then echo "Error: VariantFiltration (SNP) failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: VariantFiltration (SNP) failed"
fi
exit 1
fi

# バリアントフィルタリング（INDEL）
gatk --java-options "-Xmx32g" VariantFiltration -R ${index_dir}/oryCun2.fa -V ${s}_indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "INDEL_FILTER" -O ${s}_indel_filtered.vcf
if [ $? -ne 0 ]; then echo "Error: VariantFiltration (INDEL) failed"
if [ -f ~/automessage.sh ]; then 
    source ~/automessage.sh "Error: VariantFiltration (INDEL) failed"
fi
exit 1
fi

echo "Pipeline completed successfully! Sample: ${s}"
# source ~/autoreminder.sh を、~/autoreminder.sh があれば実行する
if [ -f ~/automessage.sh ]; then
    source ~/automessage.sh $s
fi
