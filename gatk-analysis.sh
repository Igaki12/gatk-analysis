# Description: This script is used to run GATK analysis on the data.
#!/bin/bash

# 引数1として、サンプル名を受け取る
# この時に、引数がない場合はエラーを出力して終了する
if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample> <index_dir>"
    exit 1
fi
s=$1
#引数２として、indexのフォルダを受け取る
if [ $# -ne 2 ]; then
    echo "Usage: $0 <sample> <index_dir>"
    exit 1
fi
index_dir=$2
#module load を行う。エラーが出た場合は終了する
modules=("fastqc" "bwa" "samtools" "picard" "gatk")
for m in ${modules[@]}; do
    module load $m
    if [ $? -ne 0 ]; then
        echo "Error: failed to load $m"
        exit 1
    fi
done
# module load fastqc
# module load bwa
# module load samtools
# module load picard
# module load gatk



bwa mem -t 16 -R "@RG\tID:${s}\tSM:Dutch" ${index_dir}/oryCun2.fa ${s}_1.fastq ${s}_2.fastq > ${s}.sam;
samtools view -bS ${s}.sam > ${s}.bam;
samtools sort -T ${s}.sort -o ${s}_sort.bam ${s}.bam;
samtools index ${s}_sort.bam;
picard MarkDuplicates -I ${s}_sort.bam -M ${s}_marked_dup_metrics -O ${s}_md.bam --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --SORTING_COLLECTION_SIZE_RATIO 0.1;
samtools view -bh -F 1036 ${s}_md.bam > ${s}_rmdup.bam;
samtools index ${s}_rmdup.bam;
gatk --java-options "-Xmx32g" HaplotypeCaller -R ${index_dir}/oryCun2.fa -I ${s}_rmdup.bam -O ${s}.vcf;
gatk --java-options "-Xmx32g" SelectVariants -R ${index_dir}/oryCun2.fa -V ${s}.vcf -select-type SNP -O ${s}_snp.vcf;
gatk --java-options "-Xmx32g" SelectVariants -R ${index_dir}/oryCun2.fa -V ${s}.vcf -select-type INDEL -O ${s}_indel.vcf;
gatk --java-options "-Xmx32g" VariantFiltration -R ${index_dir}/oryCun2.fa -V ${s}_snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "SNP_FILTER" -O ${s}_snp_filtered.vcf;
gatk --java-options "-Xmx32g" VariantFiltration -R ${index_dir}/oryCun2.fa -V ${s}_indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "INDEL_FILTER" -O ${s}_indel_filtered.vcf;
# GitHub Copilot: `--filter-expression`に`SOR>10.0`を追加しない理由について説明します。

# `SOR`（Strand Odds Ratio）は、バリアントのストランドバイアスを評価する指標です。`SOR>10.0`というフィルターを追加すると、ストランドバイアスが高いバリアントがフィルタリングされます。しかし、以下の理由でこのフィルターを追加しないことがあります：

# 1. **データの特性**: データセットによっては、ストランドバイアスが高くても信頼できるバリアントが含まれている場合があります。このような場合、`SOR>10.0`でフィルタリングすると有用なバリアントが除外される可能性があります。

# 2. **他のフィルターとのバランス**: 既に設定されている`QD < 2.0`、`FS > 200.0`、`ReadPosRankSum < -20.0`のフィルターが十分に厳しい場合、さらに`SOR>10.0`を追加することで過剰にフィルタリングされるリスクがあります。

# 3. **プロジェクトの要件**: プロジェクトや解析の目的によっては、ストランドバイアスを重視しない場合もあります。このような場合、`SOR`に基づくフィルターを追加する必要はありません。

# したがって、`SOR>10.0`を追加しない理由は、データの特性や他のフィルターとのバランス、プロジェクトの要件に依存します。

