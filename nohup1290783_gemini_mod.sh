#!/bin/bash

# スクリプトの説明
# このスクリプトはGATK解析を実行するためのものです。
# Dutch rabbit SRR1290783 10sample分含まれているらしい

s="SRR1290783"
index_dir="." # 参照ゲノムとインデックスファイルがあるディレクトリ (カレントディレクトリ)

# --- ここから修正 ---

# 必要なインデックスファイルの確認
# GATK等で必要なファイル + BWAの主要インデックスファイル
index_files=(
    "oryCun2.fa"
    "oryCun2.fa.fai"
    "oryCun2.dict"
    "oryCun2.fa.amb" # BWA index
    "oryCun2.fa.ann" # BWA index
    "oryCun2.fa.bwt" # BWA index
    "oryCun2.fa.pac" # BWA index
    "oryCun2.fa.sa"  # BWA index
)

for file in "${index_files[@]}"; do
    if [ ! -f "${index_dir}/${file}" ]; then
        echo "Error: Required index file ${file} not found in ${index_dir}."
        # BWAインデックスファイルが見つからない場合、作成方法のヒントを表示
        if [[ "$file" == *".fa.amb" || "$file" == *".fa.ann" || "$file" == *".fa.bwt" || "$file" == *".fa.pac" || "$file" == *".fa.sa" ]]; then
            echo "This is a BWA index file. Please make sure you have run 'bwa index ${index_dir}/oryCun2.fa'."
        fi
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Required index file ${file} not found in ${index_dir}."
        fi
        exit 1
    fi
done

# --- 修正ここまで ---

# モジュールのロード
modules=("fastqc" "bwa" "samtools" "picard" "gatk")
for m in "${modules[@]}"; do
    # module load を試みる前に存在確認を行う方がより親切かもしれないが、元の挙動を維持
    module load $m
    if [ $? -ne 0 ]; then
        echo "Error: Failed to load module $m"
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Failed to load module $m"
        fi
        exit 1
    fi
done

# 必要な入力ファイルの確認
input_files=("${s}_1.fastq" "${s}_2.fastq")
for file in "${input_files[@]}"; do
    if [ ! -s "$file" ]; then # -s はファイルが存在し、かつサイズが0より大きいことを確認
        echo "Error: Input file $file is missing or empty."
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Input file $file is missing or empty."
        fi
        exit 1
    fi
done

# --- 以降の処理は変更なし ---

# BWAによるマッピング
echo "Starting BWA mapping for ${s}..."
bwa mem -t 16 -R "@RG\tID:${s}\tSM:Dutch" ${index_dir}/oryCun2.fa ${s}_1.fastq ${s}_2.fastq > ${s}.sam
if [ $? -ne 0 ]; then
    echo "Error: BWA mapping failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: BWA mapping failed for ${s}"
    fi
    exit 1
fi
echo "BWA mapping completed."

# SAMをBAMに変換
echo "Converting SAM to BAM for ${s}..."
samtools view -bS ${s}.sam > ${s}.bam
if [ $? -ne 0 ]; then
    echo "Error: SAM to BAM conversion failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: SAM to BAM conversion failed for ${s}"
    fi
    # 中間ファイル削除の考慮 (オプション)
    # rm -f ${s}.sam
    exit 1
fi
rm -f ${s}.sam # 成功したらSAMファイルを削除
echo "SAM to BAM conversion completed."

# BAMのソート
echo "Sorting BAM for ${s}..."
samtools sort -@ 4 -T ${s}.sort_temp -o ${s}_sort.bam ${s}.bam # -@でスレッド数、-Tで一時ファイルプレフィックスを指定
if [ $? -ne 0 ]; then
    echo "Error: BAM sorting failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: BAM sorting failed for ${s}"
    fi
    # 中間ファイル削除の考慮 (オプション)
    # rm -f ${s}.bam
    exit 1
fi
rm -f ${s}.bam # 成功したらソート前BAMファイルを削除
echo "BAM sorting completed."

# BAMのインデックス作成
echo "Indexing sorted BAM for ${s}..."
samtools index ${s}_sort.bam
if [ $? -ne 0 ]; then
    echo "Error: BAM indexing failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: BAM indexing failed for ${s}"
    fi
    exit 1
fi
echo "BAM indexing completed."

# 重複のマーク
echo "Marking duplicates for ${s}..."
# メモリ使用量を考慮して JAVA_OPTS を設定 (環境変数でも可)
export _JAVA_OPTIONS="-Xmx8g" # 例: 8GBに設定 (picardはGATKよりメモリ要求が低い場合が多い)
picard MarkDuplicates I=${s}_sort.bam M=${s}_marked_dup_metrics.txt O=${s}_md.bam ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false TMP_DIR=./picard_tmp # 引数指定方法を修正, VALIDATION_STRINGENCYをLENIENTに, TMP_DIRを指定
if [ $? -ne 0 ]; then
    echo "Error: MarkDuplicates failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: MarkDuplicates failed for ${s}"
    fi
    # 中間ファイル削除の考慮 (オプション)
    # rm -f ${s}_sort.bam ${s}_sort.bam.bai
    exit 1
fi
unset _JAVA_OPTIONS # 設定したJavaオプションを解除
rm -f ${s}_sort.bam ${s}_sort.bam.bai # 成功したら入力ファイルを削除
echo "Marking duplicates completed."

# MarkDuplicatesで重複除去は行わない設定にしたため、このステップは不要 (コメントアウト推奨)
# もしMarkDuplicatesで REMOVE_DUPLICATES=true にしない場合は、以下のフィルタリングが必要
# --- 重複除去 (samtools view を使う場合) ---
# echo "Removing duplicates (using samtools view) for ${s}..."
# samtools view -bh -F 1024 -F 4 ${s}_md.bam > ${s}_rmdup.bam # -F 1024 (重複) と -F 4 (unmapped) を除く
# if [ $? -ne 0 ]; then
#     echo "Error: Removing duplicates failed for ${s}"
#     if [ -f ~/automessage.sh ]; then
#         source ~/automessage.sh "Error: Removing duplicates failed for ${s}"
#     fi
#     exit 1
# fi
# rm -f ${s}_md.bam # 成功したらMarkDuplicatesの出力ファイルを削除
# echo "Removing duplicates completed."
# mv ${s}_md.bam ${s}_rmdup.bam # MarkDuplicatesで除去しない場合は、ファイルをリネームして後続処理に渡す
# -----------------------------------------
# MarkDuplicates(REMOVE_DUPLICATES=false)の出力をそのまま使う
mv ${s}_md.bam ${s}_rmdup.bam
echo "Using MarkDuplicates output (duplicates marked, not removed) as ${s}_rmdup.bam"

# BAMのインデックス作成 (rmdup後)
echo "Indexing final BAM for ${s}..."
samtools index ${s}_rmdup.bam
if [ $? -ne 0 ]; then
    echo "Error: Final BAM indexing failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: Final BAM indexing failed for ${s}"
    fi
    exit 1
fi
echo "Final BAM indexing completed."

# GATK HaplotypeCaller
echo "Running GATK HaplotypeCaller for ${s}..."
gatk --java-options "-Xmx32g" HaplotypeCaller \
    -R ${index_dir}/oryCun2.fa \
    -I ${s}_rmdup.bam \
    -O ${s}_raw.vcf.gz # 出力を圧縮形式(.gz)にするのを推奨
if [ $? -ne 0 ]; then
    echo "Error: HaplotypeCaller failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: HaplotypeCaller failed for ${s}"
    fi
    exit 1
fi
echo "HaplotypeCaller completed."

# バリアントの選択（SNPとINDEL）
echo "Selecting SNPs for ${s}..."
gatk --java-options "-Xmx4g" SelectVariants \
    -R ${index_dir}/oryCun2.fa \
    -V ${s}_raw.vcf.gz \
    --select-type-to-include SNP \
    -O ${s}_snp.vcf.gz # 出力を圧縮形式(.gz)にするのを推奨
if [ $? -ne 0 ]; then
    echo "Error: SelectVariants (SNP) failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: SelectVariants (SNP) failed for ${s}"
    fi
    exit 1
fi
echo "SNP selection completed."

echo "Selecting INDELs for ${s}..."
gatk --java-options "-Xmx4g" SelectVariants \
    -R ${index_dir}/oryCun2.fa \
    -V ${s}_raw.vcf.gz \
    --select-type-to-include INDEL \
    -O ${s}_indel.vcf.gz # 出力を圧縮形式(.gz)にするのを推奨
if [ $? -ne 0 ]; then
    echo "Error: SelectVariants (INDEL) failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: SelectVariants (INDEL) failed for ${s}"
    fi
    exit 1
fi
echo "INDEL selection completed."
# 元のraw VCFは不要になったら削除
rm -f ${s}_raw.vcf.gz ${s}_raw.vcf.gz.tbi

# バリアントフィルタリング（SNP） - GATK推奨のハードフィルタリング条件例
echo "Filtering SNPs for ${s}..."
gatk --java-options "-Xmx4g" VariantFiltration \
    -R ${index_dir}/oryCun2.fa \
    -V ${s}_snp.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
    --filter-name "SNP_HardFilter" \
    -O ${s}_snp_filtered.vcf.gz # 出力を圧縮形式(.gz)にするのを推奨
if [ $? -ne 0 ]; then
    echo "Error: VariantFiltration (SNP) failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: VariantFiltration (SNP) failed for ${s}"
    fi
    exit 1
fi
echo "SNP filtering completed."
rm -f ${s}_snp.vcf.gz ${s}_snp.vcf.gz.tbi # フィルタ前ファイルを削除

# バリアントフィルタリング（INDEL） - GATK推奨のハードフィルタリング条件例
echo "Filtering INDELs for ${s}..."
gatk --java-options "-Xmx4g" VariantFiltration \
    -R ${index_dir}/oryCun2.fa \
    -V ${s}_indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
    --filter-name "INDEL_HardFilter" \
    -O ${s}_indel_filtered.vcf.gz # 出力を圧縮形式(.gz)にするのを推奨
if [ $? -ne 0 ]; then
    echo "Error: VariantFiltration (INDEL) failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: VariantFiltration (INDEL) failed for ${s}"
    fi
    exit 1
fi
echo "INDEL filtering completed."
rm -f ${s}_indel.vcf.gz ${s}_indel.vcf.gz.tbi # フィルタ前ファイルを削除

echo "Pipeline completed successfully for sample: ${s}"
# autoreminder.sh ではなく automessage.sh を使うように修正
if [ -f ~/automessage.sh ]; then
    source ~/automessage.sh "Pipeline completed successfully for sample: ${s}"
fi

exit 0 # 正常終了