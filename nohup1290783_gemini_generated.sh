#!/bin/bash

# スクリプトの説明
# このスクリプトはGATK解析を実行するためのものです。
# Dutch rabbit SRR1290783 10sample分含まれているらしい


s="SRR1290783"
index_dir="." # 現在のディレクトリを指す
ref_genome="oryCun2.fa"
ref_base="${index_dir}/${ref_genome}" # 参照ゲノムのフルパス（ベース）

# --- 修正箇所 開始 ---
# 必要なインデックスファイルの確認
# GATK/Picard 用 + BWA 用
required_files=(
    "${ref_base}"            # FASTA本体
    "${ref_base}.fai"        # FASTA index (samtools faidx)
    "${ref_base%.fa}.dict"   # Sequence Dictionary (Picard CreateSequenceDictionary) ※拡張子.faを除いた名前に.dictをつける
    "${ref_base}.amb"        # BWA index
    "${ref_base}.ann"        # BWA index
    "${ref_base}.bwt"        # BWA index (重要)
    "${ref_base}.pac"        # BWA index
    "${ref_base}.sa"         # BWA index
)

echo "Checking for required reference and index files..."
all_files_found=true
missing_files=() # 見つからないファイルを格納する配列
for file in "${required_files[@]}"; do
    # -f は通常ファイルが存在するかを確認
    if [ ! -f "${file}" ]; then
        # echo "Error: Required file ${file} not found." # 個別のエラー表示は抑制
        all_files_found=false
        missing_files+=("${file}") # 見つからないファイルリストに追加
    fi
done

# すべてのファイルが見つからなかった場合に終了し、メッセージを送信
if [ "$all_files_found" = false ]; then
    echo "Error: One or more required reference/index files are missing."
    echo "Missing files:"
    # 見つからなかったファイル名をリスト表示
    for missing in "${missing_files[@]}"; do
        echo "  - ${missing}"
    done
    echo "Please ensure all files exist in the directory: ${index_dir}"
    echo "You might need to run:"
    echo "  bwa index ${ref_base}"
    echo "  samtools faidx ${ref_base}"
    echo "  picard CreateSequenceDictionary R=${ref_base} O=${ref_base%.fa}.dict" # Picardコマンドの例
    if [ -f ~/automessage.sh ]; then
        # エラーメッセージを組み立てる
        error_msg="Error: Missing required reference/index files for ${ref_genome}. Missing: $(IFS=, ; echo "${missing_files[*]}")"
        source ~/automessage.sh "${error_msg}"
    fi
    exit 1
fi
echo "All required reference and index files found."
# --- 修正箇所 終了 ---


# モジュールのロード
modules=("fastqc" "bwa" "samtools" "picard" "gatk")
for m in "${modules[@]}"; do
    echo "Loading module: $m" # ロード中のモジュールを表示（デバッグ用）
    module load $m
    if [ $? -ne 0 ]; then
        echo "Error: Failed to load $m"
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Failed to load $m"
        fi
        exit 1
    fi
done
echo "All required modules loaded successfully." # モジュールロード完了メッセージ

# 必要な入力ファイルの確認
input_files=("${s}_1.fastq" "${s}_2.fastq")
for file in "${input_files[@]}"; do
    if [ ! -s "$file" ]; then # -s オプションでファイルが存在し、かつ空でないことを確認
        echo "Error: Input file $file is missing or empty."
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Input file $file is missing or empty."
        fi
        exit 1
    fi
done
echo "Input FASTQ files found: ${input_files[*]}" # 入力ファイル確認完了メッセージ

# --- ここから先の処理は変更なし ---

# BWAによるマッピング
echo "Starting BWA mapping for ${s}..."
bwa mem -t 16 -R "@RG\tID:${s}\tSM:Dutch" ${ref_base} ${s}_1.fastq ${s}_2.fastq > ${s}.sam
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
    # エラー発生時は不要な中間ファイルを削除する方が良い場合がある
    # rm -f ${s}.sam
    exit 1
fi
echo "SAM to BAM conversion completed."
rm -f ${s}.sam # 変換成功後、不要なSAMファイルを削除

# BAMのソート
echo "Sorting BAM file for ${s}..."
samtools sort -@ 8 -T ${s}.sort_temp -o ${s}_sort.bam ${s}.bam # -@ でスレッド数指定、-T で一時ファイルプレフィックス指定
if [ $? -ne 0 ]; then
    echo "Error: BAM sorting failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: BAM sorting failed for ${s}"
    fi
    # rm -f ${s}.bam
    exit 1
fi
echo "BAM sorting completed."
rm -f ${s}.bam # ソート成功後、不要なソート前BAMファイルを削除

# BAMのインデックス作成
echo "Indexing sorted BAM file for ${s}..."
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
picard MarkDuplicates -I ${s}_sort.bam -M ${s}_marked_dup_metrics.txt -O ${s}_md.bam --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --CREATE_INDEX true # --CREATE_INDEX true で出力BAMのインデックスも同時に作成
if [ $? -ne 0 ]; then
    echo "Error: MarkDuplicates failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: MarkDuplicates failed for ${s}"
    fi
    exit 1
fi
echo "Marking duplicates completed."
# Picardがインデックスを作成するので、samtools index は不要（ただし明示的に実行しても問題ない）
# rm -f ${s}_sort.bam ${s}_sort.bam.bai # MarkDuplicates成功後、不要な入力ファイルを削除

# 重複の除去 (MarkDuplicates 後の処理としては一般的ではない場合もある。重複リードを除外したい場合)
# HaplotypeCaller は重複マークされたリードを内部で扱うため、このステップは必須ではないことが多い
# echo "Removing duplicates for ${s}..."
# samtools view -bh -F 1024 -F 4 ${s}_md.bam > ${s}_rmdup.bam # -F 1024 (PCR or optical duplicate), -F 4 (unmapped) を除外
# if [ $? -ne 0 ]; then
#     echo "Error: Removing duplicates failed for ${s}"
#     if [ -f ~/automessage.sh ]; then
#         source ~/automessage.sh "Error: Removing duplicates failed for ${s}"
#     fi
#     exit 1
# fi
# echo "Removing duplicates completed."
# rm -f ${s}_md.bam ${s}_md.bam.bai # rmdup成功後、不要な入力ファイルを削除

# 重複除去後のBAMのインデックス作成 (上記ステップを実行する場合に必要)
# echo "Indexing rmdup BAM file for ${s}..."
# samtools index ${s}_rmdup.bam
# if [ $? -ne 0 ]; then
#     echo "Error: BAM indexing (rmdup) failed for ${s}"
#     if [ -f ~/automessage.sh ]; then
#         source ~/automessage.sh "Error: BAM indexing (rmdup) failed for ${s}"
#     fi
#     exit 1
# fi
# echo "BAM indexing (rmdup) completed."

# GATK HaplotypeCaller (重複マークされたBAMを入力とする)
input_bam_for_hc="${s}_md.bam" # MarkDuplicatesの出力を使う
# input_bam_for_hc="${s}_rmdup.bam" # 重複除去した場合の入力
echo "Running GATK HaplotypeCaller for ${s} on ${input_bam_for_hc}..."
gatk --java-options "-Xmx32g" HaplotypeCaller -R ${ref_base} -I ${input_bam_for_hc} -O ${s}.vcf.gz --native-pair-hmm-threads 8 # 出力は .vcf.gz を推奨, HMM計算のスレッド指定
if [ $? -ne 0 ]; then
    echo "Error: HaplotypeCaller failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: HaplotypeCaller failed for ${s}"
    fi
    exit 1
fi
echo "GATK HaplotypeCaller completed."

# バリアントの選択（SNPとINDEL）
echo "Selecting SNPs for ${s}..."
gatk --java-options "-Xmx4g" SelectVariants -R ${ref_base} -V ${s}.vcf.gz -select-type SNP -O ${s}_snp.vcf.gz # gVCFではないのでメモリは少なめで良いかも
if [ $? -ne 0 ]; then
    echo "Error: SelectVariants (SNP) failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: SelectVariants (SNP) failed for ${s}"
    fi
    exit 1
fi
echo "SNP selection completed."

echo "Selecting INDELs for ${s}..."
gatk --java-options "-Xmx4g" SelectVariants -R ${ref_base} -V ${s}.vcf.gz -select-type INDEL -O ${s}_indel.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error: SelectVariants (INDEL) failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: SelectVariants (INDEL) failed for ${s}"
    fi
    exit 1
fi
echo "INDEL selection completed."

# バリアントフィルタリング（SNP） - GATKベストプラクティス推奨値を参考に
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
echo "Filtering SNPs for ${s}..."
gatk --java-options "-Xmx4g" VariantFiltration \
    -R ${ref_base} \
    -V ${s}_snp.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD_lt_2" \
    -filter "FS > 60.0" --filter-name "FS_gt_60" \
    -filter "MQ < 40.0" --filter-name "MQ_lt_40" \
    -filter "SOR > 3.0" --filter-name "SOR_gt_3" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum_lt_m12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_lt_m8" \
    -O ${s}_snp_filtered.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error: VariantFiltration (SNP) failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: VariantFiltration (SNP) failed for ${s}"
    fi
    exit 1
fi
echo "SNP filtering completed."

# バリアントフィルタリング（INDEL） - GATKベストプラクティス推奨値を参考に
echo "Filtering INDELs for ${s}..."
gatk --java-options "-Xmx4g" VariantFiltration \
    -R ${ref_base} \
    -V ${s}_indel.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD_lt_2" \
    -filter "FS > 200.0" --filter-name "FS_gt_200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum_lt_m20" \
    -filter "SOR > 10.0" --filter-name "SOR_gt_10" \
    -O ${s}_indel_filtered.vcf.gz
if [ $? -ne 0 ]; then
    echo "Error: VariantFiltration (INDEL) failed for ${s}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: VariantFiltration (INDEL) failed for ${s}"
    fi
    exit 1
fi
echo "INDEL filtering completed."

# 中間VCFファイルの削除（オプション）
# rm -f ${s}.vcf.gz ${s}.vcf.gz.tbi ${s}_snp.vcf.gz ${s}_snp.vcf.gz.tbi ${s}_indel.vcf.gz ${s}_indel.vcf.gz.tbi

echo "Pipeline completed successfully for sample: ${s}"
# スクリプト正常終了の通知
if [ -f ~/automessage.sh ]; then
    source ~/automessage.sh "Pipeline completed successfully for sample: ${s}"
fi

exit 0 # 正常終了