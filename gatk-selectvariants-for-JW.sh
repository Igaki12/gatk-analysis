#!/bin/bash

# --- 設定 ---
# VCFファイルが格納されているディレクトリ
VCF_DIR=~/rabbit2025/4_gatk
# 参照ゲノムファイルのフルパス (!!! 実際のパスに変更してください !!!)
REF_GENOME=~/data/rabbit2025/Dutch/oryCun2.fa
# GATK実行時のJavaメモリ設定
GATK_JAVA_OPTS="-Xmx32g"
# --- 設定ここまで ---

# --- 事前チェック ---
# GATKモジュールのロード
echo "Loading GATK module..."
module load gatk
if [ $? -ne 0 ]; then
    echo "Error: Failed to load GATK module."
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: Failed to load GATK module."
    fi
    exit 1
fi
echo "GATK module loaded."

# 参照ゲノム関連ファイルの存在確認
if [ ! -f "${REF_GENOME}" ]; then
    echo "Error: Reference genome file not found: ${REF_GENOME}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: Reference genome file not found: ${REF_GENOME}"
    fi
    exit 1
fi
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Error: Reference genome index file (.fai) not found for ${REF_GENOME}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: Reference genome index file (.fai) not found for ${REF_GENOME}"
    fi
    exit 1
fi
# GATKは .dict ファイルも参照することが多い
REF_DICT=$(dirname ${REF_GENOME})/$(basename ${REF_GENOME} .fa).dict
if [ ! -f "${REF_DICT}" ]; then
    # .dict ファイルがない場合は、faファイルと同じディレクトリにあるか確認
    REF_DICT_ALT="${REF_GENOME%.*}.dict" # .fa を .dict に置換
    if [ ! -f "${REF_DICT_ALT}" ]; then
        echo "Warning: Reference genome dictionary file (.dict) not found at ${REF_DICT} or ${REF_DICT_ALT}"
        echo "GATK might require this file. Consider creating it using 'gatk CreateSequenceDictionary -R ${REF_GENOME}'."
    fi
fi

# 入力ディレクトリの存在確認
if [ ! -d "${VCF_DIR}" ]; then
    echo "Error: Input VCF directory not found: ${VCF_DIR}"
    if [ -f ~/automessage.sh ]; then
        source ~/automessage.sh "Error: Input VCF directory not found: ${VCF_DIR}"
    fi
    exit 1
fi

echo "Starting VCF selection and filtering in directory: ${VCF_DIR}"
echo "Using reference genome: ${REF_GENOME}"

# 各bgzip圧縮済みVCFファイルについて処理
processed_count=0
error_count=0
for file in "$VCF_DIR"/*.vcf.gz; do
    # ファイルが存在しない場合 (globが何もマッチしなかった場合) スキップ
    if [ ! -e "$file" ]; then
        echo "Warning: No .vcf.gz files found matching pattern in ${VCF_DIR}"
        break # ループを抜ける
    fi
    # 既に処理済みのファイル（_snp_filtered.vcf.gz など）はスキップ
    if [[ "$file" == *"_snp.vcf.gz" || "$file" == *"_indel.vcf.gz" || "$file" == *"_filtered.vcf.gz" ]]; then
        echo "Skipping already processed or intermediate file: $file"
        continue
    fi

    # サンプル名（ファイル名から拡張子除く）を取得
    # 注意: 入力ファイル名がサンプル名でない可能性もあるため、ここではファイルベース名を使用
    base_name=$(basename "$file" .vcf.gz)
    echo "--------------------------------------------------"
    echo "Processing base file: ${base_name} (Input VCF: ${file})"
    echo "--------------------------------------------------"

    # 入力VCFに対応するインデックス(.tbi)の存在確認
    if [ ! -f "${file}.tbi" ]; then
        echo "Error: Index file (.tbi) not found for ${file}. Skipping this file."
        echo "Please index the VCF file using 'gatk IndexFeatureFile -I ${file}' or 'tabix -p vcf ${file}'."
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: Index file (.tbi) not found for ${file}. Skipping ${base_name}."
        fi
        error_count=$((error_count + 1))
        continue # 次のファイルへ
    fi

    # 出力ファイル名の定義 (入力ファイルと同じディレクトリに出力)
    SNP_VCF="${VCF_DIR}/${base_name}_snp.vcf.gz"
    INDEL_VCF="${VCF_DIR}/${base_name}_indel.vcf.gz"
    SNP_FILTERED_VCF="${VCF_DIR}/${base_name}_snp_filtered.vcf.gz"
    INDEL_FILTERED_VCF="${VCF_DIR}/${base_name}_indel_filtered.vcf.gz"

    # --- バリアントの選択（SNP） ---
    echo "[${base_name}] Selecting SNPs..."
    gatk --java-options "${GATK_JAVA_OPTS}" SelectVariants \
        -R "${REF_GENOME}" \
        -V "${file}" \
        --select-type-to-include SNP \
        -O "${SNP_VCF}"
    if [ $? -ne 0 ]; then
        echo "Error: SelectVariants (SNP) failed for ${base_name}"
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: SelectVariants (SNP) failed for ${base_name}"
        fi
        error_count=$((error_count + 1))
        echo "Skipping further processing for ${base_name} due to SNP selection error."
        continue # 次のファイルへ
    fi
    echo "[${base_name}] SNP selection completed."

    # --- バリアントの選択（INDEL） ---
    echo "[${base_name}] Selecting INDELs..."
    gatk --java-options "${GATK_JAVA_OPTS}" SelectVariants \
        -R "${REF_GENOME}" \
        -V "${file}" \
        --select-type-to-include INDEL \
        -O "${INDEL_VCF}"
    if [ $? -ne 0 ]; then
        echo "Error: SelectVariants (INDEL) failed for ${base_name}"
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: SelectVariants (INDEL) failed for ${base_name}"
        fi
        # 中間SNPファイルを削除しておく
        rm -f "${SNP_VCF}" "${SNP_VCF}.tbi"
        error_count=$((error_count + 1))
        echo "Skipping further processing for ${base_name} due to INDEL selection error."
        continue # 次のファイルへ
    fi
    echo "[${base_name}] INDEL selection completed."

    # --- バリアントフィルタリング（SNP） ---
    echo "[${base_name}] Filtering SNPs..."
    gatk --java-options "${GATK_JAVA_OPTS}" VariantFiltration \
        -R "${REF_GENOME}" \
        -V "${SNP_VCF}" \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" \
        --filter-name "SNP_HardFilter" \
        -O "${SNP_FILTERED_VCF}"
    if [ $? -ne 0 ]; then
        echo "Error: VariantFiltration (SNP) failed for ${base_name}"
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: VariantFiltration (SNP) failed for ${base_name}"
        fi
        # 中間ファイルを削除しておく
        rm -f "${SNP_VCF}" "${SNP_VCF}.tbi"
        rm -f "${INDEL_VCF}" "${INDEL_VCF}.tbi"
        error_count=$((error_count + 1))
        echo "Skipping INDEL filtering for ${base_name} due to SNP filtering error."
        continue # 次のファイルへ
    fi
    echo "[${base_name}] SNP filtering completed."
    # フィルタ前のSNPファイルを削除
    echo "[${base_name}] Removing intermediate SNP VCF: ${SNP_VCF}"
    rm -f "${SNP_VCF}" "${SNP_VCF}.tbi"

    # --- バリアントフィルタリング（INDEL） ---
    echo "[${base_name}] Filtering INDELs..."
    gatk --java-options "${GATK_JAVA_OPTS}" VariantFiltration \
        -R "${REF_GENOME}" \
        -V "${INDEL_VCF}" \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
        --filter-name "INDEL_HardFilter" \
        -O "${INDEL_FILTERED_VCF}"
    if [ $? -ne 0 ]; then
        echo "Error: VariantFiltration (INDEL) failed for ${base_name}"
        if [ -f ~/automessage.sh ]; then
            source ~/automessage.sh "Error: VariantFiltration (INDEL) failed for ${base_name}"
        fi
        # 中間ファイルを削除しておく
        rm -f "${INDEL_VCF}" "${INDEL_VCF}.tbi"
        error_count=$((error_count + 1))
        continue # 次のファイルへ
    fi
    echo "[${base_name}] INDEL filtering completed."
    # フィルタ前のINDELファイルを削除
    echo "[${base_name}] Removing intermediate INDEL VCF: ${INDEL_VCF}"
    rm -f "${INDEL_VCF}" "${INDEL_VCF}.tbi"

    echo "[${base_name}] Successfully processed."
    processed_count=$((processed_count + 1))

done

echo "=================================================="
echo "All VCF processing finished."
echo "Successfully processed files: ${processed_count}"
echo "Files skipped due to errors: ${error_count}"
echo "Filtered files generated in: ${VCF_DIR}"
echo "=================================================="

# 完了メッセージ
if [ -f ~/automessage.sh ]; then
    if [ ${error_count} -eq 0 ]; then
        source ~/automessage.sh "VCF selection and filtering completed successfully for ${processed_count} files in ${VCF_DIR}."
    else
        source ~/automessage.sh "VCF selection and filtering finished with ${error_count} errors. ${processed_count} files processed successfully in ${VCF_DIR}."
    fi
fi

# エラーがあった場合は終了ステータスを1に、なければ0にする
if [ ${error_count} -gt 0 ]; then
    exit 1
else
    exit 0
fi