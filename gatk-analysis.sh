# [igaki@emb02 Dutch]$ bwa mem -t 16 -R "@RG\tID:${f}\tSM:Dutch" ~/rabbit2025/index_bwa/oryCun2.fa ~/data/rabbit2025/Dutch/SRR1290783_1.fastq ~/data/rabbit2025/Dutch/SRR1290783_2.fastq > ${f}.sam; echo $f;
# > samtools view -bS ${s}.sam > ${s}.bam;
# [igaki@emb02 bwa]$ for s in $sample ; do samtools sort -T ${s}.sort -o ${s}_sort.bam ~/rabbit2025/bwa/${s}.bam; samtools index ${s}_sort.bam; echo $s; done;source ~/autoreminder.sh
# for s in $sample ; do picard MarkDuplicates I=../2_samtools/${s}_sort.bam M=${s}_marked_dup_metrics.txt O=${s}_md.bam ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT; echo ${s}; done;source ~/autoreminder.sh

# [igaki@emb02 3_picard_retry3--SORTING_COLLECTION_SIZE_RATIO]$ for s in $sample ; do picard MarkDuplicates -I ../2_samtools/${s}_sort.bam -M ${s}_marked_dup_metrics.txt -O ${s}_md.bam --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --SORTING_COLLECTION_SIZE_RATIO 0.1 ; echo ${s}; done;source ~/autoreminder.sh
# > samtools view -bh -F 1036 ${f}_md.bam > ${f}_rmdup.bam;
# > samtools index ${f}_rmdup.bam
# gatk --java-options "-Xmx4g" HaplotypeCaller -R ~/data/rabbit2025/oryCun2.fa -I ~/rabbit2025/3_picard_retry3--SORTING_COLLECTION_SIZE_RATIO/${f}_rmdup.bam -O ${f}.vcf;