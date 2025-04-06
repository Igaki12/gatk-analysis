Gemini 2.5 pro　による解析

Methilome解析者です。解析中にエラーが発生しました。
bwaのこのエラーについて、解説と改善案を示して。

[igaki@emb02 Dutch]$ vi nohup1290783.sh 
[igaki@emb02 Dutch]$ s="SRR1290783"
[igaki@emb02 Dutch]$ index_dir="/home/igaki/data/rabbit2025/Dutch"
[igaki@emb02 Dutch]$ 
[igaki@emb02 Dutch]$ # ?~E?~A?~A??~B??~C??~C~G?~C~C?~B??~B??~C~U?~B??~B??~C??~A?確?~M
[igaki@emb02 Dutch]$ index_files=("oryCun2.fa" "oryCun2.fa.fai" "oryCun2.dict")
[igaki@emb02 Dutch]$ for file in "${index_files[@]}"; do
>     if [ ! -f "${index_dir}/${file}" ]; then
>         echo "Error: Required index file ${file} not found in ${index_dir}."
>         if [ -f ~/automessage.sh ]; then
>             source ~/automessage.sh "Error: Required index file ${file} not found in ${index_dir}."
>         fi
>         exit 1
>     fi
> done
[igaki@emb02 Dutch]$ 
[igaki@emb02 Dutch]$ vi nohup1290783.sh 
[igaki@emb02 Dutch]$ bwa mem -t 16 -R "@RG\tID:${s}\tSM:Dutch" ${index_dir}/oryCun2.fa ${s}_1.fastq ${s}_2.fastq > ${s}.sam
bash: bwa: command not found...
Install package 'bwa' to provide command 'bwa'? [N/y] N 

[igaki@emb02 Dutch]$ module load bwa
[igaki@emb02 Dutch]$ bwa mem -t 16 -R "@RG\tID:${s}\tSM:Dutch" ${index_dir}/oryCun2.fa ${s}_1.fastq ${s}_2.fastq > ${s}.sam
[E::bwa_idx_load_from_disk] fail to locate the index files