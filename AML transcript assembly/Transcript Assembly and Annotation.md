# Transcript Assembly and Annotation

## 1.StringTie3

### 1.1 ONT Reads Mapping
```bash
minimap2 -t 40 -ax splice --MD GRCh38.p14.genome.fa ${id}_NeoA_ONT_TR_full_length_reads.fastq.gz | \
samtools view -@ 40 -bS - | \
samtools sort -@ 40 -o ${id}.ont.sorted.bam -
```

### 1.2 NGS Reads Mapping
```bash
hisat2 -t -x ../reference/hisat2/GRCh38.p14_snp_tran_index \
-1 ${fastq}/${id}_1_clean.fq.gz \
-2 ${fastq}/${id}_2_clean.fq.gz \
-p 40 | \
samtools view -@ 40 -bS - | \
samtools sort -@ 40 -o ${id}.sorted.bam -
```

### 1.3 GERCG Transcript Assembly
```bash
gffcompare -C -o GERCG \
gencode.v47.chr_patch_hapl_scaff.annotation.gtf \
ensembl.gtf \
refseq.gtf \
chess.gtf \
gtex_flair_filter_transcripts.gtf
```

### 1.4 StringTie Transcript Assembly
```bash
mkdir ${id} && cd ${id}
stringtie -p 4 --mix -G GERCG.combined.gtf -o stringtie.out.gtf ../${id}.sorted.bam ../${id}.ont.sorted.bam
```

## 2. TALON

### 2.1 Index Build
```bash
talon_initialize_database \
--f GERCG.combined.gtf \
--a GERCG \
--g GRCh38 \
--l 0 \
--idprefix GERCG \
--5p 500 \
--3p 300 \
--o GERCG
```

### 2.2 Label Reads
```bash
talon_label_reads \
--f=${id}.ont.sam \
--g=GRCh38.p14.genome.fa \
--t=20 \
--ar=20 \
--deleteTmp \
--o=${id}
```

### 2.3 TALON Running
```bash
echo "${id},${id},PromethION24,${id}_labeled.sam" >> configure.csv

talon \
--f configure.csv \
--db GERCG.db \
--threads 40 \
--build GRCh38 \
--nsg \
--o ${id}
```

### 2.4 Filter
```bash
talon_filter_transcripts \
--db GERCG.db \
--datasets ${id} \
-a GERCG \
--o ${id}_filtered_transcripts.csv
```

### 2.5 Create GTF
```bash
talon_create_GTF \
--db GERCG.db \
--whitelist ${id}_filtered_transcripts.csv \
-a GERCG \
--build GRCh38 \
--o ${id}
```

## 3. FLAIR

### 3.1 Align
```bash
source activate flair

flair align -g GRCh38.p14.genome.fa \
-r ${id}_NeoA_ONT_TR_full_length_reads.fastq.gz \
--threads 40 \
--output ${id}
```

### 3.2 Junctions Processing
```bash
junctions_from_sam -s ${id}.sorted.bam --unique
sed -i 's/chrG/G/g; s/chrJ/J/g; s/chrK/K/g; s/chrML/ML/g; s/chrMU/MU/g' junctions_from_sam_junctions.bed
```

### 3.3 Correct
```bash
flair correct -q ${id}.bed \
--threads 40 \
-f ${gtf} \
-j junctions_from_sam_junctions.bed \
-g ${genome} \
-o ${id}.correct
```

### 3.4 Collapse
```bash
flair collapse -g ${genome} \
-q ${id}.bed \
-r ${id}_NeoA_ONT_TR_full_length_reads.fastq.gz \
--threads 40 \
--gtf ${gtf} \
--temp_dir ./ \
--output ${id} \
--generate_map \
--annotation_reliant generate \
--support 10 \
--stringent \
--check_splice
```

## 4. Transcript Merge
```bash
gffcompare -C -o AML_raw -i gtf.list -p AML_tmp
```

## 5. Transcript Filter

### 5.1 SQANTI3 Quality Control
```bash
python sqanti3_qc.py \
AML_raw.combined.gtf \
GERCG.combined.gtf GRCh38.p14.genome.fa \
--skipORF \
--expression ${id}.tsv \
-c STAR/${id}SJ.out.tab \
--CAGE_peak ${cage_peak} \
-o ${id} -d ${id} \
--cpus 40 --report skip
```

### 5.2 Filter Configuration (filter_min_cov10.json)
```json
{
    "full-splice_match": [{
        "perc_A_downstream_TTS": [0,59]
    }],
    "incomplete-splice_match": [{
        "perc_A_downstream_TTS": [0,59]
    }],
    "novel_in_catalog": [{
        "perc_A_downstream_TTS": [0,59],
        "min_cov": 10
    }],
    "novel_not_in_catalog": [{
        "perc_A_downstream_TTS": [0,59],
        "min_cov": 10
    }],
    "rest": [{
        "perc_A_downstream_TTS": [0,59]
    }]
}
```

### 5.3 Run Filter
```bash
sqanti3_filter.py rules \
-j filter_min_cov10.json \
-o cov10 ${id}_classification.txt
```

### 5.4 PolyA QC
```bash
# 生成polyA reads
python extract_polyA_reads.py ${short}/${id}_1_clean.fq.gz ${id}_nsg_1_polyA.fq.gz
python extract_polyA_reads.py ${short}/${id}_2_clean.fq.gz ${id}_nsg_2_polyA.fq.gz
python extract_polyA_reads.py ${ont}/${id}_NeoA_ONT_TR_full_length_reads.fastq.gz ${id}_ont_polyA.fq.gz

# 比对和过滤
hisat2 -t -x genome -1 ${id}_nsg_1_polyA.fq.gz -2 ${id}_nsg_2_polyA.fq.gz -S ${id}.ngs.polyA.sam -p 10
samtools view -bS ${id}.ngs.polyA.sam | samtools sort -o - | samtools view -h -q 1 -F 4 -F 256 -F 1024 - | grep -v XA:Z | grep -v SA:Z | samtools view -b - > ${id}.ngs.polyAdedup.unique.bam

# 组装polyA
stringtie -p 4 -G GERCG.combined.gtf -o polyA.out.gtf ${id}.ngs.polyAdedup.unique.bam
```

## 6. ORF Finding
```bash
orfanage --reference GRCh38.p14.genome.fa \
--query AML_filtered.gtf \
--threads 10 \
--output AML_filtered.CDS.gtf \
gencode.v47.chr_patch_hapl_scaff.annotation.gtf \
ensembl.gtf \
refseq.gtf \
chess.gtf \
gtex_flair_filter_transcripts.gtf
```

## 7. Annotation

### 7.1 EggNOG Annotation
`http://eggnog-mapper.embl.de`

### 7.2 InterProScan
```bash
interproscan-5.65-97.0/interproscan.sh \
-i aml.faa \
-cpu 40 \
-f tsv \
-dp \
-goterms \
-pa \
--tempdir $PWD
```
