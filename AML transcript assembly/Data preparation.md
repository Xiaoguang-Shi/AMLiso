# NGS raw data QC
```bash
trimmomatic PE -threads 10 -phred33 \
       $read1 $read2 \
       trimmomatic/${id}_1_trim.fq.gz trimmomatic/${id}_1_unpaired.fq.gz \
       trimmomatic/${id}_2_trim.fq.gz trimmomatic/${id}_2_unpaired.fq.gz \
       ILLUMINACLIP:trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50

fastp  -threads 10 \
    -i trimmomatic/${id}_1_trim.fq.gz  \
    -I trimmomatic/${id}_2_trim.fq.gz \
    -o trimmomatic/${id}_1_clean.fq.gz \
    -O trimmomatic/${id}_2_clean.fq.gz \
    --detect_adapter_for_pe -q 20 -l 50 \
    -j ${id}.fastp.json \
    -h ${id}.fastp.html
```
# ONT basecalling

```bash
guppy_basecaller \
  --input_path fast5_pass \
  --save_path basecall \
  --config dna_r9.4.1_450bps_sup.cfg \
  --device cuda:0
```
# ONT QC
```
NanoStat --fastq  ${id}_NeoA_ONT_TR.fastq.gz -n $id.stat.txt
```
# full legth reads slect
```
pychopper -t 40 -k PCS109 -m phmm ${id}_NeoA_ONT_TR.fastq.gz  ${id}_NeoA_ONT_TR_full_length_reads.fastq
```
