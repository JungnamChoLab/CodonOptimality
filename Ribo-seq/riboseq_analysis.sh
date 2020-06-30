# Ribo-seq analysis

#For the raw reads, cutadpt (version 1.12) was used to trim adapters; fastq_quality_filter (fastx_toolkit-0.0.14) discards the poor-quality reads (>50% bases with a Phred score <20);
#bowtie (version 1.0.1) (-l 20 -v 0) discard mycoplasma reads	and rRNA Reads. we use the clean reads from the company directly for the following analysis.
for i in `cat name`;do tophat2 -p 15 -o $i -G ~/genome/TAIR10_GFF3_genes_exons.gtf ~/genome/tair10seq ../clean_data/${i}.fq.gz>>mappinglog 2>&1;done
for i in `cat name`;do cufflinks -p 15 -G TAIR10_GFF3_genes_exons.gtf -o $i -u --library-type fr-firststrand ../bam/${i}.bam;done

# RNA-seq analysis
# The raw data were firstly processed through in-house perl scripts to remove reads containing adapter, ploy-N and low quality by the company. Clean reads



# Disome analysis of ddm1
for i in `cat name`;do echo $i;java -jar ~/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 -phred33 ../raw_data/${i}_1.fq.gz ../ra
w_data/${i}_2.fq.gz -baseout ${i}.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 >> ../all_log 2>&
1;done
rm *_2*

cat RISDR1_1P.fq RISDR1_1U.fq > RISDR1.fq #gzip RISDR1.fq
cat RISDR2_1P.fq RISDR2_1U.fq > RISDR2.fq #gzip RISDR2.fq

perl filter_fq_length.pl RISDR1.fq > RISDR1_40_65.fq
perl filter_fq_length.pl RISDR2.fq > RISDR2_40_65.fq

cat RISDR1_40_65.fq RISDR2_40_65.fq > ddm1_40_65.fq
tophat -p 8 -o ddm1_40_65 -G ~/genome/TAIR10_GFF3_genes_transposons.gff ~/genome/tair10seq ../raw_data_trimmed/ddm1_40_65.fq
cufflinks -u -g ~/genome/TAIR10_GFF3_genes_transposons.gff ../2_mapping/ddm1_40_65/accepted_hits.bam


# Visulalizing of Riboseq data
# We 
bowtie-build ~/genome/TAIR10_cdna.fa tair_trans
bowtie -p 12 -l 23 tair_trans ../clean_data/RISCR1_clean.fastq.gz -S RISCR1.sam
bowtie -p 12 -l 23 tair_trans ../clean_data/RISCR2_clean.fastq.gz -S RISCR2.sam
bowtie -p 12 -l 23 tair_trans ../clean_data/RISDR1_clean.fastq.gz -S RISDR1.sam
bowtie -p 12 -l 23 tair_trans ../clean_data/RISDR2_clean.fastq.gz -S RISDR2.sam
bowtie -p 12 -l 23 ../bam/tair_trans ../raw_data_trimmed/ddm1_40_65.fq -S ddm1_40_65_transcripts.sam
samtools view -bS -F 4 -@ 5 RISCR1.sam > RISCR1.bam
samtools view -bS -F 4 -@ 5 RISCR2.sam > RISCR2.bam
samtools view -bS -F 4 -@ 5 RISDR1.sam > RISDR1.bam
samtools view -bS -F 4 -@ 5 RISDR2.sam > RISDR2.bam
samtools view -bS -F 4 -@ 5 ddm1_40_65_transcripts.sam > ddm1_40_65_transcripts.bam
bowtie2-build ~/genome/TAIR10_cdna.fa tair_transord
#bowtie2 -p 12 -x tair_transord -1 ../rnaseq_clean_data/RSCR1_1.fq.gz -2 ../rnaseq_clean_data/RSCR1_2.fq.gz -S RSCR1_trans.sam;samtools view -bS -F4 -@ 10 RSCR1_trans.sam >RSCR1_trans.bam
bowtie2 -p 12 -x tair_transord -1 ../rnaseq_clean_data/RSCR2_1.fq.gz -2 ../rnaseq_clean_data/RSCR2_2.fq.gz -S RSCR2_trans.sam;samtools view -bS -F4 -@ 10 RSCR2_trans.sam >RSCR2_trans.bam
bowtie2 -p 12 -x tair_transord -1 ../rnaseq_clean_data/RSDR1_1.fq.gz -2 ../rnaseq_clean_data/RSDR1_2.fq.gz -S RSDR1_trans.sam;samtools view -bS -F4 -@ 10 RSDR1_trans.sam >RSDR1_trans.bam
bowtie2 -p 12 -x tair_transord -1 ../rnaseq_clean_data/RSDR2_1.fq.gz -2 ../rnaseq_clean_data/RSDR2_2.fq.gz -S RSDR2_trans.sam;samtools view -bS -F4 -@ 10 RSDR2_trans.sam >RSDR2_trans.bam
