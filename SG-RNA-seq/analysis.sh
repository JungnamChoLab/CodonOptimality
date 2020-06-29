#first in the clean data
ls | sed /name/d > name
cp name ../mapping;cp name ../stringtie;
#md5check
for i in `cat name`;do cd $i; md5sum -c MD5_${i}.txt 
#fastqc, we first check the data quality
~/softwares/FastQC/fastqc -t 10 -o ../0_fastqc *.fq.gz
#For the clean data got from the company, it has no adapters need to be trimmed and the data quality is very good, so we can use it directly for the following analysis. 
#mapping
for i in `cat name`;do echo $i;hisat2 -p 15 --dta --rg-id $i --rg SM:$i --rna-strandness RF --fr -x ~/genome/tair10seq -1 ../clean_data/${i}_1.clean.fq.gz -2 ../clean_data/${i}_2.clean.fq.gz -S ${i}.sam >> mapping-record 2>&1;samtools view -bSF 4 ${i}.sam > ${i}.us.bam;samtools sort -@ 15 -o ${i}.bam ${i}.us.bam;samtools index ${i}.bam;done
#use stringtie got the fpkm of genes and transposons
for i in `cat name`;do stringtie ../mapping/${i}.bam  -e -B -A ${i}.gene.abundance -G ~/genome/TAIR10_GFF3_genes_exons.gtf -p 15 -o $i/${i}.gtf --rf;done
#use featurecounts got the raw read counts of gene and transposon
cp ~/genome/TAIR10_GFF3_genes_exons.gtf TAIR10_GFF3_genes_exons_modified.gtf
vi TAIR10_GFF3_genes_exons_modified.gtf  ####  %s#pseudogenic_exon#exon#
featureCounts -s 2 -p -M -T 10 -a TAIR10_GFF3_genes_exons_modified.gtf -o stress_granule_raw_read_count ../mapping/G5_FRRB190324159-1a.bam ../mapping/G5_FRRB190324159-2a.bam ../mapping/G5T_FRRB190324160-1a.bam ../mapping/G5T_FRRB190324160-2a.bam
