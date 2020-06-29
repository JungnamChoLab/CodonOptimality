mkdir 0_fastqc 3_trim 4_hisat2 5_stringtie 6_featurecounts
cd 2.cleandata;ls | sed /name/d > name
cp name ../3_trim;cp name ../4_hisat2;cp name ../5_stringtie;
#md5check;mv data;
for i in `cat name`;do cd $i; md5sum -c MD5_${i}.txt >> ../../alllog 2>&1;mv *.fq.gz ..;cd ..;done
#fastqc
~/softwares/FastQC/fastqc -t 10 -o ../0_fastqc *.fq.gz;cd ../3_trim;
#trim
for i in `cat name`;do echo $i;java -jar ~/softwares/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 10 -phred33 ../2.cleandata/${i}_1.clean.fq.gz ../2.cleandata/${i}_2.clean.fq.gz -baseout ${i}.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:1 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 >> ../alllog 2>&1;done
#mapping
conda activate python27;cd ../4_hisat2;
for i in `cat name`;do echo $i;hisat2 -p 10 --dta --rg-id $i --rg SM:$i --rna-strandness RF --fr -x ~/genome/tair10seq -1 ../3_trim/${i}_1P.fq.gz -2 ../3_trim/${i}_2P.fq.gz -S ${i}.sam >> ../alllog 2>&1;samtools view -bS ${i}.sam > ${i}.bam;samtools sort -@ 10 ${i}.bam ${i}.sorted;samtools index ${i}.sorted.bam;rm ${i}.sam;rm ${i}.bam;done
conda activate base
cd ../5_stringtie;
#count,methods1:stringtie
for i in `cat name`;do echo $i;stringtie ../4_hisat2/${i}.sorted.bam  -e -B -A ${i}.gene.abundance -G ~/genome/TAIR10_GFF3_genes_transposons.gff -p 10 -o ${i}/{$i}.gtf;sort ${i}.gene.abundance -o ${i}.gene.abundance;done
#count,methods2:featurecounts
cd ../6_featurecounts
featureCounts -p -M -T 10 -a ~/genome/TAIR10_GFF3_genes_exons.gtf -o featurecounts ../4_hisat2/*.sorted.bam
