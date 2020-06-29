#The fly, zebrafish, human and mouse transposon sequence were downloaded from UCSC(http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=759454897_6o3yAdmKXNcJR2vLsHwFC6b1ATAo&clade=virus&org=0&db=0&hgta_group=varRep&hgta_track=microsat&hgta_table=microsat&hgta_regionType=genome&position=&hgta_outputType=sequence&hgta_outFileName=mouse_mm10_ucsc)
#RepeatMasker track was selected and named by species_version_ucsc.fa (fly_dm6_ucsc.fa, zebrafish_danRer11_ucsc.fa, human_hg38_ucsc.fa, mouse_mm10_ucsc.fa).
#Gene CDS sequence were downloaded from Ensemble release 97.
#For transposons most species do not have proper CDS sequence, so we use TransDecoder(https://github.com/TransDecoder/TransDecoder/wiki) to identify candidate coding regions.

TransDecoder.LongOrfs -S -t fly_dm6_ucsc.fa
blastp -query fly_dm6_ucsc.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > fly_dm6_ucsc.fa_longest_orfs_blastp.outfmt
hmmscan --cpu 8 --domtblout fly_dm6_ucsc.fa.pfam.domtblout ~/Pfam-A.hmm fly_dm6_ucsc.fa.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t fly_dm6_ucsc.fa --retain_pfam_hits fly_dm6_ucsc.fa.pfam.domtblout --retain_blastp_hits fly_dm6_ucsc.fa_longest_orfs_blastp.outfmt

TransDecoder.LongOrfs -S -t zebrafish_danRer11_ucsc.fa
blastp -query zebrafish_danRer11_ucsc.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > zebrafish_danRer11_ucsc.fa_longest_orfs_blastp.outfmt
hmmscan --cpu 8 --domtblout zebrafish_danRer11_ucsc.fa.pfam.domtblout ~/Pfam-A.hmm zebrafish_danRer11_ucsc.fa.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t fly_dm6_ucsc.fa --retain_pfam_hits zebrafish_danRer11_ucsc.fa.pfam.domtblout --retain_blastp_hits zebrafish_danRer11_ucsc.fa_longest_orfs_blastp.outfmt

TransDecoder.LongOrfs -S -t human_hg38_ucsc.fa
blastp -query human/human_hg38_ucsc.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > human_hg38_ucsc.fa_longest_orfs_blastp.outfmt
hmmscan --cpu 8 --domtblout human_hg38_ucsc.fa.pfam.domtblout ~/Pfam-A.hmm human/human_hg38_ucsc.fa.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t human_hg38_ucsc.fa --retain_pfam_hits human_hg38_ucsc.fa.pfam.domtblout --retain_blastp_hits human_hg38_ucsc.fa_longest_orfs_blastp.outfmt

TransDecoder.LongOrfs -S -t mouse_mm10_ucsc.fa
blastp -query mouse/mouse_mm10_ucsc.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > mouse_mm10_ucsc.fa_longest_orfs_blastp.outfmt
hmmscan --cpu 8 --domtblout mouse_mm10_ucsc.fa.pfam.domtblout ~/Pfam-A.hmm mouse/mouse_mm10_ucsc.fa.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t mouse_mm10_ucsc.fa --retain_pfam_hits mouse_mm10_ucsc.fa.pfam.domtblout --retain_blastp_hits mouse_mm10_ucsc.fa_longest_orfs_blastp.outfmt

#Both for genes and transposons, we selected the longest CDS for each.
#fly gene
vi Drosophila_melanogaster.BDGP6.22.cds.all.fa  %s# #_# three times
seqkit fx2tab -j 30 -l  -n -i -H Drosophila_melanogaster.BDGP6.22.cds.all.fa | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' > fly_gene.cds_len
vi fly_gene.cds_len  1dd  :%s#^\S*_gene:## :%s#_cds\S*##
sort -k1,1 -k3,3nr fly_gene.cds_len | awk '!a[$1]++{print}' > fly_gene.cds_longest_len
#fly TE
seqkit fx2tab -j 30 -l  -n -i -H fly_dm6_ucsc.fa.transdecoder.cds | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' > fly_transdecoder.cds_len
vi fly_transdecoder.cds_len :1dd :%s#\.p\d\+##
sort -k1,1 -k3,3nr fly_transdecoder.cds_len |awk '!a[$1]++{print}' > fly_transdecoder.longest.cds_len

#zebrafish gene
vi Danio_rerio.GRCz11.cds.all.fa  %s# #_# three times
seqkit fx2tab -j 30 -l  -n -i -H Danio_rerio.GRCz11.cds.all.fa | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' |sed '1d' |sed 's#^\S*_gene:##'| sed 's#_cds\S*##' | sort -k1,1 -k3,3nr | awk '!a[$1]++{print}' > zebrafish_gene.cds_longest_len
#zebrafish TE
seqkit fx2tab -j 30 -l  -n -i -H zebrafish_danRer11_ucsc.fa.transdecoder.cds | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' > zebrafish_transdecoder.cds_len
vi zebrafish_transdecoder.cds_len  :1dd :%s#\.p\d\+##
sort -k1,1 -k3,3nr zebrafish_transdecoder.cds_len |awk '!a[$1]++{print}' > zebrafish_transdecoder.longest.cds_len


#human gene
vi Homo_sapiens.GRCh38.cds.all.fa  %s# #_# three times
seqkit fx2tab -j 30 -l  -n -i -H Homo_sapiens.GRCh38.cds.all.fa | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' |sed '1d' |sed 's#^\S*_gene:##'| sed 's#_cds\S*##' | sort -k1,1 -k3,3nr | awk '!a[$1]++{print}' > human_hg38_gene.cds_longest_len
#human TE
seqkit fx2tab -j 30 -l  -n -i -H human_hg38_ucsc.fa.transdecoder.cds | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' > human_transdecoder.cds_len
vi human_transdecoder.cds_len :1dd :%s#\.p\d\+##
sort -k1,1 -k3,3nr human_transdecoder.cds_len |awk '!a[$1]++{print}' > human_transdecoder.longest.cds_len
#mouse gene
vi Mus_musculus.GRCm38.cds.all.fa  %s# #_# three times
seqkit fx2tab -j 30 -l  -n -i -H Mus_musculus.GRCm38.cds.all.fa | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' |sed '1d' |sed 's#^\S*_gene:##'| sed 's#_cds\S*##' | sort -k1,1 -k3,3nr | awk '!a[$1]++{print}' > mouse_mm10_gene.cds_longest_len
#mouse TE
seqkit fx2tab -j 30 -l  -n -i -H mouse_mm10_ucsc.fa.transdecoder.cds | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' > mouse_transdecoder.cds_len
vi mouse_transdecoder.cds_len :1dd :%s#\.p\d\+##
sort -k1,1 -k3,3nr mouse_transdecoder.cds_len |awk '!a[$1]++{print}' > mouse_transdecoder.longest.cds_len
