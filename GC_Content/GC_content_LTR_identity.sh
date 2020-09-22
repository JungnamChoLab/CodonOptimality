###Copia (RLC) and Gypsy (RLG) Class I Retroelements of O. sativa ssp. japonica were downloaded from http://www.genome.arizona.edu/rite/. And the download link is http://www.genome.arizona.edu/cgi-bin/rite/download.cgi?id=618&key=953f7b45c68adc50e2663f4d31620e08.
###The downloaded file was named "rice_ltc_ltg.fa". Total have 5160 sequences. 

cp rice_ltc_ltg.fa rice_ltc_ltg_.fa   
vi rice_ltc_ltg_.fa                                            ### The fasta ids were shorted.       %s#RLC_\#LTR/## %s#RLG_\#LTR/## 
perl ../split_fa.pl ../rice_ltc_ltg_.fa                        ### The fasta sequences were renamed and splited into individual sequences.
### We use blast to get the identity of left and right LTR.
for i in {1..5160};do makeblastdb -in $i -dbtype nucl -parse_seqids -out $i; done
for i in {1..5160};do blastn -query $i -out ${i}.blast -db $i -outfmt 6;done
for i in {1..5160};do cat ${i}.blast >> all.blast;done
seqkit fx2tab -j 30 -l  -n -i -H rice_ltc_ltg.fa | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' > rice_ltc_ltg_len
vi rice_ltc_ltg_len                                           ### %s#RLC_\#LTR/## %s#RLG_\#LTR/##
perl ~/bin/merge.pl rice_ltc_ltg_len 1 all.blast 1 o          ### Add fasta length on the blast result
awk '$8 != $10 && $7 !=$9' o |cut -f 2-12,15 > all.result     ### Filter the sequences blasted to itself
sort -k1,1 -k 6,6n all.result |awk '!a[$1]++{print}' > all.identity      ###Select the most left LTR as the identity


### We use TransDecoder to predict the CDS of TE and select the longest CDS to calculate GC content.
TransDecoder.LongOrfs -S -t rice_ltc_ltg.fa
blastp -query rice/rice_ltc_ltg.fa.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 15 > rice/rice_ltc_ltg.fa_longest_orfs_blastp.outfmt
hmmscan --cpu 20 --domtblout rice_ltc_ltg.fa.pfam.domtbloutq ~/Pfam-A.hmm rice_ltc_ltg.fa.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t rice_ltc_ltg.fa  --retain_pfam_hits rice_ltc_ltg.fa.pfam.domtblout --retain_blastp_hits rice_ltc_ltg.fa_longest_orfs_blastp.outfmt
perl filter_N_contain_fa.pl rice_ltc_ltg.fa.transdecoder.cds > rice_ltc_ltg.fa.transdecoder.cds.filteredN
seqkit fx2tab -j 30 -l -n -i -H rice_ltc_ltg.fa.transdecoder.cds.filteredN | cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' > rice_ltc_ltg.fa.transdecoder.cds.filteredN_length
sort -k1,1 -k 3,3nr rice_ltc_ltg.fa.transdecoder.cds.filteredN_length |awk '!a[$1]++{print}' > rice_ltc_ltg.fa.transdecoder.cds.filteredN_longest_length




