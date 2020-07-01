# Rice CDS data was downloaded from MSU7.0 (http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/) named all.cds.txt
# The longest transcript was selected for each gene
seqkit fx2tab -j 30 -l  -n -i -H all.cds.txt |cut -f 1,4 |awk '{print$1"\t"$1"\t"$2}' > all.cds_len
vi all.cds_len  #:1dd        :%s#\(\d\+\)\.\d\+#\1#
sort -k1,1 -k3,3nr all.cds_len |awk '!a[$1]++{print}' > all.cds.longest_len
cut -f 2 all.cds.longest_len > longest_cds_trans_id
perl get_fa.pl longest_cds_trans_id all.cds.txt > longest_cds.txt
grep ">" longest_cds.txt |wc -l

# We calculate the suboptimal codon frequency in a window length of 20 condons, each window move 4 codons. 
# The first 200 codons and last 200 codons of each gene (length >= 1200bp) were calculated.

perl codon_frequency_first600.pl longest_cds.txt first600
perl codon_frequency_last600.pl longest_cds.txt last600
