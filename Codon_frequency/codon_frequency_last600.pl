
use strict;
use warnings;
open IN,"< $ARGV[0]";
open OUT,"> $ARGV[1]";
my($tile,%id,$id,$seq,$length,$codon,%sub_codon,$fre,$count,@fre,%sub_count);

%sub_codon =(
'AAA' => 'S','AAT' => 'S','ACA' => 'S','ACT' => 'S',
'AGA' => 'S','AGT' => 'S','ATA' => 'S','ATT' => 'S',
'CAA' => 'S','CAT' => 'S','CCA' => 'S','CCT' => 'S',
'CTA' => 'S','CTT' => 'S','GAA' => 'S','GAT' => 'S',
'GCA' => 'S','GCT' => 'S','GTA' => 'S','GTT' => 'S',
'TAT' => 'S','TCA' => 'S','TCT' => 'S','TGT' => 'S',
'TTA' => 'S','TTG' => 'S','TTT' => 'S');

#open OUT,"> $ARGV[2]";
local $/ = ">"; 
<IN>;
while (my $line=<IN>) {
	chomp $line;
	($id,$seq) = split "\n", $line ,2;
	$id=~ s/ .*//;
	$seq =~ s/\s//g;
	my $length=length($seq);
	my $select=600;
	if ($length < $select*2) {
	next;
	}else{
		my $pos=$length-$select;
		$seq=substr($seq,$pos,$select);
#print $seq;
my $k=3; my $tile=60; my $step=12;
for(my $i = 0; $i <= $select - $tile;$i=$i+$step){
	my $subseq= substr($seq, $i, $tile);
#	 print  "\n",$subseq,"\n";
	# my $count = 0;
	for(my $j = 0; $j <= $tile - $k; $j = $j + 3) {
   		my $codon = substr($subseq, $j, $k);
#			print $codon,"\t";
		if(exists $sub_codon{$codon}) {
 	       		 $sub_count{$i}++;
#		print $sub_count{$i},",";
			}
	}
	$fre = $sub_count{$i}*300/$tile;
	push @fre,$fre;
	}
 print OUT $id,"\t","@fre","\n";  
%sub_count=();
@fre=();
}
}
$/ = "\n";
close IN;
close OUT;
#print "Done fq2fa successfully\n";


