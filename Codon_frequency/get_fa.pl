use strict;
use warnings;
open ID,"< $ARGV[0]";
open IN,"< $ARGV[1]";
my(%id,$id,$id1,$qual,$seq,$jia,$desc,$catchseq_seqio_obj1,$seq_obj,$name);

#open OUT,"> $ARGV[2]";
while (<ID>) {
  chomp;
  my @tmp=split/\t/,$_;
  $id{$tmp[0]}=1;
#print $tmp[0];
}
close ID;
local $/ = ">"; 
<IN>;
while (my $line=<IN>) {
chomp $line;
 if ($line=~/^(\S+) |/ && exists $id{$1}){
#print $1;
 print ">",$line;  
}else{
next;
}
}

$/ = "\n";
close IN;
#close OUT;
#print "Done fq2fa successfully\n";
