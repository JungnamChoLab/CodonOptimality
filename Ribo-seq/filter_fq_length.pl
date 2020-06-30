use strict;
use warnings;
### if reads content N,then throw it;
### useage: perl filter_fq_length.pl fq 18-26nt_fq 13nt_fq

open IN,"< $ARGV[0]";
my(%hash,$id,$qual,$seq,$jia,$desc,$length,$name,$n,$use,$nt13);

#open OUT,"> $ARGV[1]";
#open OUT1,"> $ARGV[2]";

local $/ = "S"; 
<IN>;
while (<IN>) {
s/\n\@S//;
($id,$seq,$jia,$qual) = split "\n", $_ , 4;
$name = '@S' . "$id";
$length = length($seq);
if (/\nN/) {
$n ++;
next;
#} elsif ($length == 13) {
#$nt13 ++;
#print OUT1 "$name\n$seq\n$jia\n$qual\n";
} elsif ($length >=40 && $length <= 65) {
$use ++;
#print OUT "$name\n$seq\n$jia\n$qual\n";
print  "$name\n$seq\n$jia\n$qual\n";
}
}
#my $total = $n + $nt13 + $use;
#my $total = $n + $nt13 + $use;
$/ = "\n";
#print "total reads\tN reads\tuse reads\t13nt reads\n$total\t$n\t$use\t$nt13\n";
close IN;
#close OUT;
#close OUT1;
#print "Done successfully\n";
