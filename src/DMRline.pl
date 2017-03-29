use strict;
my $dmrFile=$ARGV[0];
my $outFile="$dmrFile.line";
open DMR,"$dmrFile"||die"$!";
open DMRline,">$outFile"||die"$!";
print DMRline "track name=DMR description=\"DMRs fragment on genome\" color=255,0,0\n";
while(<DMR>){
	chomp;
	my @line=split(/\t/,$_);
	print DMRline "$line[0]\t$line[1]\t$line[2]\n";
}
close(DMR);
close(DMRline);
