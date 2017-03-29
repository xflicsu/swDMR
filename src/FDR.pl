use strict;
unless($ARGV[0]){
	die "perl $0 dmrfile R_BIN fdr_value\n";
}
my $data_file=$ARGV[0];
my $R=$ARGV[1];
my $fdr;
$fdr=$ARGV[2] if $ARGV[2];
$fdr||=0.05;
#-----------------------------#
#FDR Benjamini-Hochberg Method
#http://genomics.princeton.edu/storeylab/qvalue/manual.pdf
#-----------------------------#
open FDR,">$data_file.FDR.R" || die"$!";
print FDR "data=read.table(\"$data_file\",sep=\"\\t\");p=data[,5];qval=p.adjust(p,method=\"BH\",n=length(p));out=cbind(data,qval);write.table(out,file=\"$data_file.fdr\",row.names=F,col.names=F,quote=FALSE,sep=\"\\t\")";
`$R CMD BATCH $data_file.FDR.R`;
`rm ./*.FDR.Rout`;
close(FDR);
#_____________________________#
#select DMR satisfy FDR_value
#_____________________________#
open DMR, "$data_file.fdr"||die "$!";
open DMR_L, ">$data_file.fdr.low"||die "$!";
open DMR_H, ">$data_file.fdr.high"||die "$!";
while(<DMR>){
	chomp;
	my @line=split(/\s+/,$_);
	print DMR_L "$_\n" if @line[-1] <=$fdr;
	print DMR_H "$_\n" if @line[-1] > $fdr;
}

