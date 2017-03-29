#!/usr/bin/perl -w
use strict;
use FindBin;
unless(@ARGV==4){
	print "perl $0 name data outdir Rbin\n";
	exit;
}
#print $Rlib;
my $name=$ARGV[0];
my $data=$ARGV[1];
my $outdir=$ARGV[2];
my $Rbin=$ARGV[3];
my @NAME=split(/,/,$name);
my $num=@NAME;
my $index=5+$num-1;
mkdir "$outdir" unless -d "$outdir";
print "$outdir\n"; 
open OUT, ">$outdir/swDMR.boxplot.R"||die"$!";
print OUT "
pdf(\"$outdir/swDMR.boxplot.pdf\")
rt<-read.table(\"$data\")
par(oma=c(1.5,2,2,2));
boxplot(";
foreach(0..$num-1){
	$_+=5;
	print OUT "rt[,$_],";
}
print OUT "main=\"Boxplot of swDMR\",boxwex=0.1,at=1:$num,ylab=\"Methylation level\",names=c(";
foreach(@NAME[0..$num-2]){
	print OUT "\"$_\",";
}
print OUT "\"$NAME[$num-1]\"";
print OUT "),col=rainbow($num),las=3)
legend(\"topright\",c(";
foreach(@NAME[0..$num-2]){
	print OUT "\"$_\",";
}
print OUT "\"$NAME[$num-1]\"";
print OUT ")
fill=rainbow($num);
dev.off()";

`$Rbin CMD BATCH $outdir/swDMR.boxplot.R&`;
`rm swDMR.boxplot.Rout`;
