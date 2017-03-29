#!/usr/bin/perl -w
use strict;
use FindBin;
unless(@ARGV==4){
	print "perl $0 name data outdir Rbin\n";
	exit;
}
my $Rlib="$FindBin::Bin/../lib/R-packages";
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
open OUT, ">$outdir/swDMR.heatmap.R"||die"$!";
print OUT "
library(\"bitops\",lib.loc=\"$Rlib\")
library(\"caTools\",lib.loc=\"$Rlib\")
library(\"gdata\",lib.loc=\"$Rlib\")
library(\"gtools\",lib.loc=\"$Rlib\")
library(\"gplots\",lib.loc=\"$Rlib\")
pdf(\"$outdir/swDMR.heatmap.pdf\")
rt<-read.table(\"$data\")
par(oma=c(1.5,2,2,2));
Lab.pallete<-colorRampPalette(c(\"lightgrey\",\"yellow\",\"blue\"),space=\"Lab\")
heatmap.2(as.matrix(rt[,5:$index]),col=Lab.pallete(300),Colv=T,Rowv=T,dendrogram=\"both\",key=T,keysize=1.2,labCol=c(";
foreach(@NAME[0..$num-2]){
	print OUT "\"$_\",";
}
print OUT "\"$NAME[$num-1]\"";
print OUT "),cexCol=1.2,cexRow=1.2,labRow=\"\",trace=\"none\",colsep=1:$num,sepwidth=rep(0.0005,times=$num),main=\"Heatmap of DMRs\")
dev.off()";

`$Rbin CMD BATCH $outdir/swDMR.heatmap.R&`;
`rm swDMR.heatmap.Rout`;
