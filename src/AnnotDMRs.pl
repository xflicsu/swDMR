use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/module/lib64/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib64/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib/perl";
use lib "$FindBin::Bin/../lib/module/share/perl";

use File::Basename;
use Annotation;
unless(@ARGV>=5){
	print STDERR "perl $0 DMRfile bedtools annotation outdir Rbin\n";
	exit;
}
my $dmrFile=$ARGV[0];
my $bedtools=$ARGV[1];
my $annotation=$ARGV[2];
my $outdir=$ARGV[3];
my $Rbin=$ARGV[4];
### overlap with annotation file
my ($Format,$FeatureIndex)=TestFormat(\$annotation);
my $baseDMR=basename($dmrFile);
my $baseAnno=basename($annotation);
`$bedtools -b $dmrFile -a $annotation -wb -wa >$outdir/Annotation/$baseDMR\VS$baseAnno.ove`;
my @element;
my $awkIndex;
if($Format eq "resistance"){
	print STDERR "Format is not gff/gtf/bed, please adjust it again!\n";
	exit;
}else{
	open IN,"$annotation"or die"$!";
	my %hash;
	while(<IN>){
		chomp;
		my @line=split(/\t/,$_);
		$hash{$line[$FeatureIndex]}=1;
	}
	close(IN);
	@element=keys %hash;
}
my %eleNum;
if(@element<=1){
	exit;
}else{
	open IN,"$outdir/Annotation/$baseDMR\VS$baseAnno.ove"or die"$!";
	while(<IN>){
		chomp;
		my @line=split(/\t/,$_);
		$eleNum{$line[$FeatureIndex]}+=1;
	}
	close(IN);
	open Table,">$outdir/Annotation/$baseDMR\VS$baseAnno.ove.table.xls";
	foreach(keys %eleNum){
		print Table "$_\t$eleNum{$_}\n";
	}
	close(Table);
	open R,">$outdir/Annotation/$baseDMR\VS$baseAnno.ove.table.xls.R";
	print R "dat=read.table(\"$outdir/Annotation/$baseDMR\VS$baseAnno.ove.table.xls\",sep=\"\\t\");
pdf(\"$outdir/Annotation/$baseDMR\VS$baseAnno.ove.table.pie.pdf\")
ratio=sprintf(\"%.3f\",100*dat[,2]/sum(dat[,2]))
ratio=paste(ratio,\"%\",sep=\"\")
a=rep(\"(\",4)
b=rep(\")\",4)
c=rep(\" \",4)
d=prettyNum(dat[,2],big.mark = \",\")
label=paste(dat[,1],a,d,c,ratio,b,sep=\"\")
pie(dat[,2],col=2:(length(dat[,1])+1),main=\"Pie for different element enrichment\",border=\"purple\",labels=label,font=1,cex=0.7)
dev.off()
";
	close(R);
	#print "$Rbin CMD BATCH $outdir/Annotation/$baseDMR\VS$baseAnno.ove.table.xls.R\n";
	`$Rbin CMD BATCH $outdir/Annotation/$baseDMR\VS$baseAnno.ove.table.xls.R`;
	`rm $baseDMR\VS$baseAnno.ove.table.xls.Rout`;
}
