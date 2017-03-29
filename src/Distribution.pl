#!/usr/bin/perl-w
use strict;
# distribution of swDMRs
die  "\nUsage: $0 <filter dmr><Outdir><pre> <Rbin>\n" unless (@ARGV >= 3);
open A,"$ARGV[0]"  || die "$!" ;
my $outdir = $ARGV[1];
my $pre = $ARGV[2] || "lowThanFDR";
my $Rbin=$ARGV[3];
my %hash;
my %num;
my $len;
my $flank=100;
my ($line,$lef,$rig,$left,$right,$chr);
my @info;
my %region;
chomp ($line=  <A>);
@info= split /\s+/,$line;
$chr = $info[0];
$lef = $info[1];
$rig = $info[2];
$len = ($rig-$lef+1);
$len = int($len/$flank);
$hash{$len}=1;
$num{$chr}=1;
#print $chr.$lef.$rig;
while(1) 
{ 
	last if (!$line);
	chomp ($line=  <A>);
	last if (!$line);
	@info= split /\s+/,$line;
	$len= $info[2]-$info[1]+1;
	$len = int($len/$flank);
	$hash{$len}++;
	$num{$info[0]}++;
	if($info[0] eq $chr){
		$left = $info[1];
		$right = $info[2];
		if ($rig<$left){
#			print OUT "$chr\t$lef\t$rig\t".($rig-$lef+1)."\n";
			$region{$chr} += ($rig-$lef+1);
			$lef = $left;
			$rig = $right;
		}
		elsif($rig>=$left && $rig<=$right){
			$rig = $right;
		}
		elsif($rig>$right){}
	}
	else{
#		print OUT "$chr\t$lef\t$rig\t".($rig-$lef+1)."\n";
		$region{$chr} += ($rig-$lef+1);
		$chr = $info[0];
		$lef = $info[1];
		$rig = $info[2];
	}
		
		
}
$region{$chr} += ($rig-$lef+1);
close A;
my $totalnum=0;
my $totalregion=0;
open OUT , ">$outdir/$pre.num_region.xls" or die "$!";
print OUT "chromosome\t# of swDMRs\tlength of swDMRs region\n";

foreach my $chr(sort keys %num){
	$totalnum += $num{$chr};
	$totalregion += $region{$chr};
	print OUT "$chr\t$num{$chr}\t$region{$chr}\n";
}
print OUT "total\t$totalnum\t$totalregion\n";
close OUT;
my $ylim=0;

open OUT ,">$outdir/$pre.length.dis" or die "$!";
foreach my $len(sort {$a<=>$b} keys %hash){
	my $key = $len * $flank;
	my $rate = $hash{$len}/$totalnum*100;
	if($ylim < $rate){
		$ylim = $rate;
	}
	print OUT "$key\t$hash{$len}\t$rate\n";
}
$ylim++;
$ylim = int $ylim;
close OUT;
open R, ">$outdir/$pre.length.dis.R" or die "$!";
print R "pdf(\"$outdir/$pre.length.dis.pdf\",height=6,width=8)
par(font.lab=1,font.axis=1,cex.lab=1.5,cex.axis=1.5,mar=c(3.4,3.4,1,0.5),mgp=c(2,0.5,0))
a<-read.table(\"$outdir/$pre.length.dis\")
plot(V3~V1,data=a,ylab=\"Percentage\",xlab=\"Length\",type=\"l\",col=\"red\",lwd=2,ylim=c(0,$ylim))
dev.off()
";
close R;
`$Rbin CMD BATCH $outdir/$pre.length.dis.R`;
`rm -f $pre.length.dis.Rout`;

