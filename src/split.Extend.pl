=head Discription:
	Identification DMRs through sliding window from whole-genome bilsulfit sequnce data.
	Then extend the regions
=head Edit log:
	Sun Feb 26 19:45:51 CST 2012
=head Change log:

=head Commend example
	$0 $dmrdir $sam_num $len $statistics
		$dmrdir : Slinding result
 		$sam_num: samples
		$len	: join length
		$statistics method
=cut
#use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/module/lib64/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib64/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib/perl";
use lib "$FindBin::Bin/../lib/module/share/perl";
use Parallel::ForkManager;
my $extendScript="$FindBin::Bin/Extend.2.pl";

unless(@ARGV==4){
	print "Usage: split potential DMR into little pieces to parallel test ( split by distance larger than defined length )
 	perl $0 <dmrDir> <distance> <processes> <parameters>\n";
	exit;
}
my $dmrdir=$ARGV[0];
my $distance=$ARGV[1];
my $num=$ARGV[2];
my $processes=$num;
my $parameters=$ARGV[3];

open DMR,"$dmrdir/swDMR.tmp.fdr.low"||die"$!";
my $rows=`less $dmrdir/swDMR.tmp.fdr.low|wc -l`;
system(mkdir "$dmrdir/tmpdmr");
my $averageRow=$rows/$num;
my $n=0;
my $m=0;
my $flag=1;
my @tmp;
open LIST,">$dmrdir/tmpdmr/dmr.list"||die"$!";
while(<DMR>){
	$n+=1;
	chomp;
	my $line=$_;
	if($n>=$averageRow){
		my $chr=(split(/\t/,$_))[0];
		my $start=(split(/\t/,$_))[1];
		my $end=(split(/\t/,$_))[2];
		if(($chr ne $tmp[0])||($chr eq $tmp[0] && $tmp[2]-$end>$distance)){
			$flag=1;
			$n=0;
		}
	}

	if($flag==1){
		close(OUT);
		$m+=1;
		open OUT,">$dmrdir/tmpdmr/dmr.$m";
		print LIST "$dmrdir/tmpdmr/dmr.$m\n";
		$flag=0;
		print ">$dmrdir/tmpdmr/dmr.$m\n";
	}
	print OUT "$line\n";
	@tmp=($chr,$start,$end);
}

### parallel work
open LIST,"$dmrdir/tmpdmr/dmr.list"||die"$!";
my $pm = new Parallel::ForkManager($processes);
my $pid;
my $number=0;
print "Extend work start\n";

while(<LIST>){
	chomp;
	$number+=1;
	$pid = $pm->start and next;
	my $Commond="perl $extendScript $parameters $_ ";
	system("$Commond");
	print "$Commond\n";
	$pm->finish;
}
$pm->wait_all_children;

close(LIST);

open LIST,"$dmrdir/tmpdmr/dmr.list"||die"$!";
my $catCMD="cat ";
while(<LIST>){
	chomp;
	$catCMD.=" $_.Extend ";
}
$catCMD.=" >$dmrdir/swDMR.fdr.low.Extend";
close(LIST);
system($catCMD);
system("rm -fr $dmrdir/tmpdmr/");
system("sort -k1,1 -k2,2n $dmrdir/swDMR.fdr.low.Extend -o $dmrdir/swDMR.fdr.low.Extend");
open ChrList,"$dmrdir/chr.list"||die"$!";
while(<ChrList>){
	chomp;
	`grep -w $_ $dmrdir/swDMR.fdr.low.Extend >>$dmrdir/swDMR.fdr.low.Extend.new;`;
}
`mv $dmrdir/swDMR.fdr.low.Extend.new $dmrdir/swDMR.fdr.low.Extend`;

