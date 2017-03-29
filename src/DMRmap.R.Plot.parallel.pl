=head Discription:
	plot DMR maps in parallel
=cut
#use strict;
print `date`;
use FindBin;
use File::Basename;
use Parallel::ForkManager;
unless(@ARGV==3){
	print "perl $0 <DMRmap.R.list> <Rbin> <noOfProcesses>\n";
	exit;
}
my $list=$ARGV[0];
my $Rbin=$ARGV[1];
my $noOfProcesses=$ARGV[2];
open FRAG,"$list"||die"$!";
my $pm = new Parallel::ForkManager($noOfProcesses);
my $pid;
while(<FRAG>){
	chomp;
	my $R=$_;
$pid = $pm->start and next;
	`$Rbin CMD BATCH $R`;
	$pm->finish;
}
$pm->wait_all_children;
`rm *.Rout`;
exit;
