use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/module/lib64/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib64/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib/perl";
use lib "$FindBin::Bin/../lib/module/share/perl";

use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(10);
my $pid;
my $number=0;
print "work start\n";
foreach (@ARGV){
	$number+=1;
	chomp;
	$pid = $pm->start and next;
	my $Commond="$_";
	system("$Commond");
	print "$Commond\n";
	$pm->finish;
}
$pm->wait_all_children;
print "work done\n";
