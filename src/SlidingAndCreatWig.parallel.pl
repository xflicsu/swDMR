=head Discription:
	1. Identification DMRs through sliding window from whole-genome bilsulfit sequnce data.
	2. Create wig(methylation level of whole-genome by step and span).
=head Edit log:
	Tue Mar  6 09:34:44 CST 2012
=head Change log:
	Sat May  5 13:11:12 CST 2012
	add condition:
		max_level - min_level >= threshold
=cut
#use strict;
print `date`;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/module/lib64/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib64/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib/perl";
use lib "$FindBin::Bin/../lib/module/share/perl";
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;
use statistics;

#**************************************#
my ($type,$sw,$step,$outdir,$dis,$samples,$name);
my ($p_threshold,$statistics,$coverage,$points,$diff);
my ($chrom,$pos,$methy,$unmethy,$ctype,$fold,$noOfProcesses);
#**************************************#
GetOptions(
	"type=s" => \$type, "samples=s" => \$samples,"chromosome=i" => \$chrom,
    "position=i" => \$pos,"methy=i" => \$methy,"unmethy=i" => \$unmethy,"ctype=i" => \$ctype,
   	"window=i" => \$sw, "step=i" => \$step,	"p_threshold=s" => \$p_threshold, "distance=i" => \$dis,
   	"outdir=s" => \$outdir,"points=i" => \$points,"fold=f" => \$fold,"diff=f" => \$diff,
	"statistics=s" => \$statistics, "coverage=i" => \$coverage,"name=s" => \$name,"Processes=i" => \$noOfProcesses,
);
unless($samples){
	print "Usage: you must input <samples>:\n";
	print "perl $0 --type CG --samples 1,2,3 --name 1,2,3 --chromosome 1 --position 2 --ctype 3 --methy 4 --unmethy 5 --window 1000 --points 5 --step 10 --p_threshold 0.05 --outdir ./ --statistics Fisher --coverage 4 --fold 2 --diff 0.1\n";
	exit;
}
#*****************#
$chrom-=1;$pos-=1;
$methy-=1;$unmethy-=1;
$ctype-=1;
#print "33333333333\t$chrom,$pos,$methy,$unmethy,$ctype\n";
#default value
mkdir "$outdir/tmp" unless -d "$outdir/tmp";
#*****************#
$type||="CG";$sw||=1000;
$points||=1;$fold||=2;
$p_threshold||=0.05;$step||=1;
$noOfProcesses||=1;
#print "22222222222222222222\t$step\n";
#**************************#
#open files
#**************************#
my @fileSam=split(/,/,$samples);
my $sam_num=@fileSam;
our (%fhand,$hand,%outhand);
our @NAME=split(/,/,$name);
open COLOR,"$FindBin::Bin/../etc/ColorRGBvalue.txt" or die "$!";
my @color;
while(<COLOR>){
	chomp;
	my @line=split(/\s+/,$_);
	chomp($line[0]);
	push @color, $line[0];
}
open_files(\@fileSam);
sub open_files{
	my $files=shift;
	my $num=@$files;
	for(0..$num-1){
		if(@$files[$_]=~/\S+.gz$/){open $fhand{$_},"<:gzip","@$files[$_]" or die $!;
		}else{open $fhand{$_},"@$files[$_]" or die "$!";}
		open $outhand{"$_"},">$outdir/$NAME[$_].wig" or die "$!";
	}
}
open ChrList,">$outdir/chr.list"||die"$!";
open FRAG,">$outdir/swDMR.tmp.frag"||die"$!";
#**************************#
#programe is running
#**************************#
my @TotalMethyReads;
my @position;
my $flag==0;
my $tmp_chr;
my $tmp_chr2;
my $exit=0;
my $min=1000;
my $trans;
my $transBefore;
my @samps_line;
my $long_step_flag=0;
my ($chr,$tmp_position);
my $flagN=0;
while(1){
	$flagN++;
#	print "1111111111\t",`date`,"\t111111111111\n";
	last if $exit==1;
	$long_step_flag=0,goto LONG_STEP if $long_step_flag==1;
	for(my $k=0;$k<$sam_num;$k++){
		my $h=$fhand{$k};
		my $tmp_line=<$h>;
		@{$samps_line[$k]}=split(/\s+/,$tmp_line);
		$exit=1,last if $tmp_line eq "";
#		print "$samps_line[$k][$methy]+$samps_line[$k][$unmethy]\n";
		$min=$samps_line[$k][$methy]+$samps_line[$k][$unmethy] if $min>=$samps_line[$k][$methy]+$samps_line[$k][$unmethy];
	}###read a sample line

	########################my test
	##############################################

	$tmp_position=$samps_line[0][$pos];
	last if $exit==1;
	#print "$tmp_position\n";
	$chr=$samps_line[0][$chrom];
	if($chr ne $tmp_chr){
		$flag=0;
		for(my $k=0;$k<$sam_num;$k++){
			$min=$samps_line[$k][$methy]+$samps_line[$k][$unmethy] if $min>=$samps_line[$k][$methy]+$samps_line[$k][$unmethy];
			my $handle=$outhand{"$k"};
			print $handle "track type=wiggle_0 name=\"$NAME[$k]\" description=\"methylation level mC/(mC+Un_mC)\" visibility=full autoScale=off viewLimits=0.0:1.0 color=$color[$k] yLineMark=11.76 yLineOnOff=on priority=10\nvariableStep chrom=$chr\n";
		}
		print ChrList "$chr\n";
		$tmp_chr = $chr;
	}
	if($flag==0){
		$transBefore=$trans;
		$trans=$tmp_position;
		$flag++;
	}
	$tmp_chr2=$chr if $flagN==1;
	if($type){
		if($type ne "C" && $samps_line[0][$ctype] ne "$type"){
			$min=1000;
			next;
		}
	}
	if($min<$coverage){
		$min=1000;
		next;
	}else{
		for(my $k=0;$k<$sam_num;$k++){
			my $rate=$samps_line[$k][$methy]/($samps_line[$k][$methy]+$samps_line[$k][$unmethy]);
			my $handle=$outhand{"$k"};
			print $handle "$samps_line[$k][$pos]\t$rate\n";
		}
	}

	LONG_STEP:
	
##########get points in a sliding window
	if($samps_line[0][$pos]-$trans<2*$sw&&$tmp_chr2 eq $chr){
		if($samps_line[0][$pos]-$trans+1<=$sw){
			for(my $k=0;$k<$sam_num;$k++){
				push @{$TotalMethyReads[$k][0]},$samps_line[$k][$methy];
				push @{$TotalMethyReads[$k][1]},$samps_line[$k][$unmethy];
			}
			push @position,$tmp_position;
			$long_step_flag=0;
		}else{
##### chose statistics methods to test DMRs
			my $pointNumber=@position;
			my @input;
			for(my $k=0;$k<$sam_num;$k++){
				push @input,[@{$TotalMethyReads[$k][0]}],[@{$TotalMethyReads[$k][1]}];
			}
############## test 
			if($pointNumber>=$points){
				my ($mean_rate,$maxRate,$minRate)=&statistics::make_mean_rate(\@input);
				my $foldMethy;
				if($minRate!=0){
					$foldMethy=$maxRate/$minRate;
				}else{
					$foldMethy=$fold;
				}
				my $diff_sample=$maxRate-$minRate;
				if($foldMethy>=$fold&&$diff_sample>=$diff){
					my $site2=$trans+$sw-1;
					my $out.="$chr\t$trans\t$site2\t$pointNumber\t";
					for(my $k=0;$k<$sam_num;$k++){
						for(my $m=0;$m<$pointNumber;$m++){
							$out.="$position[$m]\/$TotalMethyReads[$k][0][$m]\/$TotalMethyReads[$k][1][$m]";
							$out.="," if $m < $pointNumber-1;
						}
						$out.="\t" if $k<$sam_num-1;
					}
					print FRAG "$out\n";
				}
			}
#####shift steps
			my $index=0;
			foreach(@position){
				if($trans+$step>$_){
					$index++;
				}else{last;}
			}
			for(my $k=0;$k<$sam_num;$k++){
				@{$TotalMethyReads[$k][0]}=@{$TotalMethyReads[$k][0]}[$index..$pointNumber-1];
				@{$TotalMethyReads[$k][1]}=@{$TotalMethyReads[$k][1]}[$index..$pointNumber-1];
			}
			@position=@position[$index..$pointNumber-1];
			$trans+=$step;
			$long_step_flag=1;
		}
	}else{
		if($tmp_chr2 ne $chr){
			my $pointNumber=@position;
			my @input;
			for(my $k=0;$k<$sam_num;$k++){
				push @input,[@{$TotalMethyReads[$k][0]}],[@{$TotalMethyReads[$k][1]}];
			}
############## test 
			if($pointNumber>=$points){
				my ($mean_rate,$maxRate,$minRate)=&statistics::make_mean_rate(\@input);
				my $foldMethy;
				if($minRate!=0){
					$foldMethy=$maxRate/$minRate;
				}else{
					$foldMethy=$fold;
				}
				my $diff_sample=$maxRate-$minRate;
				if($foldMethy>=$fold&&$diff_sample>=$diff){
					my $site2=$transBefore+$sw-1;
					my $out.="$tmp_chr2\t$transBefore\t$site2\t$pointNumber\t";
					for(my $k=0;$k<$sam_num;$k++){
						for(my $m=0;$m<$pointNumber;$m++){
							$out.="$position[$m]\/$TotalMethyReads[$k][0][$m]\/$TotalMethyReads[$k][1][$m]";
							$out.="," if $m < $pointNumber-1;
						}
						$out.="\t" if $k<$sam_num-1;
					}
					print FRAG "$out\n";
				}
			}
		}
		@position=();
		@TotalMethyReads=();
		for(my $k=0;$k<$sam_num;$k++){
			$TotalMethyReads[$k][0][0]=$samps_line[$k][$methy];
			$TotalMethyReads[$k][1][0]=$samps_line[$k][$unmethy];
		}
		$position[0]=$tmp_position;
		$trans=$tmp_position;
		$long_step_flag=0;
		$tmp_chr2 = $chr;
	}
}

### the last one for test ###
my $pointNumber=@position;
my @input;
for(my $k=0;$k<$sam_num;$k++){
	push @input,[@{$TotalMethyReads[$k][0]}],[@{$TotalMethyReads[$k][1]}];
}
############## test 
if($pointNumber>=$points){
	my ($mean_rate,$maxRate,$minRate)=&statistics::make_mean_rate(\@input);
	my $foldMethy;
	if($minRate!=0){
		$foldMethy=$maxRate/$minRate;
	}else{
		$foldMethy=$fold;
	}
	my $diff_sample=$maxRate-$minRate;
	if($foldMethy>=$fold&&$diff_sample>=$diff){
		my $site2=$trans+$sw-1;
		my $out.="$chr\t$trans\t$site2\t$pointNumber\t";
		for(my $k=0;$k<$sam_num;$k++){
			for(my $m=0;$m<$pointNumber;$m++){
				$out.="$position[$m]\/$TotalMethyReads[$k][0][$m]\/$TotalMethyReads[$k][1][$m]";
				$out.="," if $m < $pointNumber-1;
			}
			$out.="\t" if $k<$sam_num-1;
		}
		print FRAG "$out\n";
	}
}


close(FRAG);
@position=();

open FRAG,"$outdir/swDMR.tmp.frag"||die"$!";
open Fraglist,">$outdir/swDMR.tmp.frag.list"||die"$!";
my $pm = new Parallel::ForkManager($noOfProcesses);
my $pid;
my $n=0;
my %repeatRegion;
while(<FRAG>){
	$n++;
	chomp;
	my $line=$_;
	last if $line eq "";
	my @dmrs=split(/\s+/,$line);
	my $chr=$dmrs[0];
	my $trans=$dmrs[1];
	my $site2=$dmrs[2];
	my $pointNumber=$dmrs[3];
	my @temp_join=@dmrs[4..(4+$sam_num-1)];
	my @input;
	my @TotalMethyReads;
	my @pos=split(/,/,$temp_join[0]);
	my $points=@pos;
	my @position_2;
	my @sites=split(/,/,$dmr[4]);
	my @bgsite=split(/\//,$sites[0]);
	my @edsite=split(/\//,$sites[-1]);
	if(exists $repeatRegion{"$chr-$pointNumber-$bgsite[0]-$edsite[0]"}){
		next;
	}else{
		%repeatRegion=();
		$repeatRegion{"$chr-$pointNumber-$bgsite[0]-$edsite[0]"}=1;
	}
	foreach(@pos){
		my @tmp=split(/\//,$_);
		push @position_2,$tmp[0];
	}
	for(my $k=0;$k<$sam_num;$k++){
		my @ttmp=split(/,/,$temp_join[$k]);
		for(my $i=0;$i<$points;$i++){
			my @site_methy=split(/\//,$ttmp[$i]);
			$TotalMethyReads[$k][0][$i]=$site_methy[1];
			$TotalMethyReads[$k][1][$i]=$site_methy[2];
		}
	}
	for(my $k=0;$k<$sam_num;$k++){
		push @input,[@{$TotalMethyReads[$k][0]}],[@{$TotalMethyReads[$k][1]}];
	}
	print Fraglist "$outdir/tmp/tmp.$n\n";
	
$pid = $pm->start and next;
	my ($p_value)=&{$statistics}(\@input);
#	if($p_value<=$p_threshold){
		my $out.="$chr\t$trans\t$site2\t$pointNumber\t$p_value\t";
		for(my $k=0;$k<$sam_num;$k++){
			for(my $m=0;$m<$pointNumber;$m++){
				$out.="$position_2[$m]\/$TotalMethyReads[$k][0][$m]\/$TotalMethyReads[$k][1][$m]";
				$out.="," if $m < $pointNumber-1;
			}
			$out.="\t" if $k<$sam_num-1;
		}
		open DMR,">$outdir/tmp/tmp.$n"||die"$!";
		print DMR "$out\n";
		close(DMR);
#	}
#	print "$p_value\n";
#	print "abc\n";
	$pm->finish;
}
$pm->wait_all_children;
close(ChrList);
sleep 2;
close(Fraglist);
open Fraglist,"$outdir/swDMR.tmp.frag.list"||die"$!";
`rm -f $outdir/swDMR.tmp`;
while(<Fraglist>){
	chomp;
	`cat $_ >>$outdir/swDMR.tmp`;
}
`sort -k1,1 -k2,2n $outdir/swDMR.tmp|uniq >$outdir/swDMR.tmp.2;mv $outdir/swDMR.tmp.2 $outdir/swDMR.tmp`;
`rm $outdir/swDMR.tmp.frag.list;rm $outdir/swDMR.tmp.frag;rm -rf $outdir/tmp`;
exit;
