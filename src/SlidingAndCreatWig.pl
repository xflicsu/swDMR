=head Discription:
	1. Identification DMRs through sliding window from whole-genome bilsulfit sequnce data.
	2. Create wig(methylation level of whole-genome by step and span).
=head Edit log:
	Tue Mar  6 09:34:44 CST 2012
=head Change log:
	
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
use statistics;

#**************************************#
my ($type,$sw,$step,$outdir,$dis,$samples,$name);
my ($p_threshold,$statistics,$coverage,$points);
my ($chrom,$pos,$methy,$unmethy,$ctype,$fold);
#**************************************#
GetOptions(
	"type=s" => \$type, "samples=s" => \$samples,"chromosome=i" => \$chrom,
    "position=i" => \$pos,"methy=i" => \$methy,"unmethy=i" => \$unmethy,"ctype=i" => \$ctype,
   	"window=i" => \$sw, "step=i" => \$step,	"p_threshold=s" => \$p_threshold,
   	"distance=i" => \$dis, "outdir=s" => \$outdir,"points=i" => \$points,"fold=i" => \$fold,
	"statistics=s" => \$statistics, "coverage=i" => \$coverage,"name=s" => \$name,
);
unless($samples){
	print "Usage: you must input <samples>:\n";
	print "perl $0 --type CG --samples 1,2,3 --name 1,2,3 --chromosome 1 --position 2 --ctype 3 --methy 4 --unmethy 5 --window 1000 --points 5 --step 10 --p_threshold 0.05 --outdir ./ --statistics fisher --coverage 4 --fold 2\n";
	exit;
}
#*****************#
$chrom-=1;$pos-=1;
$methy-=1;$unmethy-=1;
$ctype-=1;
#print "33333333333\t$chrom,$pos,$methy,$unmethy,$ctype\n";
#default value
#*****************#
$type||="CG";$sw||=1000;
$points||=1;$fold||=2;
$p_threshold||=0.05;$step||=1;
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
open DMR,">$outdir/swDMRs.tmp"||die"$!";
#**************************#
#programe is running
#**************************#
my @TotalMethyReads;
my @position;
my $flag==0;
my $tmp_chr;
my $exit=0;
my $min=1000;
my $trans;
while(1){
#	print "1111111111\t",`date`,"\t111111111111\n";
	last if $exit==1;
	my @samps_line;
	for(my $k=0;$k<$sam_num;$k++){
		my $h=$fhand{$k};
		my $tmp_line=<$h>;
		@{$samps_line[$k]}=split(/\s+/,$tmp_line);
		$exit=1,last if $tmp_line eq "";
#		print "$samps_line[$k][$methy]+$samps_line[$k][$unmethy]\n";
		$min=$samps_line[$k][$methy]+$samps_line[$k][$unmethy] if $min>=$samps_line[$k][$methy]+$samps_line[$k][$unmethy];
	}###read a sample line
	my $tmp_position=$samps_line[0][$pos];
	last if $exit==1;
	#print "$tmp_position\n";
	my $chr=$samps_line[0][$chrom];
	if($chr ne $tmp_chr){
		@TotalMethyReads=();
		$flag==0;
		for(my $k=0;$k<$sam_num;$k++){
#			$min=$samps_line[$k][$methy]+$samps_line[$k][$unmethy] if $min>=$samps_line[$k][$methy]+$samps_line[$k][$unmethy];
			my $handle=$outhand{"$k"};
			print $handle "track type=wiggle_0 name=\"$NAME[$k]\" description=\"variableStep format\" visibility=full autoScale=off viewLimits=0.0:1.0 color=$color[$k] yLineMark=11.76 yLineOnOff=on priority=10\nvariableStep chrom=$chr\n";
		}
	}
	if($flag==0){
		$trans=$tmp_position;
		$flag++;
	}
	$tmp_chr=$chr;
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

##########get points in a sliding window
	if($samps_line[0][$pos]-$trans+1<=$sw){
		for(my $k=0;$k<$sam_num;$k++){
			push @{$TotalMethyReads[$k][0]},$samps_line[$k][$methy];
			push @{$TotalMethyReads[$k][1]},$samps_line[$k][$unmethy];
		}
		push @position,$tmp_position;
	}else{
##### chose statistics methods to test DMRs
#		print "2222222222\t",`date`,"\t222222222222\n";
		my $pointNumber=@position;
		my @input;
		for(my $k=0;$k<$sam_num;$k++){
			push @input,[@{$TotalMethyReads[$k][0]}],[@{$TotalMethyReads[$k][1]}];
		}
############## test 
		#print "@input[0]";
		if($pointNumber>=$points){
			my ($mean_rate,$maxRate,$minRate)=&statistics::make_mean_rate(\@input);
			my $foldMethy;
			if($minRate!=0){
				$foldMethy=$maxRate/$minRate;
			}else{
				$foldMethy=$fold;
			}
			if($foldMethy>=$fold){
				my ($p_value)=&{$statistics}(\@input);
#				if($p_value<=$p_threshold){
#					my $point_num=@{$TotalMethyReads[0][0]};
					my $site2=$trans+$sw-1;
					my $out.="$chr\t$trans\t$site2\t$pointNumber\t$p_value\t";
					for(my $k=0;$k<$sam_num;$k++){
						for(my $m=0;$m<$pointNumber;$m++){
							$out.="$position[$m]\/$TotalMethyReads[$k][0][$m]\/$TotalMethyReads[$k][1][$m]";
							$out.="," if $m < $pointNumber-1;
						}
						$out.="\t" if $k<$sam_num-1;
					}
					print DMR "$out\n";
#				}
			}
		}
#####shift steps
#		print "33333333333333\t",`date`,"\t333333333333333\n";
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
		if($samps_line[0][$pos]-$trans+1<=$sw){
			for(my $k=0;$k<$sam_num;$k++){
				push @{$TotalMethyReads[$k][0]},$samps_line[$k][$methy];
				push @{$TotalMethyReads[$k][1]},$samps_line[$k][$unmethy];
			}
			push @position,$tmp_position;
		}else{
			@position=();
			@TotalMethyReads=();
			for(my $k=0;$k<$sam_num;$k++){
				$TotalMethyReads[$k][0][0]=$samps_line[$k][$methy];
				$TotalMethyReads[$k][1][0]=$samps_line[$k][$unmethy];
			}
			$position[0]=$tmp_position;
			$trans=$tmp_position;
		}
#		print "444444444444444\t",`date`,"\t444444444444444444\n";
	}
}
print `date`;
