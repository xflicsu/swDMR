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
use statistics;
unless(@ARGV>=4){
	print "Usage: you should input more than three args
 	perl $0 dmrDir samNumber statistics length fold diff p_threshold point_threshold\n";
	exit;
}
my $dmrdir=$ARGV[0];
my $sam_num=$ARGV[1];
my $statistics=$ARGV[2];
my $len=$ARGV[3];
my $fold=$ARGV[4]||2;
my $diff=$ARGV[5]||0.1;
my $p_threshold=$ARGV[6]||0.01;
my $points_threshold=$ARGV[7]||5;
open DMR,"$dmrdir/swDMR.tmp.fdr.low"||die"$!";
open DMROUT,">$dmrdir/swDMR.tmp.fdr.low.Extend"||die"$!";
open ChrList,"$dmrdir/chr.list"||die"$!";
#### produce the first vector of DMR ####
my $line=<DMR>;
chomp($line);
my $old=$line;
my @dmrs=split(/\s+/,$line);
my $tmp_chr=$dmrs[0];
my @methy_first=@dmrs[5..(5+$sam_num-1)];
my @temp_join=@methy_first;

my @site_tmp=split(/,/,$temp_join[0]);
my @site_first=split(/\//,$site_tmp[0]);
my $start=$site_first[0];
my @site_end=split(/\//,$site_tmp[-1]);
my $end=$site_end[0];
my $point_tmp;

#### define some variables ####
my $point_site;my $point_site_tmp;my $pval;
my $flag2=0;
my $methy_levels_before;
my $mean_dep_before;
my $leng_tmp;
my @temp_join_last;
my $chr;
while(<DMR>){
#### read second vector of DMR ####
	chomp;
	$line=$_;
	last if $line eq "";
	@dmrs=split(/\s+/,$line);
	$chr=$dmrs[0];
	my @methy_second=@dmrs[5..(5+$sam_num-1)];
	@site_tmp=split(/,/,$temp_join[0]);
	@site_first=split(/\//,$site_tmp[0]);
	@site_end=split(/\//,$site_tmp[-1]);
	my @methy_1=split(/,/,$temp_join[0]);
	my @temp_join_bak=@temp_join;
	my $point_site_1=@methy_1;
	my @methy_1_site=split('/',$methy_1[-1]);
	$start=$site_first[0];
	$end=$site_end[0];
	
	my @methy_2=split(/,/,$methy_second[0]);
	my @methy_2_site=split('/',$methy_2[0]);
	
##########################################
#### join vector with overlap ####
	if($methy_2_site[0]-$end+1<=$len && $chr eq $tmp_chr){
#### find the point with same coordinate ####
		$point_site=@methy_2;
		for(my $i=0;$i<$point_site;$i++){
			my @methy_2_site=split('/',$methy_2[$i]);
			if($methy_2_site[0]==$methy_1_site[0]){
				$point_site_tmp=$i+1;
				last;
			}elsif($i==0 && $methy_1_site[0]<$methy_2_site[0]){
				$point_site_tmp=0;
				last;
			}
		}
#### create new vectors in potential DMR ####
		for(my $k=0;$k<$sam_num;$k++){
			my @temp1=split(/,/,$temp_join[$k]);
			my @temp2=split(/,/,$methy_second[$k]);
			my @tmp;
			$leng_tmp=@temp2;
			if($point_site_tmp!=$leng_tmp-1){
				@tmp=(@temp1,@temp2[$point_site_tmp..$point_site-1]);
				$temp_join_last[$k]=join ",",@temp2[$point_site_tmp..$leng_tmp-1];
				$temp_join[$k]=join ",",@tmp;
			}else{
				@tmp=@temp1;
				$temp_join_last[$k]="";
				$temp_join[$k]=join ",",@tmp;
			}
			$point_tmp=@tmp;
		}
		
##########################################
#### prepare eligible data to do statistics test ####
		
		if($point_tmp>=$points_threshold){
			my @input;
			my @TotalMethyReads;
			my $points=split(/,/,$temp_join[0]);
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
			my ($methy_levels,$maxRate,$minRate,$mean_dep);
			($methy_levels,$maxRate,$minRate,$mean_dep)=&statistics::make_mean_rate(\@input);
			my $foldMethy;
			if($minRate!=0){
				$foldMethy=$maxRate/$minRate;
			}else{
				$foldMethy=$fold;
			}
		
			my $diff_sample=$maxRate-$minRate;
			$pval=&{$statistics}(\@input);
			unless($foldMethy>=$fold&&$diff_sample>=$diff&&$pval<=$p_threshold){
#### if all condition is ok, keep going; if not, print the before DMR ####
				my @input_bak;
				my @TotalMethyReads_bak;
				my @points_bak_tmp=split(/,/,$temp_join_bak[0]);
				my $points_bak=@points_bak_tmp;
				if($points_bak>=$points_threshold){
					for(my $k=0;$k<$sam_num;$k++){
						my @ttmp=split(/,/,$temp_join_bak[$k]);
						for(my $i=0;$i<$points;$i++){
							my @site_methy=split(/\//,$ttmp[$i]);
							$TotalMethyReads_bak[$k][0][$i]=$site_methy[1];
							$TotalMethyReads_bak[$k][1][$i]=$site_methy[2];
						}
					}
					for(my $k=0;$k<$sam_num;$k++){
						push @input_bak,[@{$TotalMethyReads_bak[$k][0]}],[@{$TotalMethyReads_bak[$k][1]}];
					}
					my ($methy_levels_bak,$maxRate_bak,$minRate_bak,$mean_dep_bak);
					($methy_levels_bak,$maxRate_bak,$minRate_bak,$mean_dep_bak)=&statistics::make_mean_rate(\@input_bak);
					my $foldMethy_bak;
					if($minRate_bak!=0){
						$foldMethy_bak=$maxRate_bak/$minRate_bak;
					}else{
						$foldMethy_bak=$fold;
					}
					my $diff_sample_bak=$maxRate_bak-$minRate_bak;
					$pval=&{$statistics}(\@input);
					if($foldMethy_bak>=$fold&&$diff_sample_bak>=$diff&&$pval<=$p_threshold){
						print DMROUT "$dmrs[0]\t$start\t$end\t$points_bak\t";
						print DMROUT join "\t", @$methy_levels_bak;
						print DMROUT "\t";
						print DMROUT join "\t", @$mean_dep_bak;
						print DMROUT "\t$pval\n";
					}
				}
#### change the last points	to prepare to do other DMR calling ####
				for(my $k=0;$k<$sam_num;$k++){
					$temp_join[$k]=$temp_join_last[$k];
				}
			}
		}
#########################################
	}else{
#### there is no overlap region in both first and second vector data set, so do statistics test ####
		if($point_site_1>=$points_threshold){
			my @input;
			my @TotalMethyReads;
			my $points=split(/,/,$temp_join[0]);
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
	
			my ($methy_levels,$maxRate,$minRate,$mean_dep);
			($methy_levels,$maxRate,$minRate,$mean_dep)=&statistics::make_mean_rate(\@input);
			my $foldMethy;
			if($minRate!=0){
				$foldMethy=$maxRate/$minRate;
			}else{
				$foldMethy=$fold;
			}
			my $diff_sample=$maxRate-$minRate;
			$pval=&{$statistics}(\@input);
			if($foldMethy>=$fold&&$diff_sample>=$diff&&$pval<=$p_threshold){
				print DMROUT "$tmp_chr\t$start\t$end\t$points\t";
				print DMROUT join "\t", @$methy_levels;
				print DMROUT "\t";
				print DMROUT join "\t", @$mean_dep;
				print DMROUT "\t$pval\n";
			}
		}
		$tmp_chr=$chr;
		@temp_join=@methy_second;
#########################################		
	}
}

my @TotalMethyReads;
my $points=split(/,/,$temp_join[0]);

my @site_tmp=split(/,/,$temp_join[0]);
my @site_first=split(/\//,$site_tmp[0]);
my $start=$site_first[0];
my @site_end=split(/\//,$site_tmp[-1]);
my $end=$site_end[0];

if($points>=$points_threshold){
	for(my $k=0;$k<$sam_num;$k++){
		my @ttmp=split(/,/,$temp_join[$k]);
		for(my $i=0;$i<$points;$i++){
			my @site_methy=split(/\//,$ttmp[$i]);
			$TotalMethyReads[$k][0][$i]=$site_methy[1];
			$TotalMethyReads[$k][1][$i]=$site_methy[2];
		}
	}
	my @input;
	for(my $k=0;$k<$sam_num;$k++){
		push @input,[@{$TotalMethyReads[$k][0]}],[@{$TotalMethyReads[$k][1]}];
	}
	my ($methy_levels,$maxRate,$minRate,$mean_dep);
	($methy_levels,$maxRate,$minRate,$mean_dep)=&statistics::make_mean_rate(\@input);
	my $foldMethy;
	if($minRate!=0){
		$foldMethy=$maxRate/$minRate;
	}else{
		$foldMethy=$fold;
	}
	my $diff_sample=$maxRate-$minRate;
	$pval=&{$statistics}(\@input);
	if($foldMethy>=$fold&&$diff_sample>=$diff&&$pval<=$p_threshold){
		print DMROUT "$dmrs[0]\t$start\t$end\t$points\t";
		print DMROUT join "\t", @$methy_levels;
		print DMROUT "\t";
		print DMROUT join "\t", @$mean_dep;
		print DMROUT "\t$pval\n";
	}
}
close(DMROUT);

#### sort by the list of DMR ####
open DMROUT,">$dmrdir/swDMR.tmp.fdr.low.Extend.new"||die"$!";
`sort -k1,1 -k2,2n $dmrdir/swDMR.tmp.fdr.low.Extend -o $dmrdir/swDMR.tmp.fdr.low.Extend`;
while(<ChrList>){
	chomp;
	`grep -w $_ $dmrdir/swDMR.tmp.fdr.low.Extend >>$dmrdir/swDMR.tmp.fdr.low.Extend.new;`;
}
`mv $dmrdir/swDMR.tmp.fdr.low.Extend.new $dmrdir/swDMR.fdr.low.Extend`;
`rm $dmrdir/chr.list;rm $dmrdir/swDMR.tmp`;
close(DMROUT);
