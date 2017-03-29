#use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/module/lib64/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib64/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5";
use lib "$FindBin::Bin/../lib/module/lib/perl5/site_perl";
use lib "$FindBin::Bin/../lib/module/lib/perl";
use lib "$FindBin::Bin/../lib/module/share/perl";
use File::Basename qw(dirname);
use Getopt::Long;
use swDMRsMap;
use PerlIO::gzip;

my ($help,$dmrFile1,$dmrFile2,$wigData,$outdir,$noOfProcesses,$gene,$CGI,$color,$left,$right);
my ($samples,$type,$chromo,$ctype,$methy,$unmethy,$pos,$chrom,$coverage,$samplesName,$Rbin);
GetOptions(
	"h|help|?" => \$help,
	"dmr1=s" => \$dmrFile1,
	"dmr2=s" => \$dmrFile2,
	"outdir=s" => \$outdir,
	"gene=s" => \$gene,
	"CGI=s" => \$CGI,
	"color=s" => \$color,
	"left=i" => \$left,
	"right=i" => \$right,
	"samples=s" => \$samples,
	"addName=s" => \$samplesName,
	"type=s" => \$type,
	"ctype=i" => \$ctype,
	"methy=i" => \$methy,
	"unmethy=i" => \$unmethy,
	"pos=i" => \$pos,
	"chromo=i" => \$chromo,
	"cov=s" => \$coverage,
	"Rbin=s" => \$Rbin,
);

$chromo-=1;$pos-=1;
$methy-=1;$unmethy-=1;
$ctype-=1;
if($help||!$dmrFile1){
	print "\nperl $0 -dmr1 <DMR flie> -dmr2 <DMR file> -s <sample1,sample2,sample3,,,> <-a sampleNames> -type CG -methy 6 -unmethy 7 -pos 2 -chromo 1 -o <outdir> -g <gene.bed> -C <CGI.bed> -c <etc/color> -l <left> -right <right> --Rbin /R/Path/R\n\n";
	exit
}
my @samName=split(/,/,$samplesName);
my $Script="$FindBin::Bin";
my $dir="$outdir/DMRsMap";
mkdir "$outdir/DMRsMap" unless -d "$outdir/DMRsMap";

### draw map of DMR (SVG) with parallel processes ###
open DMR,"$dmrFile1"||die"$!";
### read the color palette ###
open COLOR,"$FindBin::Bin/../etc/ColorRGBvalue.txt" or die "$!";
my @color;
while(<COLOR>){
	chomp;
	my @line=split(/\s+/,$_);
	chomp($line[2]);
	push @color, $line[2];
}

my @fileSam=split(/,/,$samples);
my $sam_num=@fileSam;
my %fhand;my $hand;
for(0..$sam_num-1){
	if($fileSam[$_]=~/\S+.gz$/){
		open $fhand{$_}, "<:gzip", "$fileSam[$_]" or die $!;
	}else{
		open $fhand{$_}, "$fileSam[$_]" or die "$!";
	}
}

####data get
my @methy_y;
my @density;
my $exit=0;
my $min=1000;
my $flagDMR=1;
my $dmrBeforeFlag=0;
my $flagMethy=1;
my @samps_line;
my @beforeDMR;
my @dmr;
my $tmpPoint;
my $tmpChr;
open DMRmapList,">$dir/DMRmap.R.list"||die"$!";
while(1){
	last if $exit==1;
	if($flagMethy==1){
		$tmpPoint=$samps_line[0][$pos];
		for(my $k=0;$k<$sam_num;$k++){
			my $hand=$fhand{$k};
			my $tmp_line=<$hand>;
			@{$samps_line[$k]}=split(/\s+/,$tmp_line);
			$exit=1,last if $samps_line[$k][0] eq "";
			$min=$samps_line[$k][$methy]+$samps_line[$k][$unmethy] if $min>=$samps_line[$k][$methy]+$samps_line[$k][$unmethy];
		}
		last if $exit==1;
		if($type){
			if($type ne "C" && $samps_line[0][$ctype] ne "$type"){
				$flagMethy=1;
				$min=1000;
				next;
			}
		}
		$min=1000,$flagMethy=1,next if $min<$coverage;
	}
	### dmr data ###
	if($flagDMR==1){
		$dmrBeforeFlag+=1;
		$beforeDMR=@dmr;
		my $dmrLine=<DMR>;
		last if $dmrLine eq "";
		@dmr=split(/\t/,$dmrLine);

		### get overlap points ###
		if($dmrBeforeFlag!=1&&$dmrBefore[0] eq $dmr[0]){
			if($dmr[1]-$left<=$dmrBefore[2]+$right){	
				my $point=@{$methy_y[0]};
				#print "overlap region\t$point\n";
				my $index;
				for(my $ip=0;$ip<$point;$ip++){
					if($methy_y[0][$ip][0]>=$dmr[1]-$left){
						$index=$ip;
						last;
					}
				}
				for(my $k=0;$k<$sam_num;$k++){
					@{$methy_y[$k]}=@{$methy_y[$k]}[$index..$point-1];
				}
			}else{
				@methy_y=();
			}
		}else{
			@methy_y=();
		}
		@dmrBefore=@dmr;
	}

	if($dmr[0] eq $samps_line[0][$chromo]){
		if($dmr[1]-$left>$samps_line[0][$pos]){
			$flagMethy=1;
			$flagDMR=0;
			@methy_y=();
			$min=1000;
		}elsif($dmr[1]-$left<=$samps_line[0][$pos] && $dmr[2]+$right>=$samps_line[0][$pos]){
			$flagDMR=0;
			$flagMethy=1;
			for(my $k=0;$k<$sam_num;$k++){
				my $methylation_level=$samps_line[$k][$methy]/($samps_line[$k][$methy]+$samps_line[$k][$unmethy]);
				push @{$methy_y[$k]}, [$samps_line[$k][$pos],$methylation_level];
			}
			$min=1000;
		}elsif($dmr[2]+$right<$samps_line[0][$pos]){
			#### draw map ####
			my $position1=$dmr[1];
			my $position2=$dmr[2];
			my $positionString="$dmr[0]:$position1-$position2";
			my ($name,$mapPosition,$number1,$number2)=&swDMRsMap::PositionNormalized($positionString,$left,$right);	
			print DMRmapList "$dir/$name.R\n";
			my @methy_y_tmp;
			my @methy_y_x_tmp;
			my $n=@{$methy_y[0]};
			for(my $i=0;$i<$n;$i++){
				if($methy_y[0][$i][1] != -1){
					for(my $k=0;$k<$sam_num;$k++){
						push @{$methy_y_tmp[$k]},$methy_y[$k][$i][1];
					}
				}
				push @methy_y_x_tmp,$methy_y[0][$i][0];
			}
			###	Draw DMR map
			&swDMRsMap::DrawDMRmap(\@methy_y_x_tmp,\@methy_y_tmp,\@color,\@samName,$name,$mapPosition,$dir,$dmrFile2,$gene,$CGI,$dmr[0],$number1,$number2,$sam_num,$Rbin);
			$flagDMR=1;
			$flagMethy=0;
		}
		$tmpChr=$samps_line[0][$chromo];
	}else{
		$flagDMR=0;
		$flagMethy=1;
		if($tmpChr ne $samps_line[0][$chromo]){
		#### draw map ####
			my $position1=$dmr[1];
			my $position2=$dmr[2];
			my $positionString="$dmr[0]:$position1-$position2";
			my ($name,$mapPosition,$number1,$number2)=&swDMRsMap::PositionNormalized($positionString,$left,$right);
			print DMRmapList "$dir/$name.R\n";
       		my @methy_y_tmp;
	        my @methy_y_x_tmp;
    	    my $n=@{$methy_y[0]};
        	for(my $i=0;$i<$n;$i++){
        		if($methy_y[0][$i][1] != -1){
            		for(my $k=0;$k<$sam_num;$k++){
                		push @{$methy_y_tmp[$k]},$methy_y[$k][$i][1];
	        		}
    	        }
       		 	push @methy_y_x_tmp,$methy_y[0][$i][0];
      		}
        	### Draw DMR map
	        &swDMRsMap::DrawDMRmap(\@methy_y_x_tmp,\@methy_y_tmp,\@color,\@samName,$name,$mapPosition,$dir,$dmrFile2,$gene,$CGI,$dmr[0],$number1,$number2,$sam_num,$Rbin);
    	    $flagMethy=0;
			$flagDMR=1;
			$tmpChr=$samps_line[0][$chromo];
		}
		@methy_y=();
		$min=1000;
	}
}
close(DMRmapList);
#`perl $Script/DMRmap.R.Plot.parallel.pl $dir/DMRmap.R.list $Rbin`;
