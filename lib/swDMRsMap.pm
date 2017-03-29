#!/usr/local/bin/perl -w
package swDMRsMap;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/module/lib64/perl5/site_perl";
use lib "$FindBin::Bin/module/lib64/perl5";
use lib "$FindBin::Bin/module/lib/perl5";
use lib "$FindBin::Bin/module/lib/perl5/site_perl";
use lib "$FindBin::Bin/module/lib/perl";
use lib "$FindBin::Bin/module/share/perl";
use File::Basename;
use PerlIO::gzip;
my $Rlib="$FindBin::Bin/../lib/R-packages";
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(PositionNormalized DrawDMRmap);
sub PositionNormalized{
	my $position=shift;
	my $left=shift;
	my $right=shift;
	$position=~/(\S+?):(\d+?)-(\d+?)$/i;
	my $chr=$1;my $number1=$2;my $number2=$3;
	chomp($chr);chomp($number1);chomp($number2);
	my $namePosition1=Normalize($number1);
	my $namePosition2=Normalize($number2);
	my $nameNormalize="$chr:$namePosition1-$namePosition2";
	if($number1>$left){
		$number1-=$left;
	}else{
		$number1=1;
	}
	$number2+=$right;
	my $position1=Normalize($number1);
	my $position2=Normalize($number2);
	my $NormalizePosition="$chr:$position1-$position2";
	return ($nameNormalize,$NormalizePosition,$number1,$number2);
}

sub Normalize{
	my $number=shift;
	my $nlen=rindex $number."\$", "\$";
	my $fNumber;
	while ( $nlen>3 ) { 
		$fNumber = ",". substr($number,$nlen-3,3).$fNumber;
		$number = substr($number,0,-3);
		$nlen=rindex $number."\$", "\$";
	}
	if ( $nlen <= 3 ) {
		$fNumber = $number.$fNumber;
	}
	return $fNumber;
}

sub DrawDMRmap{

	my ($methy_y_x_tmp,$methy_y_tmp,$color_tmp,$sam_name_tmp,$name,$title,$outdir,$dmr,$gene,$CGI,$chr,$start,$end,$sam_num,$Rbin)=@_;
	
	
	####data get
	my @methy_y_x=@{$methy_y_x_tmp};
	my @methy_y=@{$methy_y_tmp};
	my @color=@{$color_tmp};
	my @sam_name=@{$sam_name_tmp};
	my $number=@sam_name;
	my @ownColor=@color[0..$number-1];
	#print join "\t",@{$methy_y[0]},"\n";
	my $samName=join "\",\"",@sam_name;
	my $colorName=join "\",\"",@ownColor;
	my $methyX=join ",",@methy_y_x;
	open OUT,">$outdir/$name.R"||die"$!";
	#######################################################################
	print OUT "
library(\"pspline\",lib.loc=\"$Rlib\")
pdf(\"$outdir/$name.pdf\",width=10,height=4)
par(mar=c(0.2,5,0,1),oma=c(3,1,1,1))
m=matrix(c(1:2),2,1)
layout(m,heights=c(2,1))
";
	my $methyY;
	for(my $k=0;$k<$sam_num;$k++){
		$methyY=join ",",@{$methy_y[$k]};
		print OUT "
plot(c($methyX),c($methyY)*100,pch=20,col=\"$color[$k]\",ylab=\"Methylation level (%)\",xlim=c($start,$end),ylim=c(0,100),xlab=\"\",xaxt=\"n\")
lines(smooth.Pspline(c($methyX),c($methyY)*100,method=2),col=\"$color[$k]\",lwd=4)
";
		if($k!=$sam_num-1){
			print OUT "
par(new=T)\n";
		}
	}
	
	print OUT "legend(\"topright\",c(\"$samName\"),pch=16,col=c(\"$colorName\"))
box(col=\"gray\")\n";
	print OUT "plot(c($methyX),c($methyY)*100,pch=20,type=\"n\",xlim=c($start,$end),ylim=c(0,100),ylab=\"\",yaxt=\"n\")\n";
	#### dmr tick
	open DMR,"$dmr"||die "$!";
	while(<DMR>){
		chomp;
		my @dmrs=split;
		next if $dmrs[0] ne "$chr";
		next if $dmrs[2] <$start;
		last if $dmrs[1] >$end;
		if($dmrs[1]>=$start){
			print OUT "segments($dmrs[1],90,$dmrs[2],90,col=\"#D52B2B\",lwd=5)\n";
		}
	}

	#### CGI tick
	open CGI, "$CGI"||die"$!";
	while(<CGI>){
		chomp;
		my @cgi=split;
		my $chrIndex=0;
		my $startIndex=1;
		my $endIndex=2;
		if(!$cgi[0]=~/^chr/i){
			$chrIndex+=1;
			$startIndex+=1;
			$endIndex+=1;
		}
		next if $cgi[$chrIndex] ne $chr;
		next if $cgi[$startIndex] < $start;
		last if $cgi[$endIndex] > $end;
		if($cgi[1]>=$start){
			print OUT "segments($cgi[$startIndex],70,$cgi[$endIndex],70,col=\"#44BB8C\",lwd=5)\n";
		}
	}
	close(CGI);
	#### gene tick 
	open GENE,"$gene"||die"$!";
	my $first=0;
	my $tmp=0;
	while(<GENE>){
		chomp;
		my @gene=split;
		next if $gene[0] ne "$chr";
		next if $gene[2] < $start;
		last if $gene[1] > $end;
		if($first!=0){
			$tmp-=8;
		}
		my $start_tmp=$start;my $end_tmp=$end;
		if($gene[1]>=$start){
			$start_tmp=$gene[1];
		}
		if($gene[2]<=$end){
			$end_tmp=$gene[2];
		}
		my $middle=($end_tmp-$start_tmp)/2+$start_tmp;
		my $piece=($end_tmp-$start_tmp)/40;
		
		### with arrows
		if($gene[3] eq "+"){
			print OUT "segments($gene[1],50+$tmp,$gene[2],50+$tmp,col=\"gray\",lwd=5)
for(i in 0:39){
	arrows($start_tmp,50+$tmp,$end_tmp-i*$piece,50+$tmp,col=\"#27408B\",length=0.06)
}
text($middle,50+$tmp,labels=\"$gene[5]\")
";
		}else{
			print OUT "segments($gene[1],50+$tmp,$gene[2],50+$tmp,col=\"gray\",lwd=5)
for(i in 0:39){
	arrows($end_tmp,50+$tmp,$start_tmp+i*$piece,50+$tmp,col=\"#27408B\",length=0.06)
}
text($middle,50+$tmp,labels=\"$gene[5]\")
";
		}
		$first+=1;
	}
	print OUT "laby=c(\"DMR\",\"CGI\",\"Gene\")
axis(2,at=c(90,70,50),labels=laby,cex=2,cex.lab=2,font=2,las=2)
";
	print OUT "box(col=\"gray\")\n";
	close(OUT);
	close(GENE);
#	`$Rbin CMD BATCH $outdir/$name.R;rm $name.Rout`;
}
1;
