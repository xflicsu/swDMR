#!/usr/local/bin/perl
package Annotation;
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/module/lib64/perl5/site_perl";
use lib "$FindBin::Bin/module/lib64/perl5";
use lib "$FindBin::Bin/module/lib/perl5";
use lib "$FindBin::Bin/module/lib/perl5/site_perl";
use lib "$FindBin::Bin/module/lib/perl";
use lib "$FindBin::Bin/module/share/perl";
use File::Basename;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(TestFormat );

sub TestFormat{
	my $file=shift;
#	print "55555\t$$file\n";
	open IN,"$$file"or die"$!";
	my $head=<IN>;
	chomp($head);
	my @line=split(/\s+/,$head);
	my $format;
	my $FeatureIndex=0;
	if($line[0]=~/chr\S+/i&&$line[1]=~/\d+/&&$line[2]=~/\d+/){
		$format="bed";
		$FeatureIndex=3;
	}elsif($line[0]=~/chr\S+/i&&$line[3]=~/\S+/&&$line[4]=~/\d+/){
		$format="gffOrgtf";
		$FeatureIndex=2;
	}else{
		$format="resistance";
	};
	return($format,$FeatureIndex);
}
1;

