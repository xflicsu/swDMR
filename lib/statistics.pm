#!/usr/local/bin/perl
package statistics;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(T_test Kolmogorov Fisher ChiSquare Wilcoxon ANOVA Kruskal XX make_mean_rate);
use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/module/lib64/perl5/site_perl";
use lib "$FindBin::Bin/module/lib64/perl5";
use lib "$FindBin::Bin/module/lib/perl5";
use lib "$FindBin::Bin/module/lib/perl5/site_perl";
use lib "$FindBin::Bin/module/lib/perl";
use lib "$FindBin::Bin/module/share/perl";
use Statistics::R;

=pod
=head1 Function
   Create all the function used in cwDMR program;
=head1 Input
   Data inputed in the sub function should be like the following format
        : @S1MethyReads @S1UnmethyReads @S2MethyReads @S2UnmethyReads @S3MethyReads @S3UnmethyReads
   Statistics methods:
	Two samples statistics methods: T_test, Kolmogorov, Fisher, Wilcox, ChiSquare.
	Miltiple samples statisitcs methods: Kruskal, ANOVA.
=head1 Problem
	ANOVA and Kolmogorov
=cut

our $VERSION = '0.1';

#*********T Test********#
sub T_test{
	my ($methyTmp) = shift;
	my $num=@{$methyTmp};
	my @methy=@{$methyTmp};
	my @point_rate;
	for(my $i=0;$i<$num;$i+=2){
		my @rate=PointMethyRate(\@{$methy[$i]},\@{$methy[$i+1]});
		push @point_rate,[@rate];
	}
	my $R = Statistics::R->new();
	$R -> set('x',[@{$point_rate[0]}]);
	$R -> set('y',[@{$point_rate[1]}]);
	my $cmds=<<"EOF";
	result=t.test(x,y)\$p.value
EOF
	my $array=$R -> get('x');
	$R -> run($cmds);
	my $pVal = $R -> get('result');
	$R->stop();
	return ($pVal);
}
sub XX{
	my $a=shift;
	print "oooooooooooo\t$a\n";
	print "1234\tit is ok!\n";
}
#*********Two-sample Kolmogorov-Smirnov test********#
sub Kolmogorov{
	my ($methyTmp) = shift;
	my $num=@{$methyTmp};
	my @methy=@{$methyTmp};
	my @point_rate;
	for(my $i=0;$i<$num;$i+=2){
		my @rate=PointMethyRate(\@{$methy[$i]},\@{$methy[$i+1]});
		push @point_rate,[@rate];
	}
	my $R = Statistics::R->new();
	$R -> set('x',[@{$point_rate[0]}]);
	$R -> set('y',[@{$point_rate[1]}]);
	my $cmds=<<"EOF";
		result=ks.test(x,y)\$p.value
EOF
	$R -> run($cmds);
	my $pVal = $R -> get('result');
	$R->stop();
	print $pVal,"\n";
	return ($pVal);
}

#*********Wilcoxon Test********#
sub Wilcoxon{
	my ($methyTmp) = shift;
	my $num=@{$methyTmp};
	my @methy=@{$methyTmp};
	my @point_rate;
	for(my $i=0;$i<$num;$i+=2){
		my @rate=PointMethyRate(\@{$methy[$i]},\@{$methy[$i+1]});
		push @point_rate,[@rate];
	}
	my $R = Statistics::R->new();
	$R -> set('x',[@{$point_rate[0]}]);
	$R -> set('y',[@{$point_rate[1]}]);
	my $cmds=<<"EOF";
		result=wilcox.test(x,y)\$p.value
EOF
	$R -> run($cmds);
	my $pVal = $R -> get('result');
	$R->stop();
	return ($pVal);
}
#*********ChiSquare Test********#
#### 2 X 2 table P140 Medical Statistics (Second Edition) Sun zhenqiu

sub ChiSquare{
	my ($methyTmp) = shift;
	my $num=@{$methyTmp};
	my @methy=@{$methyTmp};
	my @methy_table;
	my @un_methy_table;
	for(my $i=0;$i<$num;$i+=2){
		my $methy=sum(\@{$methy[$i]});
		my $unmethy=sum(\@{$methy[$i+1]});
		push @methy_table,$methy;
		push @un_methy_table,$unmethy;
	}
	my @MethyTable=(@methy_table,@un_methy_table);
	my $R = Statistics::R->new();
	$R -> set('x',[@MethyTable]);
	my $cmds=<<"EOF";
		result=chisq.test(matrix(x,nrow=2),simulate.p.value=TRUE)\$p.value
EOF
	$R -> run($cmds);
	my $pVal = $R -> get('result');
	$R->stop();
	my $sum_RC=sum(\@MethyTable);
	my $methyRate=sum(\@methy_table)/$sum_RC;
	my $unMethyRate=1-$methyRate;
	my $sam1=$methy_table[0]+$un_methy_table[0];
	my $sam2=$methy_table[1]+$un_methy_table[1];
	my @Ttable;
	foreach($sam1,$sam2){
		push @Ttable,$_*$methyRate;
		push @Ttable,$_*$unMethyRate;
	}
	my $min_T=min(\@Ttable);
	if(($sum_RC>=40&&$min_T>=1&&$min_T<5)||($sum_RC<40||$min_T<1)){
		$R = Statistics::R->new();
		$R -> set('x',[@MethyTable]);
		$cmds=<<"EOF";
			result=fisher.test(matrix(x,nrow=2),simulate.p.value=TRUE)\$p.value
EOF
		$R -> run($cmds);
		$pVal = $R -> get('result');
		$R->stop();
	}
	return ($pVal);
}
#*********Fisher's Exact Test********#
#### 2 X 2 table P140 Medical Statistics (Second Edition) Sun zhenqiu
sub Fisher{
	my ($methyTmp) = shift;
	my $num=@{$methyTmp};
	my @methy=@{$methyTmp};
	#my @MethyTable = RC_table(\@methy);
	my @methy_table;
	my @un_methy_table;
	for(my $i=0;$i<$num;$i+=2){
		my $methy=sum(\@{$methy[$i]});
		my $unmethy=sum(\@{$methy[$i+1]});
		push @methy_table,$methy;
		push @un_methy_table,$unmethy;
	}
	my @MethyTable=(@methy_table,@un_methy_table);
	#print join "\t55555555\t\n",@MethyTable;
	#print "44444444444444444444\n";
	my $R = Statistics::R->new();
	$R -> set('x',[@MethyTable]);
	#print "33333333\t@MethyTable\n";
	my $output_value = $R->get('x');
	my $cmds=<<"EOF";
		result=fisher.test(matrix(x,nrow=2),simulate.p.value=TRUE)\$p.value
EOF
	#print join"\n", @{$R->get('x')};
	$R -> run($cmds);
	my $pVal = $R -> get('result');
	$R->stop();
	return ($pVal);
}
#### compare groups of data through ANOVA
sub ANOVA{
	my ($methyTmp) = shift;
	my $num=@{$methyTmp};
	my @methy=@{$methyTmp};
	my @point_rate;
	for(my $i=0;$i<$num;$i+=2){
		my @rate=PointMethyRate(\@{$methy[$i]},\@{$methy[$i+1]});
		push @point_rate,[@rate];
	}
	my $R = Statistics::R->new();
	my $cmds1="class=c(";
	my $cmds2="rate=c(";
	for(my $i=0;$i<$num/2;$i++){
		$R -> set("a$i",[@{$point_rate[$i]}]);
		$cmds1.="a$i";
		$cmds2.="rep(\"a$i\",length(a$i))";
		$cmds1.=",",$cmds2.="," if $i<$num/2-1;
	}
	$cmds1.=");";
	$cmds2.=");";
	my $cmds3="com=cbind(class,rate);class=factor(com[,2]);rate=as.numeric(com[,1]);result=anova(lm(rate~class))";
	$cmds=join"\n",$cmds1,$cmds2,$cmds3;
	$R -> run($cmds);
	my $pVal = $R -> get('result');
	$R->stop();
	if(@$pVal[19] eq "<"){
		@$pVal[19]=@$pVal[20];
	}
	return (@$pVal[19]);
}

#*********Friedman test********#
sub Friedman{
	my ($samp1,$samp2)=@_;
	my @methy = @_;
	my $R = Statistics::R->new();
	$R -> set('x',[@$samp1]);
	$R -> set('y',[@$samp2]);
	my $cmds=<<"EOF";
		result=t.test(x,y)
EOF
}
sub Kruskal{
	my ($methyTmp) = shift;
	my $num=@{$methyTmp};
	my @methy=@{$methyTmp};
	my @point_rate;
	for(my $i=0;$i<$num;$i+=2){
		my @rate=PointMethyRate(\@{$methy[$i]},\@{$methy[$i+1]});
		push @point_rate,[@rate];
	}
	my $R = Statistics::R->new();
	my $cmds="result <- kruskal.test(list(";
	for(my $i=0;$i<$num/2;$i++){
		$R -> set("a$i",[@{$point_rate[$i]}]);
		$cmds.="a$i";
		$cmds.="," if $i<$num/2-1;
	}
	$cmds.="))\$p.value";
	$R -> run($cmds);
	my $pVal = $R -> get('result');
	$R->stop();
	return ($pVal);
}

#********function used in statistics methods********#
##export make_mean_rate##

sub make_mean_rate{
	my ($methyTmp) = shift;
	my $num=@{$methyTmp};
#	print "33333333333\t$num\n";
	my @methy=@{$methyTmp};
	my @mean_rate;
	my @mean_depth;
	for(my $i=0;$i<$num;$i+=2){
		my $methy=sum(\@{$methy[$i]});
		my $point=@{$methy[$i]};
		my $unmethy=sum(\@{$methy[$i+1]});
		my $rate=$methy/($methy+$unmethy);
		my $depth=($methy+$unmethy)/$point;
		push @mean_rate,$rate;
		push @mean_depth,$depth;
	}
	my ($max,$min)=max_min(\@mean_rate);
	return (\@mean_rate,$max,$min,\@mean_depth);
}

#sub make_methy_rate{
#	my @methy = @_;
#	my $num = @methy;
#	my @point_rate;
#	for(my $i=0;$i<$num;$i+=2){
#		my @rate=PointMethyRate(@{$methy[$i]},@{$methy[$i+1]});
#		push @point_rate,@rate;
#	}
#	return (@point_rate);
#}
sub sum{
	my $array = shift;
	my $sum=0;
	foreach(@{$array}){
		$sum+=$_;
	}
	return($sum);
}
##export PointMethyRate##
sub PointMethyRate{
	my ($methy,$unmethy) = @_;
	my $num = @{$methy};
	my @point_rate;
	for(my $i=0;$i<$num;$i++){
#		print "44444444\t@{$methy}[$i]\t@{$unmethy}[$i]\n";
		my $rate=@{$methy}[$i]/(@{$methy}[$i]+@{$unmethy}[$i]);
		push @point_rate,$rate;
	}
	return (@point_rate);
}
sub max_min{
	my $array = shift;
	my $min=1;
	my $max=0;
	map{if($min>=$_){$min=$_;}} (@{$array});
	map{if($max<=$_){$max=$_;}} (@{$array});
	return ($max,$min);
}
sub min{
	my $array = shift;
	my @myArray=@{$array};
	my $min=$myArray[0];
	map{if($min>=$_){$min=$_;}} (@myArray);
}
1;
