#! /usr/bin/perl
#	by Yonggan Wu (weasteam@gmail.com)
#	by Simon Rayner (raynere@wh.iov.cn)
#################################USE############################################
use warnings;
use strict;
#################################USE############################################
my $usage="perl $0 <syspmt> <seq title> <long sequence>\n\n";
my $syspmt=shift or die $usage;
my $title=shift or die $usage;
my $seq=shift or die $usage;
our %syspmt;
syspmt($syspmt);
#################################START PROGRAM##################################
#1	receive one input sequences and its name
#2	split it according to the length
print "spliting seqence for $title...\n";
my @pstart;#the array of start positions
my @pend;#$the array of end positions
my $len=length($seq);
my $splitseq;
my $priseq;
my @tmp;
#generate the start and end positions
$pstart[0]=1;
$pend[0]=$syspmt{'foldlen'};
my $tmp=1;
while ($pend[$tmp-1]<=$len){
	$pstart[$tmp]=$pstart[$tmp-1]+($syspmt{'foldlen'}-$syspmt{'overlaplen'});
	$pend[$tmp]=$pend[$tmp-1]+($syspmt{'foldlen'}-$syspmt{'overlaplen'});
	$tmp+=1;
}
$tmp-=1;
$pstart[$tmp]=$pstart[$tmp-1]+($syspmt{'foldlen'}-$syspmt{'overlaplen'});
$pend[$tmp]=$len;
if ($len<=$syspmt{'foldlen'}){
	@pstart=();
	@pend=();
	$pstart[0]=1;
	$pend[0]=$len;
}
for (my $i=1;$i<=@pstart;$i+=1){
	$splitseq=substr($seq,$pstart[$i-1]-1,$pend[$i-1]-$pstart[$i-1]+1);
	#print "$title\_$pstart[$i-1]\_$pend[$i-1]\n$splitseq\n";
	#optimize the sequence from here
	system "perl $syspmt{'share'}/miRPara_optimize_seq.pl $syspmt $title\_$pstart[$i-1]\_$pend[$i-1] $splitseq";
}
sub syspmt{
	my ($syspmt)=@_;
	open (IN,"$syspmt") or die "Error: not fount $syspmt\n";
	while (<IN>){
		next if /^#/;
		$_=~s/[\n\r]//g;
		my @tmp=split("\t",$_);
		$syspmt{$tmp[0]}=$tmp[1];
	}
	close IN;
}
