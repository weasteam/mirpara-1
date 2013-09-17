#! /usr/bin/perl
#	by Yonggan Wu (weasteam@gmail.com)
#	by Simon Rayner (raynere@wh.iov.cn)
#################################USE############################################
use warnings;
#use strict;
use Cwd qw(abs_path);
use Cwd;
use File::chdir;
#################################USE############################################
my $usage="perl $0 <syspmt> <seq title> <sequence>\n\n";
my $syspmt=shift or die $usage;
my $title=shift or die $usage;
my $seq=shift or die $usage;
our %syspmt;
syspmt($syspmt);
#update file directory
my $tmp=$syspmt;
$tmp=~s/\//,/g;
my @tmp=split(",",$tmp);
$dir=$syspmt;
$tmp=pop @tmp;
$dir=~s/\/$tmp//;
$CWD = $dir;
#################################START PROGRAM##################################
#the script will optimize the input sequences and collect the required informaiton
print "Optimizing sequences for $title...\n";
my @strdata=dofold($seq);#do UNAfolding
@tmp=split("_",$title);
my $subtitle=$tmp[0];
my %validatedseq=();
if (@strdata ne 0){#if the secondary structure could be generated
	my $line=1;
	my @tmp=();
	foreach (@strdata){
		next if /^Sequence/;
		next if /^Bases/;
		next if /dG/;
		$_=~s/\n//;
		$tmp[$line]=$_;
		if ($line eq 7){
			$line=0;
			#generate the secondary structure
			#the folding structure multiple check
			#################################################################
			my $tmp=hairpin2seq(cutstr(@tmp[3..6]));
			uniquestr:
			if (length($tmp)>=$syspmt{'prilengthlimit'}){#length filter
				#my @strdata2=dofold($tmp);
				my ($mfe,$priseq,@str)=longeststrtrim(dofold($tmp));
				if (length($priseq) ne length($tmp)){#the hairpin contain budding stem
					$tmp=$priseq;
					goto uniquestr;
				}
				##############################################################
				#generate the validated seq
				if (not exists $validatedseq{"$subtitle\_$priseq"}){
					@tmp=split("_",$title);
					$tmp[1]+=index($seq,$priseq);
					$tmp[2]=$tmp[1]+length($priseq)-1;
					$validatedseq{"$subtitle\_$priseq"}=1;
					#print "$subtitle\_$priseq\t$subtitle\:$tmp[1]-$tmp[2],$mfe,$str[0],$str[1],$str[2],$str[3]\n";
					$str[0]=~s/ /\_/g;
					$str[1]=~s/ /\_/g;
					$str[2]=~s/ /\_/g;
					$str[3]=~s/ /\_/g;
					system "perl $syspmt{'share'}/miRPara_find_candidates.pl $syspmt $subtitle\:$tmp[1]-$tmp[2],$priseq,$mfe,$str[0],$str[1],$str[2],$str[3]";
				}
			}
		}
		$line+=1;
	}
}
sub dofold{
   #	to generate the secondary structure
   #1	received a title and a cuted ~500 seq;
   #2	do unafuld and return a str
   my ($seq)=@_;
   my $rand=int(rand(9999)).int(rand(9999));
   my $fh="fold".$rand;
   open ($fh,">$rand");
   print {$fh} ">$rand\n$seq\n";
   close $fh;
   system "$syspmt{'unafold'} $rand > $rand.log";
   my @str="";
   if (-e "$rand\_1.ct"){
	  system "$syspmt{'ct2out'} <$rand\_1.ct> $rand.str";#run ct2out
	  open ($fh,"$rand.str");
	  @str=<$fh>;
	  close $fh;
   }
	system "rm $rand*";
	return @str;
}
sub hairpin2seq {#conver the second structure into a line
   #The hairpin2seq was used to generate a seq from second str
   #	receive 4 str
   #	return one seq
   my (@s)=@_;
   my $seq;#for the sequence
   my $upperseq;#for the upper strand
   my $lowerseq;#for the lower strand
   my $m=0;
   while ($m<=(length($s[1])-1)){
	  $upperseq.=substr($s[0],$m,1).substr($s[1],$m,1);#get the seq
	  $lowerseq.=substr($s[2],-($m+1),1).substr($s[3],-($m+1),1);
	  $m+=1;
   }
   $seq=$upperseq.$lowerseq;
   $seq=~ s/[\\\-\.]//g;#get rid of \\
   $seq=~ s/\s+//g;#get rid of black
   return $seq;
}
sub cutstr{
   #generate the proper strs, without budding stem;
   #receive 4 str
   #return 4 str
   #hairpin length limit
   my $hairpinlength=60;
   my (@seq)=@_;
   while ($seq[3]=~/\\/){
	  $seq[3]=substr($seq[3],index($seq[3],"\\")+1);
   }
   $seq[3]=~s/^\s+//;
   if (length($seq[3])>$hairpinlength){
		$seq[3]=substr($seq[3],-$hairpinlength);
   }
   $seq[0]=substr($seq[0],-length($seq[3]));
   $seq[1]=substr($seq[1],-length($seq[3]));
   $seq[2]=substr($seq[2],-length($seq[3]));
   return @seq;
}
sub longeststrtrim{
   #generate the longest hairpin from the given str
   #receive unafold random number
   #return 4 str and $mfe, and primary sequences
   my (@strdata)=@_;
   my $mfe;
   my $line=1;
   my @tmp;
   my @seq;
   my $tmp;
   my %str;
   my %pri;
   foreach (@strdata){
	  next if /^Sequence/;
	  next if /^Bases/;
	  $_=~s/\n//;
	  if ($_=~/dG/){
		 $mfe=substr($_,index($_,"=")+1);
		 $mfe=~s/^\s+//g;#delete the blacks
		 $mfe=substr($mfe,0,index($mfe,"\t"));
		 next;
	  }
	  $tmp[$line]=$_;
	  if ($line eq 7){
		 $line=0;
		 @seq=cutstr(@tmp[3..6]);
		 $tmp=hairpin2seq(@seq);
		 $str{length($tmp)}=join(",",@seq);
		 $pri{length($tmp)}=$tmp;
	  }
	  $line+=1
   }
   @tmp=sort { $b <=> $a;} (keys(%str));
   @seq=split(",",$str{$tmp[0]});
   return ($mfe,$pri{$tmp[0]},@seq);
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
