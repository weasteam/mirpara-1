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
my $usage="perl $0 <syspmt> <data>\n\n";
my $syspmt=shift or die $usage;
my $data=shift or die $usage;
$data=~s/\_/ /g;
our %syspmt;
syspmt($syspmt);
our %lower=prilimit("lower");#prilower band
our %upper=prilimit("upper");#priupper band
#update file directory
my $tmp=$syspmt;
$tmp=~s/\//,/g;
my @tmp=split(",",$tmp);
$dir=$syspmt;
$tmp=pop @tmp;
$dir=~s/\/$tmp//;
$CWD = $dir;
#################################START PROGRAM##################################
#	to generate the posssible candidates
#	receive title $mfe,$priseq,@str
#	process for parameters;
my %para;
my @milen;
my $rand;
@tmp=split(",",$data);
my @dat;
my $dat;
my $title=$tmp[0];
my $seq=$tmp[1];
my $mfe=$tmp[2];
my @str=($tmp[3],$tmp[4],$tmp[5],$tmp[6]);
@tmp=();
#priid,priseq,prilength,budding stem
$para{'priid'}=$title;#_, incase same to the real name.
$para{'priseq'}=$seq;
$para{'prilength'}=length($para{'priseq'});
$para{'buddingstem'}=0;#0 refer to no
$para{'primfe'}=$mfe;
($para{'pristr_1'},$para{'pristr_2'},$para{'pristr_3'},$para{'pristr_4'})=@str;
#create candidates 20-24
open (PMT,">>$syspmt{'outfilepmt'}");
open (OUT,">>$syspmt{'outfilerst'}");
for ($para{'mistart'}=1;$para{'mistart'}<=$para{'prilength'};$para{'mistart'}+=1){
  @milen=(20..24);
  foreach (@milen){
	 #milength,miend
	 $para{'milength'}=$_;
	 #----------------------------------------------------------------------
	 #debug
	 #$para{'mistart'}=17;
	 #$para{'milength'}=21;
	 #----------------------------------------------------------------------
	 $para{'miend'}=$para{'mistart'}+$para{'milength'}-1;
	 if ($para{'miend'}<=$para{'prilength'}){#incase out of range
		#miseq
		$para{'miseq'}=uc(substr($para{'priseq'},$para{'mistart'}-1,$para{'milength'}));
		#miid
		$para{'miid'}="$para{'priid'}\:$para{'mistart'}\_$para{'milength'}";
		#print "$para{'miid'}\n";
		#priseq
		$para{'priseq'}=lc(substr($para{'priseq'},0,$para{'mistart'}-1)).uc(substr($para{'priseq'},$para{'mistart'}-1,$para{'milength'})).lc(substr($para{'priseq'},$para{'miend'}));
		#basal segment
		$tmp=$para{'pristr_1'};
		if ($tmp=~/^\s/){
		   $para{'length_basalsegment'}=0;
		}
		else{
		   $tmp=~s/\s+/,/g;
		   @tmp=split(",",$tmp);
		   $para{'length_basalsegment'}=length($tmp[0]);
		}
		#strand
		$tmp=join("",$para{'pristr_1'},$para{'pristr_2'});
		$tmp=~s/[\s\\\.-]+//g;
		if ($para{'miend'} <= length($tmp)){
		   $para{'strand'}='5';
		}
		elsif($para{'mistart'} > length($tmp)){
		   $para{'strand'}='3';
		}
		else{
		   $para{'strand'}=0;#across the terminal loop
		}
		if ($para{'strand'} ne 0){
		   #upperstart
		   @tmp=();
		   my $primary;
		   if ($para{'strand'} eq 5){
			  $primary=$para{'mistart'}-1;
			  @tmp[1,2]=($para{'pristr_1'},$para{'pristr_2'});
		   }
		   else{
			  $primary=$para{'prilength'}-$para{'miend'};
			  @tmp[1,2]=($para{'pristr_3'},$para{'pristr_4'});
		   }
		   $para{'upperstart'}=$primary;
		   $tmp=$para{'pristr_1'};
		   until (length($tmp) eq $primary){
			  $tmp=substr($tmp[1],0,$para{'upperstart'}).substr($tmp[2],0,$para{'upperstart'});
			  $tmp=~s/[\-\s]+//g;
			  $para{'upperstart'}+=1;
		   }
		   #upperend
		   if ($para{'strand'} eq 5){
			  $primary=$para{'miend'};
			  @tmp[1,2]=($para{'pristr_1'},$para{'pristr_2'});
		   }
		   else{
			  $primary=$para{'prilength'}-$para{'mistart'}+1;
			  @tmp[1,2]=($para{'pristr_3'},$para{'pristr_4'});
		   }
		   $para{'upperend'}=$primary;
		   $tmp=$para{'pristr_1'};
		   until (length($tmp) eq $primary){
			  $tmp=substr($tmp[1],0,$para{'upperend'}).substr($tmp[2],0,$para{'upperend'});
			  $tmp=~s/[\-\s]+//g;
			  $para{'upperend'}+=1;
		   }
		   $para{'upperend'}-=1;
		   #the second structure of the pri-miRNA
		   if ($para{'strand'} eq 5){
			  $para{'pristr_3'}=lc($para{'pristr_3'});
		      $para{'pristr_4'}=lc($para{'pristr_4'});
		      $para{'pristr_1'}=lc(substr($para{'pristr_1'},0,$para{'upperstart'}-1)).uc(substr($para{'pristr_1'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_1'},$para{'upperend'}));
		      $para{'pristr_2'}=lc(substr($para{'pristr_2'},0,$para{'upperstart'}-1)).uc(substr($para{'pristr_2'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_2'},$para{'upperend'}));
		   }
		   else{
			  $para{'pristr_1'}=lc($para{'pristr_1'});
			  $para{'pristr_2'}=lc($para{'pristr_2'});
			  $para{'pristr_3'}=lc(substr($para{'pristr_3'},0,$para{'upperstart'}-1)).uc(substr($para{'pristr_3'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_3'},$para{'upperend'}));
			  $para{'pristr_4'}=lc(substr($para{'pristr_4'},0,$para{'upperstart'}-1)).uc(substr($para{'pristr_4'},$para{'upperstart'}-1,$para{'upperend'}-$para{'upperstart'}+1)).lc(substr($para{'pristr_4'},$para{'upperend'}));
		   }
		   #terminal loop
		   $tmp="";
		   @tmp=($para{'pristr_1'},$para{'pristr_2'},$para{'pristr_3'},$para{'pristr_4'});
		   foreach (@tmp){
			  $_=~s/\s+$//;
			  $_=~s/\s+/,/g;
			  @dat=split(",",$_);
			  $tmp.=pop @dat;
		   }
		   $tmp=~s/\\//;
		   $para{'length_terminalloop'}=length($tmp);
		   #pre
		   @tmp=($para{'pristr_1'},$para{'pristr_2'},$para{'pristr_3'},$para{'pristr_4'});
		   $tmp=length($para{'pristr_1'});
		   foreach (@tmp){
			  $_=~s/[AUCGT]/N/g;
			  if (index($_,"N")<$tmp and index($_,"N") ne -1){
				 $tmp=index($_,"N");
			  }
		   }
		   #prestr
		   $para{'prestr_1'}=substr($para{'pristr_1'},$tmp);
		   $para{'prestr_2'}=substr($para{'pristr_2'},$tmp);
		   $para{'prestr_3'}=substr($para{'pristr_3'},$tmp);
		   $para{'prestr_4'}=substr($para{'pristr_4'},$tmp);
		   #pre seq
		   $para{'preseq'}=hairpin2seq($para{'prestr_1'},$para{'prestr_2'},$para{'prestr_3'},$para{'prestr_4'});
		   $para{'prelength'}=length($para{'preseq'});
		   #premfe
		   $para{'premfe'}=$para{'primfe'}*($para{'prelength'}/$para{'prilength'});
		   #lower str
		   $dat=$tmp-$para{'length_basalsegment'};
		   if ($dat <0){
			  $dat=0
		   }
		   $para{'lowerstr_1'}=substr($para{'pristr_1'},$para{'length_basalsegment'},$dat);
		   $para{'lowerstr_2'}=substr($para{'pristr_2'},$para{'length_basalsegment'},$dat);
		   $para{'lowerstr_3'}=substr($para{'pristr_3'},$para{'length_basalsegment'},$dat);
		   $para{'lowerstr_4'}=substr($para{'pristr_4'},$para{'length_basalsegment'},$dat);
		   #mi
		   @tmp=($para{'prestr_1'},$para{'prestr_2'},$para{'prestr_3'},$para{'prestr_4'});
		   $tmp=length($para{'prestr_1'});
		   foreach (@tmp){
			  $_=~s/[atcgu]/n/g;
			  $_=~s/[AUCGT]/N/g;
			  if (index($_,"n")<$tmp and $_=~/N/){
				 $tmp=index($_,"n");
			  }
		   }
		   #mistr
		   $para{'mistr_1'}=substr($para{'prestr_1'},0,$tmp);
		   $para{'mistr_2'}=substr($para{'prestr_2'},0,$tmp);
		   $para{'mistr_3'}=substr($para{'prestr_3'},0,$tmp);
		   $para{'mistr_4'}=substr($para{'prestr_4'},0,$tmp);
		   #top str
		   $dat=length($para{'prestr_1'})-length($para{'mistr_1'})-int($para{'length_terminalloop'}/2+0.1);
		   $para{'topstr_1'}=substr($para{'prestr_1'},$tmp,$dat);
		   $para{'topstr_2'}=substr($para{'prestr_2'},$tmp,$dat);
		   $para{'topstr_3'}=substr($para{'prestr_3'},$tmp,$dat);
		   $para{'topstr_4'}=substr($para{'prestr_4'},$tmp,$dat);
		   #str length
		   $para{'length_upperstem'}=length($para{'mistr_2'});
		   $para{'length_topstem'}=length($para{'topstr_2'});
		   $para{'length_lowerstem'}=length($para{'lowerstr_2'});
		   #gc content
		   $para{'prigc'}=gc($para{'priseq'});#prigc
		   $para{'pregc'}=gc($para{'preseq'});#pregc
		   $para{'migc'}=gc($para{'miseq'});#migc
		   #pairs
		   $para{'pripairs'}=pairs($para{'pristr_2'},'loop');
		   $para{'prepairs'}=pairs($para{'prestr_2'},'loop');
		   $para{'mipairs'}=pairs($para{'mistr_2'});
		   #gu
		   $para{'prigu'}=gu($para{'pristr_2'},$para{'pristr_3'},'loop');
		   $para{'pregu'}=gu($para{'prestr_2'},$para{'prestr_3'},'loop');
		   $para{'migu'}=gu($para{'mistr_2'},$para{'mistr_3'});
		   #ntcontent
		   $para{'printcontent_a'}=ntcontent($para{'priseq'},'a');
		   $para{'printcontent_u'}=ntcontent($para{'priseq'},'u');
		   $para{'printcontent_c'}=ntcontent($para{'priseq'},'c');
		   $para{'printcontent_g'}=ntcontent($para{'priseq'},'g');
		   $para{'prentcontent_a'}=ntcontent($para{'preseq'},'a');
		   $para{'prentcontent_u'}=ntcontent($para{'preseq'},'u');
		   $para{'prentcontent_c'}=ntcontent($para{'preseq'},'c');
		   $para{'prentcontent_g'}=ntcontent($para{'preseq'},'g');
		   $para{'mintcontent_a'}=ntcontent($para{'miseq'},'a');
		   $para{'mintcontent_u'}=ntcontent($para{'miseq'},'u');
		   $para{'mintcontent_c'}=ntcontent($para{'miseq'},'c');
		   $para{'mintcontent_g'}=ntcontent($para{'miseq'},'g');
		   #firstbase
		   $para{'firstbase'}=uc(substr($para{'miseq'},0,1));
		   #pribulge#biggestbulge#bulge, number, type,size,position
		   ($para{'priinternalloopnumber'},$para{'priinternalloop'})=bulge(substr($para{'pristr_1'},$para{'length_basalsegment'}),-1);
		   ($para{'preinternalloopnumber'},$para{'preinternalloop'})=bulge($para{'prestr_1'},-1);
		   ($para{'miinternalloopnumber'},$para{'miinternalloop'})=bulge($para{'mistr_1'},0);
		   ($para{'internalloopnumber_lowerstem'},$para{'internalloop_lowerstem'})=bulge($para{'lowerstr_1'},0);
		   ($para{'internalloopnumber_topstem'},$para{'internalloop_topstem'})=bulge($para{'topstr_1'},0);
		   #unpairedbases,does not include the terminal loop, 1 refer to terminal loop, 0 refer no ternimal loop
		   $para{'priunpairedbases'}=unpairedbases(substr($para{'pristr_1'},$para{'length_basalsegment'}),substr($para{'pristr_4'},$para{'length_basalsegment'}),1);
		   $para{'preunpairedbases'}=unpairedbases($para{'prestr_1'},$para{'prestr_4'},1);
		   $para{'miunpairedbases'}=unpairedbases($para{'mistr_1'},$para{'mistr_4'},0);
		   $para{'unpairedbases_lowerstem'}=unpairedbases($para{'lowerstr_1'},$para{'lowerstr_4'},0);
		   $para{'unpairedbases_topstem'}=unpairedbases($para{'topstr_1'},$para{'topstr_4'},0);
		   #unpairedrate, not include terminal loop, the method good for budding stems--does not count budding stem
		   if (($para{'prilength'}-$para{'length_terminalloop'}) ne 0){
		   		$para{'priunpairedrate'}=sprintf '%4.4f', $para{'priunpairedbases'}/($para{'prilength'}-$para{'length_terminalloop'});
		   }
		   else{
		   		$para{'priunpairedrate'}=-1;
		   }
		   if (($para{'prelength'}-$para{'length_terminalloop'}) ne 0){
		   		$para{'preunpairedrate'}=sprintf '%4.4f', $para{'preunpairedbases'}/($para{'prelength'}-$para{'length_terminalloop'});
		   }
		   else{
		   		$para{'preunpairedrate'}=-1;
		   }
		   $tmp=$para{'mistr_1'}.$para{'mistr_2'}.$para{'mistr_3'}.$para{'mistr_4'};
		   $tmp=~s/[\s\-\.\\]+//g;
		   $para{'miunpairedrate'}=sprintf '%4.4f', $para{'miunpairedbases'}/length($tmp);
		   $tmp=$para{'lowerstr_1'}.$para{'lowerstr_2'}.$para{'lowerstr_3'}.$para{'lowerstr_4'};
		   $tmp=~s/[\s\-\.\\]+//g;
		   if (length($tmp) ne 0){
			  $para{'unpairedrate_lowerstem'}=sprintf '%4.4f', $para{'unpairedbases_lowerstem'}/length($tmp);
		   }
		   else{
			  $para{'unpairedrate_lowerstem'}=-1;
		   }
		   $tmp=$para{'topstr_1'}.$para{'topstr_2'}.$para{'topstr_3'}.$para{'topstr_4'};
		   $tmp=~s/[\s\-\.\\]+//g;
		   if (length($tmp) ne 0){
			  $para{'unpairedrate_topstem'}=sprintf '%4.4f', $para{'unpairedbases_topstem'}/length($tmp);
		   }
		   else{
			  $para{'unpairedrate_topstem'}=-1;
		   }
		   #stability
		   my $seq1;
		   my $seq2;
		   if ($para{'strand'} eq '5'){
			  $seq1=$para{'mistr_1'};
			  $seq2=$para{'mistr_2'};
		   }
		   else{
			  $seq1=$para{'mistr_3'};
			  $seq2=$para{'mistr_4'};
		   }
		   my $len1=4;#left
		   my $len2=4;#right
		   my $m="";
		   my $n="";
		   until((length($m)+length($n)) eq 4){#get the left four nts
			  $m=substr($seq1,0,$len1);
			  $n=substr($seq2,0,$len1);
			  $m=~s/[\s\-\.\\]+//g;
			  $n=~s/[\s\-\.\\]+//g;
			  $len1+=1;
		   }
		   $len1-=1;
		   $m="";
		   $n="";
		   until((length($m)+length($n)) eq 4){#get the right four nts
			  $m=substr($seq1,-$len2,$len2);
			  $n=substr($seq2,-$len2,$len2);
			  $m=~s/[\s\-\.\\]+//g;
			  $n=~s/[\s\-\.\\]+//g;
			  $len2+=1;
		   }
		   $len2-=1;
		   if ($para{'strand'} eq '5'){
			  $seq1=lc(join("",substr($para{'mistr_2'},0,$len1),substr($para{'mistr_3'},0,$len1)));
			  $seq2=lc(join("",substr($para{'mistr_2'},-$len2,$len2),substr($para{'mistr_3'},-$len2,$len2)));
		   }
		   else{
			  $seq1=lc(join("",substr($para{'mistr_2'},-$len2,$len2),substr($para{'mistr_3'},-$len2,$len2)));
			  $seq2=lc(join("",substr($para{'mistr_2'},0,$len1),substr($para{'mistr_3'},0,$len1)));
		   }
		   $seq1=~s/[\s\-\.\\]+//g;
		   $seq1=~ s/c/zz/g;
		   $seq2=~s/[\s\-\.\\]+//g;
		   $seq2=~ s/c/zz/g;
		   if (length($seq2) ne 0){#incase $seq2 eq nothing
				$para{'stability'}=sprintf '%.2f',length($seq1)/length($seq2);
		   }
		   else{
				if (length($seq1) eq 0){
					$para{'stability'}=1.00;
				}
				else{
					$para{'stability'}=length($seq1)+1;#stabily, if seq2 eq nothing, stability eq to $seq1's length plus one.
				}
		   }
		   #overhang
		   $tmp="";
		   $m=1;
		   while ((length($tmp) < 2)  and ($para{'upperend'}+$m)<=length($para{'pristr_3'})){
			  if ($para{'strand'} eq "5"){#generate the two nt
				 $tmp=join("",substr($para{'pristr_3'},$para{'upperend'},$m),
						 substr($para{'pristr_4'},$para{'upperend'},$m));
			  }
			  else{
				 $tmp=join("",substr($para{'pristr_1'},$para{'upperstart'}-$m-1,$m),
						 substr($para{'pristr_2'},$para{'upperstart'}-$m-1,$m));
			  }
			  $tmp=~ s/[\s\\-]+//g;
			  $m++;
		   }
		   if ($tmp eq ""){#no overhang
			  $para{'penultimateposition'}="NULL";
			  $para{'terminalnucleotide'}="NULL";
		   }
		   if (length($tmp) eq 1){#one overhang
			  $para{'penultimateposition'}=$tmp;
			  $para{'terminalnucleotide'}="NULL";
		   }
		   if (length($tmp) eq 2){#two overhang
			  if ($para{'strand'} eq "5"){
				 $para{'penultimateposition'}=substr($tmp,0,1);
				 $para{'terminalnucleotide'}=substr($tmp,1,1);
			  }
			  else{
				 $para{'penultimateposition'}=substr($tmp,1,1);
				 $para{'terminalnucleotide'}=substr($tmp,0,1);
			  }
		   }
		   $para{'penultimateposition'}=uc($para{'penultimateposition'});
		   $para{'terminalnucleotide'}=uc($para{'terminalnucleotide'});
		   #################################################################
		   #pripredict
		   my @pmt=pmt("range");
		   my $prirst;
		   foreach (@pmt){
			  $tmp=$_;
			  $prirst=pricompare($para{$tmp},$lower{$tmp},$upper{$tmp});
			  if ($prirst eq 0){
				 goto pripredictfail;
			  }
		   }
		   #################################################################
		   #output the parameters
		   @pmt=pmt("all");
		   my $printdata="";#create print data to avoid different cores write into the data at the same time
		   foreach (@pmt){
			  $printdata.="$para{$_}\t";
		   }
		   $printdata=~s/\t$//;
		   $printdata.="\n";
		   open (OUT,">>$syspmt{'outfilepmt'}");
		   print OUT $printdata;
		   close OUT;
		   #################################################################
		   #miRNA prediction
		   if ($syspmt{'onlypmt'} ne 1){#calculate the parameter only without prediction
		   		$printdata=~s/ /\_/g;
		   		$printdata=~s/\t/,/g;
		   		system ("perl $syspmt{'share'}/miRPara_predict.pl $syspmt $printdata");
		   }
		   pripredictfail:
		}#end if ($para{'strand'} ne 0)
	 }#end of if ($para{'miend'}<=$para{'prilength'}){#incase out of range
  }#end of foreach (@milen)
}#end of for ($para{'mistart'}=1;$para{'mistart'}<=length($para{'priseq'});$para{'mistart'}+=1)
close PMT;
close OUT;
exit;
sub gc{#parameter
	#===========The gc was used to generate gc content of given seq
	my $seq;#receive the sequences
	my $tmp;#tmp string
	my $gc;
	($seq)=@_;
	$tmp=lc($seq);
	$tmp=~s/[au]+//g;
	$gc=sprintf '%4.4f', length($tmp)/length($seq);#prigc
	return $gc;
}
sub pairs{#parameter
	#===========The pairs was used to generate the number of pairs of second str
	my $seq;
	my $loop;
	my $pairs;
	($seq,$loop)=@_;
	if (defined $loop){
		if (lc($loop) eq "loop"){
			chop($seq);
		}
	}
	$seq=~ s/\s+//g;
	$pairs=length($seq);
	return $pairs;
}
sub gu{#parameter
	#===========The gu was used to generate the number of GU wobbles
	my $seq1;
	my $seq2;
	my $tmp;#tmp string
	my $loop;
	my $gu;
	my $u;
	my $a;
	($seq1,$seq2,$loop)=@_;
	if (defined $loop){
		if (lc($loop) eq "loop"){
			chop($seq1);
			chop($seq2);
		}
	}
	$tmp=join("",$seq1,$seq2);
	$tmp=~ s/\s+//g;
	$u=lc($tmp);
	$u=~ s/[acg]+//g;
	$a=lc($tmp);
	$a=~ s/[ucg]+//g;
	$gu=length($u)-length($a);
	return $gu;
}
sub ntcontent{#parameter
   #===========The ntcontent was used to generate the nt content of four nts
   my $seq;
   my $tmp;
   my $nt;
   my $ntcontent;
   ($seq,$nt)=@_;
   $tmp=lc($seq);
   $nt=lc($nt);
   $tmp=~ s/[$nt]+//g;
   $ntcontent=sprintf '%4.4f', (length($seq)-length($tmp))/length($seq);
   return $ntcontent;
}
sub bulge{#parameter
   #the bulge was used to generate the number, the biggest size
   #	receive the first
   #	report the number of bulges and the biggest internal loop
   my ($seq,$del)=@_;
   my $number;
   my $size=0;
   my @tmp;
   my $tmp;;
   $seq=~s/\s+$//;
   $seq=~s/^\s+//;
   $seq=~s/\s+/,/g;
   @tmp=split(",",$seq);
   if ($del eq -1){
	  pop @tmp;
   }
   #number
   $number=@tmp;
   if ($number ne 0){#in case no any bulge
	  #size
	  foreach (@tmp){
		 if (length($_)>$size){
			$size=length($_);
		 }
	  }
   }
   else{
	   $number=0;
	   $size=-1;
   }
   return $number,$size;
}
sub unpairedbases{#parameter
	#===========The unpairedbases was used to generate the size of unpaired bases
	#does not include the terminal loop
	my ($seq1,$seq2,$loop)=@_;
	my $tmp;
	my @seq;
	my @tmp;
	if ($loop eq 1){
	  @seq=($seq1,$seq2);
	  foreach (@seq){
		 $_=~s/\s+$//g;
		 $_=~s/^\s+//g;
		 $_=~s/\s+/,/g;
		 @tmp=split(",",$_);
		 pop (@tmp);
		 $_=join("",@tmp);
	  }
	  $seq1=$seq[0];
	  $seq2=$seq[1];
	}
	
	$tmp=$seq1.$seq2;
	$tmp=~s/[\s\-]+//g;
	return length($tmp)
}
sub nt2number{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The nt2number was used to translate the nt to ACC number
   #===========Usage: nt2number(nt)
   #===========The number of acc will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my $value;
   ($value)=@_;
   $value=lc($value);
   if ($value eq "null"){
	  $value="-1";
   }
   elsif ($value eq ""){
	  $value="";#if nothing
   }
   elsif ($value eq "a"){
	  $value="1";
   }
   elsif ($value eq "c"){
	  $value="2";
   }
   elsif ($value eq "g"){
	  $value="3";
   }
   elsif (($value eq "u") or ($value eq "t")){
	  $value="4";
   }
   elsif ($value eq "n"){
	  $value="0";
   }
   else{
	  $value="-1";
   }
   return $value;
}
sub pricompare{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The compare was used to decide whether a value in in certain range
   #===========Usage: compare(n1,n2,n3,)
   #===========the value 1 for true or 0 for flase will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my $n1="";#given number
   my $n2="";#lower
   my $n3="";#upper
   my $result;
   ($n1,$n2,$n3)=@_;
   if (lc($n1) eq "null" or lc($n1) eq "a" or lc($n1) eq "u" or lc($n1) eq "c" or lc($n1) eq "g" or lc($n1) eq ""){
	  $n1=nt2number($n1);
   }
   if ($n1 eq "-1"){
	  $result=1;
   }
   else{
	  if ($n1>=$n2 and $n1<=$n3){
		 $result=1;
	  }
	  else{
		 $result=0;
	  }
   }
   return $result;
}
sub prilimit{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The range was used to provide range information
   #===========Usage: range(parameter)
   #			--up all parameters
   #			--down the range parameters
   #===========the corespond value will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my ($id)=@_;
   my @range;
   my %range;
   my @rangeparameter=pmt("range");#the range parameter
   my $m=0;
   if (lc($id) eq "lower"){
	  @range=("17","41","49","0","-1","17","-1","3","10","12","15","-172.2",
			  "-307.7","0.16","0.21","0","0","0","0","0.05","0.06","0.08",
			  "0.07","0.05","0.07","0.1","0.08","0","0","-1","-1","0","0",
			  "-1","-1","0","0","0","-1","-1","0","0","0","-1","-1","0",
			  "1","20","1","19","-1","-1","0.2","0","0","0");
   }
   elsif (lc($id) eq "upper"){
	  @range=("27","250.96","379","24","116.42","35","103.85","27.42","23",
			  "100.42","150","-6.16","-20","0.86","0.83","0.57","0.59",
			  "0.7","0.6","0.43","0.44","0.49","0.47","0.42","0.43","0.5",
			  "0.45","12","17.42","16","17.42","6","20.85","13.42","12","26",
			  "57.42","80.92","56.54","52.85","0.51","0.47","0.46","0.64",
			  "0.75","6","258.38","279.65","117.42","136.42","5.21","14",
			  "0.85","23.85","9","12");
   }
   for ($m=0;$m<@rangeparameter;$m+=1){
	  $range{$rangeparameter[$m]}=$range[$m];
   }
   return %range;
}
sub pmt{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The pmt was used to provide pmt files
   #===========Usage: pmt(parameter)
   #			--all all parameters
   #			--range the range parameters
   #			--svmoverall the svm parameters of all species
   #			--svmoanimal the svm parameters of animal
   #			--svmplant the svm parameters of plants
   #			--svmvirus the svm parameters of all virus
   #			--display the display parameters
   #===========the parameters will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my ($class)=@_;
   my @parameter;
   if ($class eq "Viruses"){
   		$class="Metazoa";
   }
   pmtstart:
   if (lc($class) eq "range"){
	  @parameter=("milength","prelength","prilength","length_basalsegment",
				  "length_lowerstem","length_upperstem","length_topstem",
				  "length_terminalloop","mipairs","prepairs","pripairs",
				  "premfe","primfe","migc","prigc","mintcontent_a","mintcontent_c",
				  "mintcontent_g","mintcontent_u","prentcontent_a","prentcontent_c",
				  "prentcontent_g","prentcontent_u","printcontent_a","printcontent_c",
				  "printcontent_g","printcontent_u","miinternalloop","preinternalloop",
				  "internalloop_lowerstem","internalloop_topstem","miinternalloopnumber",
				  "priinternalloopnumber","internalloopnumber_lowerstem",
				  "internalloopnumber_topstem","miunpairedbases","preunpairedbases",
				  "priunpairedbases","unpairedbases_lowerstem","unpairedbases_topstem",
				  "miunpairedrate","preunpairedrate","priunpairedrate","unpairedrate_lowerstem",
				  "unpairedrate_topstem","migu","mistart","miend","upperstart",
				  "upperend","stability","preinternalloopnumber","pregc",
				  "priinternalloop","pregu","prigu");#56
   }
   elsif (lc($class) eq "overall"){
        @parameter=("unpairedrate_lowerstem","prelength","internalloopnumber_lowerstem","length_upperstem",
        "miinternalloop","firstbase","mintcontent_a","migc","pregc","prentcontent_u",
        "prentcontent_a","internalloop_topstem","preunpairedrate","mipairs","prepairs",
        "internalloopnumber_topstem","unpairedrate_topstem","mintcontent_c","mistart",
        "miunpairedrate","mintcontent_g","terminalnucleotide","prentcontent_c",
        "prentcontent_g","mintcontent_u");#25pmt
   }
   elsif (lc($class) eq "metazoa"){
        @parameter=("internalloop_topstem","internalloopnumber_topstem","length_topstem","length_upperstem",
        "migc","miinternalloop","miinternalloopnumber","mintcontent_a","mintcontent_c","mintcontent_g",
        "mintcontent_u","mistart","miunpairedrate","penultimateposition","pregc","prelength",
        "prentcontent_a","prentcontent_c","prentcontent_g","prentcontent_u","preunpairedrate",
        "stability","unpairedrate_lowerstem","unpairedrate_topstem");#24pmt
   }
   elsif (lc($class) eq "viridiplantae"){
        @parameter=("firstbase","internalloop_topstem","internalloop_lowerstem","internalloopnumber_lowerstem",
        "length_upperstem","migc","migu","miinternalloop","miinternalloopnumber","mintcontent_a",
        "mintcontent_g","mintcontent_u","mipairs","penultimateposition","pregc","prentcontent_a",
        "prentcontent_c","prentcontent_g","prentcontent_u","preunpairedrate","stability",
        "unpairedrate_lowerstem","unpairedrate_topstem","upperstart");#24pmt
   }
   elsif (lc($class) eq "all"){
	  @parameter=("miid","priid",
	  "miseq","mistr_1","mistr_2","mistr_3","mistr_4",
	  "preseq","prestr_1","prestr_2","prestr_3","prestr_4",
	  "priseq","pristr_1","pristr_2","pristr_3","pristr_4",
	  "milength","prelength","prilength",
	  "length_basalsegment","length_lowerstem","length_upperstem","length_topstem","length_terminalloop",
	  "mipairs","prepairs","pripairs",
	  "premfe","primfe",
	  "migc","pregc","prigc",
	  "mintcontent_a","mintcontent_c","mintcontent_g","mintcontent_u",
	  "prentcontent_a","prentcontent_c","prentcontent_g","prentcontent_u",
	  "printcontent_a","printcontent_c","printcontent_g","printcontent_u",
	  "miinternalloop","preinternalloop","priinternalloop",
	  "internalloop_lowerstem","internalloop_topstem",
	  "miinternalloopnumber","preinternalloopnumber","priinternalloopnumber",
	  "internalloopnumber_lowerstem","internalloopnumber_topstem",
	  "miunpairedbases","preunpairedbases","priunpairedbases",
	  "unpairedbases_lowerstem","unpairedbases_topstem",
	  "miunpairedrate","preunpairedrate","priunpairedrate",
	  "unpairedrate_lowerstem","unpairedrate_topstem",
	  "migu","pregu","prigu",
	  "strand","firstbase","penultimateposition","terminalnucleotide",
	  "mistart","miend","upperstart","upperend","stability");#77
   }
   else{
   		print "Warning: unrecognized species group $class, use Overall instead!\n";
   		$class="Overall";
   		goto pmtstart;
	#  @parameter=();
   }
   return @parameter;
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
