#! /usr/bin/perl
#	by Yonggan Wu (weasteam@gmail.com)
#	by Simon Rayner (raynere@wh.iov.cn)
#################################USE############################################
use warnings;
#use strict;
use Cwd qw(abs_path);
use Cwd;
use File::chdir;
#################################GET PARAMETERS#################################
my $usage="perl $0 <syspmt> <parameters>\n\n";
my $syspmt=shift or die $usage;
my $para=shift or die $usage;
$para=~s/\_/ /g;
my @pmt=pmt("all");
$para=~s/\n//;
my @para=split(",",$para);
my %para;
for (my $i =0;$i<@pmt;$i++){
	$para{$pmt[$i]}=$para[$i];
}
our %syspmt;
syspmt($syspmt);
my %mirbase=mirbase();
#update file directory
my $tmp=$syspmt;
$tmp=~s/\//,/g;
my @tmp=split(",",$tmp);
$dir=$syspmt;
$tmp=pop @tmp;
$dir=~s/\/$tmp//;
$CWD = $dir;
our %organisms=organisms();#species group name
#################################USE############################################
print "Predicting for $para{'miid'}...\n";
#	predict the given result
#	receive the parameters
#	print out the correct results
@pmt=pmt($organisms{lc($syspmt{'species'})});
my $key=int(rand(999999999));
open (OUT,">$key.txt");
my $printdata="2";
my $count=1;
foreach (@pmt){
  my $value;
  if (lc($para{$_}) eq "null" or $_ eq "firstbase" or $_ eq "penultimateposition" or $_ eq "terminalnucleotide"){
	 $value=nt2number($para{$_});
  }
  else{
	 $value=$para{$_};
  }
  $printdata.="\t$count:$value";
  $count++;
}
print OUT "$printdata\n";
close OUT;
#now predict
system "$syspmt{'svmpredict'} -q -b 1 $key.txt $syspmt{'svmmodel'} $key\_result.txt";
my $rst=0;
open (IN,"$key\_result.txt");
@tmp=<IN>;
close IN;
system "rm $key*";
@tmp=split(" ",$tmp[1]);
$rst=$tmp[1];
if ($rst>=$syspmt{'cutoff'}){
	if ($para{'miseq'}=~/[actug]/){
		$para{'miseq'}=uc($para{'miseq'});
	}
	if ($para{'miseq'}=~/T/){
		$para{'miseq'}=~s/T/U/g;
	}
	if (exists $mirbase{$para{'miseq'}}){
	 $mirbase=$mirbase{$para{'miseq'}};
	}
	else{
	 $mirbase="";
	}
	open (OUT,">>$syspmt{'outfilerst'}");
	print OUT "$para{'priid'}\t$para{'priseq'}\t$para{'miid'}\t$para{'miseq'}\t$para{'strand'}P\t$rst\t$mirbase\n";
	close OUT;
}
exit;
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
sub mirbase{
	#get a list of miRNA short name and sequences.
	my $mirbase="$syspmt{'share'}/models/miRBase/current";
	my $file="mature.fa";
	my %seq;
	if (not -e "$mirbase/$file"){
		#download the most recent data from mirbase
		download($mirbase,$file);
	}
	if (-e "$mirbase/$file"){
		my $title;
		open (IN,"$mirbase/$file");
		while (<IN>){
			$_=~s/[\n\r]//g;
			if ($_=~/>/){
				my @tmp=split(" ",$_);
				$title=pop @tmp;
			}
			else{
				if (exists $seq{$_}){
					if ($seq{$_}!~/$title/){
						$seq{$_}.=";$title";
					}
				}
				else{
					$seq{$_}=$title;
				}
			}
		}
		close IN;
	}
	return %seq;
}
sub download{
	#to check the download mirbase files
	my ($folder,$name)=@_;
	#download the most recent data from mirbase
	if (not -e $folder){
		system "mkdir -p $folder";
	}
	print "Warning: Not found $name!\nGoing to download from miRBase...\n";
	system "wget -P $syspmt{'share'}/models/miRBase/current/ ftp://mirbase.org/pub/mirbase/CURRENT/$name.gz";
	system "gunzip $syspmt{'share'}/models/miRBase/current/$name.gz";
}
sub organisms{
	#get a list of miRNA short name and sequences.
	my $mirbase="$syspmt{'share'}/models/miRBase/current";
	my $file="organisms.txt";
	my %spe;
	$spe{"overall"}="Overall";
	$spe{"metazoa"}="Metazoa";
	$spe{"viruses"}="Viruses";
	$spe{"viridiplantae"}="Viridiplantae";
	$spe{"animal"}="Metazoa";
	$spe{"virus"}="Viruses";
	$spe{"plant"}="Viridiplantae";
	if (not -e "$mirbase/$file"){
		#download the most recent data from mirbase
		download($mirbase,$file);
	}
	if (-e "$mirbase/$file"){
		my $title;
		open (IN,"$mirbase/$file");
		while (<IN>){
			$_=~s/[\n\r]//g;
			my @tmp=split("\t",$_);
			my @dat=split(";",$tmp[3]);
			$spe{$tmp[0]}=$dat[0];
		}
		close IN;
	}
	return %spe;
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
