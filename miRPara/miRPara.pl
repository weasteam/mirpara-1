#! /usr/bin/perl
#	by Yonggan Wu (weasteam@gmail.com)
#	by Simon Rayner (raynere@wh.iov.cn)
#################################USE############################################
use warnings;
#use strict;
use Cwd qw(abs_path);
use Cwd;
use File::chdir;
use Getopt::Long;
##################################PARAMETER#####################################
#parameters
my $dir = getcwd;
our %syspmt=();
$syspmt{'abname'}="unknown";#update to guess from file name
$syspmt{'species'}="overall";
$syspmt{'level'}=1;#default level 1
$syspmt{'prilengthlimit'}=60;
$syspmt{'version'}="7.0";
$syspmt{'foldlen'}=521;
$syspmt{'overlaplen'}=176;
$syspmt{'cutoff'}=0.8;#cutoff probabilities
$syspmt{'onlypmt'}=0;#only calculate the pmt valude instead of do further prediction
$syspmt{'svmmodel'}="";
$syspmt{'share'}=abs_path($0);#absolute path, can not replace with $0;
my %options;
promgramCheck();#program/environment, check existance of UNAFold.pl and ct2out
################################PROGROM START###################################
GetOptions \%options, 'version|v' => sub { version('miRPara.pl') },
						'help|h' => sub { usage() },
						'species|s=s' => \$syspmt{'species'},
						'level|l=i'=>\$syspmt{'level'},
						'cutoff|c=f'=>\$syspmt{'cutoff'},
						'foldlength|f=i'=>\$syspmt{'foldlen'},
						'overlaplength|o=i'=>\$syspmt{'overlaplen'},
						'pmt'=>\$syspmt{'onlypmt'},
						'prilengthlimit|p=i' => \$syspmt{'prilengthlimit'} or die("Type perl $0 -h for help!\n");
#################################data prepare###################################
##infile
my $infile=shift or die("Error: data not specified\nRun 'perl $0 -h' for help\n");
#define the share directory and working directory
##share directory - mirpara.pl directory, for models and mirbase
##working directory - rnafold output directory, svm predict directory
$syspmt{'share'}=~s/\/$0//;
$syspmt{'share'}=~s/\n//g;
$syspmt{'share'}=~s/miRPara.pl//;
$syspmt{'share'}=~s/\/$//;
if ($infile=~/^\./){#change to full directory
	$dir=~s/\/$//;
	$infile=~s/\.\///;
	$infile="$dir/$infile";
}
if ($syspmt{'species'}=~/^\./){#change the species name, if treat as model path, before path changes
	$syspmt{'species'}=~s/\.\///;
	$syspmt{'species'}="$dir/$syspmt{'species'}";
}
if ($infile=~/\//){
	my $tmp=$infile;
	$tmp=~s/\//,/g;
	my @tmp=split(",",$tmp);
	$dir=$infile;
	$infile=pop @tmp;
	$dir=~s/\/$infile//;
}
$CWD = $dir;
##outfiles
if ($infile=~/\.pmt/){
   $syspmt{'outfilerst'}="$infile.out";
   $syspmt{'outfilerst'}=~s/\.pmt//;
}
else{
   my $tmp=$infile;
   $tmp=~s/\.fasta//;
   $tmp=~s/\.fas//;
   $tmp=~s/\.fa//;
   $syspmt{'outfilerst'}="$tmp.out";
   $syspmt{'outfilepmt'}="$tmp.pmt";
}
##species
if ($syspmt{'species'} =~/\//){#the model address provided, share, species, level should be changed
	$syspmt{'svmmodel'}=$syspmt{'species'};
	$syspmt{'species'}=~s/\//,/g;
	my @tmp=split(",",$syspmt{'species'});
	$syspmt{'species'}=pop @tmp;
	@tmp=split("_",$syspmt{'species'});
	$syspmt{'species'}=$tmp[0];
	$tmp[1]=~s/\.model//;
	$syspmt{'level'}=$tmp[1];
}
elsif ($syspmt{'species'}=~/animal/ or $syspmt{'species'}=~/plant/ or $syspmt{'species'}=~/virus/){#to check the species
	$syspmt{'species'}=lc($syspmt{'species'});
	my %spenames;
	$spenames{"animal"}="Metazoa";
	$spenames{"plant"}="Viridiplantae";
	$spenames{"virus"}="Viruses";
	$syspmt{'species'}=$spenames{$syspmt{'species'}};
}
else{
	$syspmt{'species'}=lc($syspmt{'species'});
	if ($syspmt{'species'} ne "overall"){
		#change the abname into species short name
		#abname will be assigned same to species name if the species name is provided
		$syspmt{'abname'}=$syspmt{'species'};
	}
}
#check cutoff
if ($syspmt{'cutoff'}>1){
	print "Warning: the cutoff ($syspmt{'cutoff'}) is too high!!!\n";
}
#foldlength and overlap length check
if ($syspmt{'foldlen'}<=$syspmt{'overlaplen'}){
   die "Error: The cut length should be longer than overlap length\n!";
}
if ($syspmt{'foldlen'}<=0 or $syspmt{'overlaplen'}<=0){
   die "Error: The cut length or overlap length is not valid!\n";
}
#load the models
if ($syspmt{'onlypmt'} ne 1){
	if ($syspmt{'svmmodel'} eq ""){#if the model is not provided from -s key
		$syspmt{'svmmodel'} = "$syspmt{'share'}/models/$syspmt{'species'}\_$syspmt{'level'}.model";
	}
	if (not -e $syspmt{'svmmodel'}){
		die "Error: Could not find the model($syspmt{'svmmodel'}), please type miRPara.pl -h for help!\n";
	}
}
$syspmt{'outfilerst'}=~s/\.out/_level_$syspmt{'level'}\.out/;
#save a copy of runing parameters
open (OUT,">mirpara_parameters.txt");
print OUT "#PARAMETER\tVALUE\n";
my @tmp=keys (%syspmt);
foreach (@tmp){
	print OUT "$_\t$syspmt{$_}\n";
}
close OUT;
##################################start#########################################
#1	generate the sequences;
#2	do split
#output the predict results
open (OUT,">$syspmt{'outfilerst'}");
print OUT qq{#miRNA results predicted by miRPara$syspmt{'version'}
#By Yonggan Wu and Simon Rayner
#Report bugs to Yonggan Wu (weasteam\@gmail.com) or Dr. Simon Rayner (sraynere\@wh.iov.cn)
#Bioinformatics Group, State Key Laboratory of Virology, Wuhan Institute of Virology, Chinese Academy of Sciences.
#The Department of Biological Sciences, Texas Tech University.
#Homepage: http://www.whiov.ac.cn/bioinformatics/mirpara
#Google Project: http://code.google.com/p/mirpara/
#Facebook: https://www.facebook.com/mirpara2009
#model: $syspmt{'svmmodel'}
#cutoff: $syspmt{'cutoff'}
#foldlength: $syspmt{'foldlen'}
#overlaplength: $syspmt{'overlaplen'}
#
#priid\tpriseq\tmiid\tmiseq\tstrand\tSVM_probability\tmiRBase
};
close OUT;
################################################################################
#output the pmt results
open (OUT,">$syspmt{'outfilepmt'}");
#print the titile of parameter
my @pmt=pmt("all");
my $tmp=shift @pmt;
print OUT "#$tmp";
foreach (@pmt){
  print OUT "\t$_";
}
print OUT "\n";
close OUT;
################################################################################
#pmt file prediction
if ($infile =~ /\.pmt/){
   my %para;
   my $dat;
   my @filepmt;
   open (IN,$infile);
   while (<IN>){
	  $_=~s/\n//;
	  if ($_=~/^#/){
			$_=~s/#//;
			@filepmt=split("\t",$_)
	  }
	  else{
			$_=~s/\n//;
			%para=();
			my @tmp=split("\t",$_);
			my @pmt=@filepmt;
			my $tmp=@pmt;
			for (my $dat=1;$dat<=$tmp;$dat+=1){
			   $para{$pmt[$dat-1]}=$tmp[$dat-1];
			}
			##########################################################################
			#parameter prediction
			#$pm->start and next; # do the fork
			print "Predicting for $para{'miid'}...";
			predict(%para);
			#$pm->finish;
	  }
   }
   close IN;
   #$pm->wait_all_children;
}
elsif ($infile =~/\.fa/){
	#split the long input files into short sequences
	open (IN,$infile) or die ("Error: $infile does not exist!\n");
	my $title;
	my $seq="";
	my %seq;
	while (<IN>){
		$_=~s/[\n\r]//g;
		if ($_=~/^>/){
			if ($seq ne ""){
				$seq{$title}=$seq;
			}
			$_=~s/>//;
			$_=~s/[\s\|]/,/g;
			my @tmp=split(",",$_);
			$title=$tmp[0];
			$seq="";
		}
		else{
			$_=~s/\s+//g;
			$_=lc($_);
			$_=~s/t/u/g;
			$seq.=$_;
		}
	}
	close IN;
	if ($seq ne ""){#the last sequences
		$seq{$title}=$seq;
	}
	my @seq=keys(%seq);
	foreach (@seq) {
		unless ($pid = fork) {
			unless (fork) {
				system "perl $syspmt{'share'}/miRPara_split_seq.pl $dir/mirpara_parameters.txt $_ $seq{$_}";
			}
			exit 0;
		}
		waitpid($pid, 0);
	}
}#whether only calculate the parameter
else{
   die "Error: data format unrecongnized,\nRun 'miRPara.pl -h' for help\n";
}
#calculate the parameters only
if ($syspmt{'onlypmt'} eq 1){
   system "rm $syspmt{'outfilerst'}";
}
print "\nAll Done.\nPlease wait until all backgroup perl programs were done!\n";
##################web serve#####################################################
#system "tar zcf miRPara_result.tar.gz ./*";#web serve
##################STOP HERE#####################################################
sub usage {
    print qq{
Usage: miRPara.pl [Options] file [File]

Options:
-v, --version
-h, --help
-s, --species=species group name <overall, animal, plant, virus> (defaults as overall)
              or species short name <hsa, mmu, cel, ...>
              or path of a model file: /home/xxx/miRPara/model/hsa_10.model
-c, --cutoff=<the cutoff to svm prediction probabilities> (defaults to 0.8)
-l, --levels=<1..20>(defaults to 1)
-p, --prilengthlimit=<limit to the pri-miRNA length> (defaults to 60)
-f, --foldlength=<The length to be split for folding> (defaults to 500)
-o, --overlaplength=<The overlap length of two nearby splited fragments> (defaults to 150)
--pmt, --Only calculate the parameters without further prediction

File (one of following):
*.fasta, --Only fasta format sequences were allowed
*.pmt, --Repredict from the parameter files

Report bugs to Yonggan Wu (weasteam\@gmail.com) or Dr. Simon Rayner (raynere\@wh.iov.cn)

Homepage: http://www.whiov.ac.cn/bioinformatics/mirpara
Google Project: http://code.google.com/p/mirpara/
Facebook: https://www.facebook.com/mirpara2009

};
    exit;
}
sub version ($) {
	print qq {
$_[0] $version

By Yonggan Wu (weasteam\@gmail.com) and Simon Rayner (raynere\@wh.iov.cn)
Copyright (C) 2009
Wuhan Institute of Virology, Chinese Academic of Science. Wuhan, 430071, China
Department of Biological Sciences, Texas Tech University, Texas, 33409, USA

Homepage: http://www.whiov.ac.cn/bioinformatics/mirpara
Google Project: http://code.google.com/p/mirpara/
Facebook: https://www.facebook.com/mirpara2009

};
    exit;
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
sub promgramCheck{
	my %pro;
	$pro{"unafold"}="UNAFold.pl";
	$pro{"ct2out"}="ct2out";
	$pro{"svmpredict"}="svm-predict";
	my @tmp=keys(%pro);
	foreach (@tmp){
		$syspmt{$_}=`which $pro{$_}`;
		$syspmt{$_}=~s/\n//;
		if ($syspmt{$_} eq ""){
			die "Error: $pro{$_} not found!\n";
		}
	}
}
