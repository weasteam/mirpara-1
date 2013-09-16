#! /usr/bin/perl
#	by Yonggan Wu (weasteam@gmail.com)
#	by Simon Rayner ()
#	Version: 2.0 (2012-10-02 15:44:02 )
#	version 2.1 (2012-10-03 17:52:47 ): fix the bug to download miRBase17 data
#	version 2.2 (2012-10-09 12:51:39 ): ignore empty data set, fix a bug for species group
#	version 2.3 (2012-10-12 18:30:07 ): fix a bug that not able to produce parameters for some species
#	version 3.0 (2012-11-17 09:40:06 ): use libsvm to train
#										: able to run the whole process separatly for different step
#										: able to run for one specific species
#										: multiple core support
#	version 3.1 (2012-11-23 20:21:12 ): enable multicore processing for the parameter calculation
#	version 3.2 (2012-11-24 15:32:17 ): fix a bug the terminal the program after step 3
#	The script will prepare the data for SVM training from miRBase release
#	Changes since last version: calculate the parameters first and then group
################################################################################
#use strict;
use warnings;
use File::chdir;
use Cwd;
use threads;
#use threads::shared;
#use Algorithm::SVM;
#use Algorithm::SVM::DataSet;
my $dir = getcwd;
$CWD = $dir;
my $usage=qq{
Usage: perl $0 <miRBase Version> <Full miRPara address> <steps> <cores> <species>
Example: perl $0 19 /usr/bin/miRPara.pl 1,2,3,4 12 all

Steps:
1	extract experimentally varified miRNA sequences
2	calculate the miRPara parameters
3	generate the parameters for training
4	model training

Species (one of following):
<all> - all species
<species short name> - eg. hsa

};
our $version=shift or die $usage;
if ($version<13){
	die "Error: sorry, we do not support miRBase before verison 13.0\n";
}
my $mirpara=shift or die $usage;
my $step=shift or die $usage;
$step=",".$step.",";
my $cores=shift or die $usage;
my $spename=shift or die $usage;#species names
my $folder="./miRBase/$version";
my %fh;
our %spe;
our %speg;
if ($step !~/,1,/){
	goto step2;
}
################################################################################
#step 1 extract experimental varified miRNA sequences
################################################################################
my @files=("*THIS_IS_RELEASE_*","organism*","miRNA.dat*","hairpin.fa*","mature.fa*");
#check the given folder parameter
if (-e $folder){# if found the data
	$folder=~s/\/$//;
	my @tmp=glob("$folder/$files[0]");
	print "Checking required files...\n";
	for (my $i=1;$i<@files;$i++){
		my @tmp=glob("$folder/$files[$i]");
		if (@tmp eq 0){
			print "Warning: Not found $folder/$files[$i]!\nGoing to download from miRBase...\n";
			system "wget -P ./miRBase/$version/ ftp://mirbase.org/pub/mirbase/$version/$files[$i]";
		}
	}
	print "All Found\n";
}
else{#if not found
	system "mkdir -p ./miRBase/$version/";
	foreach (@files){
		system "wget -P ./miRBase/$version/ ftp://mirbase.org/pub/mirbase/$version/$_";
	}
}
#prepare the evidence infor
#use pre-miRNA accession number
my %evi;
if (not -e "./miRBase/$version/miRNA.dat"){
	system "gunzip ./miRBase/$version/miRNA.dat.gz";
}
print "Extracting the evidence from miRNA.dat...";
open (IN,"./miRBase/$version/miRNA.dat");
my $ac;
my $evi="";
my %evi_group;
while (<IN>){
	if ($_=~/^AC/){
		$_=~s/;\n//;
		$_=~s/AC\s+//;
		$ac=$_;
	}
	if ($_=~/^FT                   \/evidence/){
		$_=~s/FT                   \/evidence//;
		$_=~s/\n//;
		$_=~s/=//;
		if ($evi eq "" or $evi ne "experimental"){#if two miRNAs, will select if one one of them is from experimental
			$evi=$_;
		}
	}
	if ($_=~/^SQ/){
		$evi{$ac}=$evi;
		$evi_group{$evi}=1;
		$evi="";
	}
}
close IN;
print "Done\n";
#read species ifnor
readspe();
#finally, output the pre-miRNA sequences
print "Outputing the experimentally varified sequences...";
if (not -e "experimental"){
	system "mkdir experimental";
}
if (not -e "./miRBase/$version/hairpin.fa"){
	system "gunzip ./miRBase/$version/hairpin.fa.gz";
}
my @evi=keys(%evi_group);
foreach (@evi){
	$evi=$_;
	open ($fh{$evi},">./experimental/miRBase_$version\_$evi\_overall.fasta");
}
open (IN,"./miRBase/$version/hairpin.fa");
my $spe;
while (<IN>){
	if ($_=~/^>/){
		my @tmp=split(" ",$_);
		if (exists $evi{$tmp[1]}){
			$evi=$evi{$tmp[1]};
		}
		else{
			$evi="11";
		}
	}
	if (exists $fh{$evi}){
		print {$fh{$evi}} $_;
	}
}
close IN;
foreach (@evi){
	$evi=$_;
	close $fh{$evi};
}
print "Done\n";
################################################################################
#step2	calculate the miRPara parameters
################################################################################
step2:
if ($step !~/,2,/){
	goto step3;
}
print "Calculating the parameters...";
system "perl $mirpara -t $cores --pmt ./experimental/miRBase_$version\_experimental_overall.fasta";
print "Done\n";
################################################################################
#step3	generate the parameters for training
################################################################################
step3:
if ($step !~/,3,/){
	goto step4;
}
print "Generating the training parameters...\n";
my $out="training_parameters_$version";
if (not -e $out){
	system "mkdir $out";
}
#unzip the data
if (not -e "./miRBase/$version/mature.fa"){
	system "gunzip ./miRBase/$version/mature.fa.gz";
}
#get the mature sequences
print "Reading the mature sequences...\n";
my %mature;
my $title;
open (IN,"./miRBase/$version/mature.fa");
while (<IN>){
	$_=~s/\n//;
	if ($_=~/^>/){
		my @tmp=split(" ",$_);
		$title=$tmp[0];
		$title=~s/>//;
		$title=lc($title);
	}
	else{
		$mature{$_}.="$title,";
	}
}
close IN;
#output the sequences
print "Done\nExtracting positive parameters...\n";
my %start=();#upper start
my %end=();#upper end
readspe();
my @spe=keys(%spe);
push @spe,keys(%speg);
%fh=();
foreach (@spe){
	$spe=$_;
	open ($fh{$spe},">./training_parameters_$version/$spe.pmt");
}
open ($fh{"overall"},">./training_parameters_$version/overall.pmt");
open (IN,"./experimental/miRBase_$version\_experimental_overall.pmt");
while (<IN>){
	$_=~s/\n//;
	if ($_=~/^#/){
		my $title=$_;
		foreach (@spe){
			print {$fh{$_}} "$title\n";
		}
		print {$fh{"overall"}} "$title\n";
	}
	else{
		#if ($_ =~ /ath-MIR156a/){
		#	print "check";
		#}
		my @tmp=split("\t",$_);
		my @dat=split("_",$tmp[0]);
		my $id=lc($dat[0]);
		@dat=split("-",$id);
		my $spe=$dat[0];
		if (exists $mature{$tmp[2]}){
			if ($mature{$tmp[2]}=~/$id/){
				print "Extracting positive parameters for $id\n";
				$start{$id}=$tmp[74];
				$end{$id}=$tmp[75];
				print {$fh{"overall"}} "$_\n";
				print {$fh{$spe}} "$_\n";
				print {$fh{$spe{$spe}}} "$_\n";
			}
		}
	}
}
close IN;
my @fh=keys (%fh);
foreach (@fh){
	close $fh{$_};
}
print "Done\n";
#create the random sequences
my $i;
my @data=();
my $mark="";
@spe=keys(%spe);
if ($spename ne "all"){
	@spe=();
	if (not exists $spe{$spename}){
		die "Error: unrecognized species name: $spename\n";
	}
	push @spe,$spename;
}
%fh=();
#foreach (@spe){
#	extractpara("1",$_);
#}
my $run;
my %thr;
for ($run=0;$run<@spe;$run=$run+$cores){
	my %thr=();
	for ($i=0;$i<$cores;$i++){
		if (($run+$i)<@spe){
			#my $tmp=$run+$i;
			#print "$tmp:$spe[$tmp]\n";
			$thr{$run+$i} = threads->new("extractpara",$run+$i,$spe[$run+$i]);
		}
	}
	my @run=keys (%thr);
	foreach (@run){
		$thr{$_} ->join;
	}
	runningcheck1:
	my @running = threads->list(threads::running);
	if (@running>0){
		my $tmp=@running;
		print "Waiting for $tmp unfinished threads...\n";
		sleep(3);
		goto runningcheck1;
	}
	my @joinable = threads->list(threads::joinable);
	if (@joinable>0){
		foreach (@joinable){
			$_->join();
		}
	}
}
################################################################################
#step4	model training
################################################################################
step4:
if ($step !~/,4,/){
	goto stepend;
}
$out="models_$version";
if (-e $out){
	check:
	print "Warning: $out is exist, do you want to remove the folder? (yes/no)";
	my $check=<>;
	$check=~s/\n//;
	if (lc($check) eq "yes" or $check eq ""){
		print "removing the $out folder...\n";
		system "rm -rf $out";
		system "mkdir $out";
	}
	elsif (lc($check) eq "no"){
		print "The same model will be overwritten!\n";
	}
	else{
		goto check;
	}
}
else{
	system "mkdir $out";
}
if (not -e "training_parameters_$version"){
	die "Error: not found training_parameters_$version\n";
}
#read species ifnor
readspe();
#collect the species that going to be created
print "collect the species that going to be created...\n";
my @group;
my %group;
@files=glob("./training_parameters_$version/*.pmt");
foreach (@files){
	#next if /random/;
	my @tmp;
	open (IN,$_);
	@tmp=<IN>;
	close IN;
	if (@tmp>1){
		my $group=$_;
		$group=~s/\.\/training_parameters_$version\///;
		$group=~s/\.pmt//;
		my @dat=split("_",$group);
		$group=$dat[0];
		$group{$group}+=1;
		#push @group,$group;
	}
}
#quick train
#delete $group{"Metazoa"};
#delete $group{"overall"};
#delete $group{"Viridiplantae"};

my @tmp=keys (%group);
foreach (@tmp){
	if ($group{$_} eq 21){#total number of files, include one positive
		push @group,$_;
	}
}
if ($spename ne "all"){
	@group=();
	push @group,$spename;
}
#get the species belonging
print "get the species belonging...\n";
$spe{"overall"}="Overall";
#start training
our @svmpmt;
foreach (@group){
	my $group=$_;
	#$group="Viruses";
	print "trainging models for $group\n";
	my @data=();
	if (exists $spe{$group}){
		@svmpmt=pmt($spe{$group});
	}
	else{
		@svmpmt=pmt($speg{$group});
	}
	if (@svmpmt eq 0){
		@svmpmt=pmt($spe{"overall"});
	}
	#start training
	#training($group,6);
	#for (my $i=1;$i<=20;$i++){
	#	training("Viridiplantae",$i);
	#}
	#$group="Viridiplantae";
	for ($run=0;$run<20;$run=$run+$cores){
		my %thr=();
		for ($i=0;$i<$cores;$i++){
			if (($run+$i)<20){
				my $level=$run+$i+1;
				$thr{$run} = threads->new("training",$group,$level);
			}
		}
		my @run=keys (%thr);
		foreach (@run){
			$thr{$_} ->join;
		}
		runningcheck2:
		my @running = threads->list(threads::running);
		if (@running>0){
			my $tmp=@running;
			print "Waiting for $tmp unfinished threads...\n";
			sleep(3);
			goto runningcheck2;
		}
		my @joinable = threads->list(threads::joinable);
		if (@joinable>0){
			foreach (@joinable){
				$_->join();
			}
		}
	}
}
################################################################################
sub extractpara{
	#the script will extract the positive and the negative paraeters for one species
	my ($fh,$spe)=@_;
	my %fh;
	$fh="pmt$fh";
	open ($fh,"./experimental/miRBase_$version\_experimental_overall.pmt");
	while (<$fh>){
		$_=~s/\n//;
		if ($_=~/^#/){
			$title=$_;
			for ($i=1;$i<=20;$i++){
				open ($fh{"$spe\_$i"},">./training_parameters_$version/$spe\_random_$i.pmt");
				print {$fh{"$spe\_$i"}} "$title\n";
				if (not -e "./training_parameters_$version/overall\_random_$i.pmt"){
					open ($fh{"overall\_$i"},">./training_parameters_$version/overall\_random_$i.pmt");
					print {$fh{"overall\_$i"}} "$title\n";
					$fh{"overall_check_$i"}=1;
				}
				if (not -e "./training_parameters_$version/$spe{$spe}\_random_$i.pmt"){
					open ($fh{"$spe{$spe}\_$i"},">./training_parameters_$version/$spe{$spe}\_random_$i.pmt");
					print {$fh{"$spe{$spe}\_$i"}} "$title\n";
					$fh{"$spe{$spe}_check_$i"}=1;
				}
			}
		}
		else{
			my @tmp=split("\t",$_);
			my @dat=split("_",$tmp[0]);
			my $id=$dat[0];
			#print "Calculating negative for $id...\n";
			#if ($id =~ /ath/){
			#	print "check";
			#}
			if ($id ne $mark and $mark ne ""){
				if ($mark=~/$spe/){
					print "Extracting negative parameters for $mark\n";
				}
				for ($i=1;$i<=20;$i++){
					my $tmp=@data;
					my @dat=split("-",$mark);
					my $s=$dat[0];
					#if ($spe{$s} eq "Viridiplantae"){
					#	print "check\n";
					#}
					#if the totally pmt is less than the levels
					if ($tmp<=$i){
						my @tmp=@data;
						foreach (@tmp){
							if ($s eq $spe){
								print {$fh{"$spe\_$i"}} "$_\n";
							}
							if (exists $fh{"overall_check_$i"}){
								if ($fh{"overall_check_$i"} eq 1){
									print {$fh{"overall\_$i"}} "$_\n";
								}
							}
							if (exists $fh{"$spe{$s}_check_$i"}){
								if ($fh{"$spe{$s}_check_$i"} eq 1){
									print {$fh{"$spe{$s}\_$i"}} "$_\n";
								}
							}
						}
					}
					else{#the the pmt rows is more than levels, get one from random
						my %random=();
						my $j;
						my $t=0;
						for ($j=1;$j<=$i;$j++){
							getrand:
							$t++;
							my $rand=int(rand($tmp));#need to check here
							if (not exists $random{$rand}){
								$random{$rand}=1;
								if ($s eq $spe){
									print {$fh{"$spe\_$i"}} "$data[$rand]\n";#need to check here
								}
								if (exists $fh{"overall_check_$i"}){
									if ($fh{"overall_check_$i"} eq 1){
										print {$fh{"overall\_$i"}} "$data[$rand]\n";
									}
								}
								if (exists $fh{"$spe{$s}_check_$i"} ){
									if ($fh{"$spe{$s}_check_$i"} eq 1){
										print {$fh{"$spe{$s}\_$i"}} "$data[$rand]\n";
									}
								}
							}
							else{
								if ($t>500){
									die "Warning: too much runs for $s at level $i, please check\n";
								}
								goto getrand;
							}
						}
					}
				}
				@data=();
			}
			if (exists $start{lc($id)}){
				if (abs($tmp[74]-$start{lc($id)})>5 and abs($tmp[75]-$end{lc($id)})>5){
					push @data,$_;
				}
			}
			$mark=$id;
		}
	}
	close $fh;
	for ($i=1;$i<=20;$i++){
		if (exists $fh{"$spe\_$i"}){
			close $fh{"$spe\_$i"};
		}
		if (exists $fh{"overall_check_$i"}){
			if ($fh{"overall_check_$i"} eq 1){
				close $fh{"overall\_$i"};
				$fh{"overall_check_$i"}=0;
			}
		}
		if (exists $fh{"$spe{$spe}_check_$i"}){
			if ($fh{"$spe{$spe}_check_$i"} eq 1){
				close $fh{"$spe{$spe}\_$i"};
				$fh{"$spe{$spe}_check_$i"}=0;
			}
		}
	}
}
################################################################################
sub training{
	#will receive one species one group, create the data, do the training
	my ($group,$level)=@_;
	my $trainingset="./models_$version/$group\_$level.txt";
	my $modelfile="./models_$version/$group\_$level.model";
	#list the parameters
	my %pmt;
	@tmp=pmt('mirbase');
	my $tmp=0;
	foreach (@tmp){
		$pmt{$_}=$tmp;
		$tmp+=1;
	}
	#read the positive data
	my $fhin="IN".$level;
	my $fhout="OUT".$level;
	open ($fhout,">$trainingset");
    open ($fhin,"training_parameters_$version/$group.pmt");
	while (<$fhin>){
		next if /^#/;
		$_=~s/\n//g;
		my @tmp=split("\t",$_);
		if (@tmp eq pmt('mirbase')){
			print $fhout "1\t";
			my $count=1;
			foreach (@svmpmt){
				my $value;
				if (lc($tmp[$pmt{$_}]) eq "null" or $_ eq "firstbase" or $_ eq "penultimateposition" or $_ eq "terminalnucleotide"){
					$value=nt2number($tmp[$pmt{$_}]);
				}
				else{
					$value=$tmp[$pmt{$_}];
				}
				print $fhout "$count:$value\t";
				$count++;
			}
			print $fhout "\n";
		}
	}
	close $fhin;
	open ($fhin,"training_parameters_$version/$group\_random\_$level.pmt");
	while (<$fhin>){
		next if /^#/;
		$_=~s/\n//g;
		my @tmp=split("\t",$_);
		if (@tmp eq pmt('mirbase')){
			print $fhout "0\t";
			my $count=1;
			foreach (@svmpmt){
				my $value;
				if (lc($tmp[$pmt{$_}]) eq "null" or $_ eq "firstbase" or $_ eq "penultimateposition" or $_ eq "terminalnucleotide"){
					$value=nt2number($tmp[$pmt{$_}]);
				}
				else{
					$value=$tmp[$pmt{$_}];
				}
				print $fhout "$count:$value\t";
				$count++;
			}
			print $fhout "\n";
		}
	}
	close $fhin;
	close $fhout;
	system "svm-train -s 0 -t 2 -b 1 -q $trainingset $modelfile";
	system "rm $trainingset";
}
################################################################################
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
   if ($class eq "mirbase"){
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
   elsif ($class eq "Overall"){
        @parameter=("unpairedrate_lowerstem","prelength","internalloopnumber_lowerstem","length_upperstem",
        "miinternalloop","firstbase","mintcontent_a","migc","pregc","prentcontent_u",
        "prentcontent_a","internalloop_topstem","preunpairedrate","mipairs","prepairs",
        "internalloopnumber_topstem","unpairedrate_topstem","mintcontent_c","mistart",
        "miunpairedrate","mintcontent_g","terminalnucleotide","prentcontent_c",
        "prentcontent_g","mintcontent_u");#25pmt
   }
#   elsif (lc($class) eq "svmanimal"){
#        @parameter=("internalloop_topstem","internalloopnumber_topstem","length_topstem","length_upperstem",
#        "migc","miinternalloop","miinternalloopnumber","mintcontent_a","mintcontent_c","mintcontent_g",
#        "mintcontent_u","mistart","miunpairedrate","penultimateposition","pregc","prelength",
#        "prentcontent_a","prentcontent_c","prentcontent_g","prentcontent_u","preunpairedrate",
#        "stability","unpairedrate_lowerstem","unpairedrate_topstem");#24pmt
#   }
   elsif ($class eq "Metazoa"){
        @parameter=("internalloop_topstem","internalloopnumber_topstem","length_topstem","length_upperstem",
        "migc","miinternalloop","miinternalloopnumber","mintcontent_a","mintcontent_c","mintcontent_g",
        "mintcontent_u","mistart","miunpairedrate","penultimateposition","pregc","prelength",
        "prentcontent_a","prentcontent_c","prentcontent_g","prentcontent_u","preunpairedrate",
        "stability","unpairedrate_topstem");#24pmt
   }
   elsif ($class eq "Viridiplantae"){#Viridiplantae
        @parameter=("firstbase","internalloop_topstem","internalloop_lowerstem","internalloopnumber_lowerstem",
        "length_upperstem","migc","migu","miinternalloop","miinternalloopnumber","mintcontent_a",
        "mintcontent_g","mintcontent_u","mipairs","penultimateposition","pregc","prentcontent_a",
        "prentcontent_c","prentcontent_g","prentcontent_u","preunpairedrate","stability",
        "unpairedrate_lowerstem","unpairedrate_topstem","upperstart");#24pmt
	}
   else{
	  @parameter=();
   }
   return @parameter;
}
################################################################################
sub nt2number{
   #>>>>>>>>>>>>>>>>>>>>>>>>>>INTRODUCTION<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   #===========The nt2number was used to translate the nt to ACC number
   #===========Usage: nt2number(nt)
   #===========The number of acc will be return
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   my ($value)=@_;
   $value=lc($value);
   if ($value eq "null"){
	  $value="-1";
   }
   elsif ($value eq ""){
	  $value=0;#if nothing
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
################################################################################
sub readspe{
	#the script will read species informations
	print "Preparing the species information from organisms.txt...";
	if (not -e "./miRBase/$version/organisms.txt"){
		system "gunzip ./miRBase/$version/organisms.txt.gz";
	}
	if (not -e "./miRBase/$version/organisms.txt"){
		system "mv ./miRBase/$version/organism.txt ./miRBase/$version/organisms.txt";
	}
	open (IN,"./miRBase/$version/organisms.txt") or die "Error: Not found organisms.txt\n";
	while (<IN>){
		$_=~s/\n//;
		my @tmp=split("\t",$_);
		my @dat=split(";",$tmp[@tmp-1]);
		$spe{$tmp[0]}=$dat[0];
		$speg{$dat[0]}=$dat[0];
	}
	close IN;
	print "Done\n";
}
stepend:
