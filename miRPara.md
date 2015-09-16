
---

# INSTALL #

---


miRPara was written in perl and developed under Linux. We suggest to run miRPara under linux. But it should be able to run under Mac OS and Windows if all third-party software packages were installed.

Please Cite:

MiRPara: a SVM-based software tool for prediction of most probable microRNA coding regions in genome scale sequences. Wu Y., Wei B., Liu H., Li T., Rayner S. BMC Bioinformatics. 2011 Apr 19; 12(1):107

**Required Packages**

_Version 6.0 or above_
  * Getopt::Long (CPAN)
  * threads (CPAN)
  * threads::shared (CPAN)
  * Cwd (CPAN)
  * File::chdir (CPAN)
  * UNAFold (<a href='http://mfold.rna.albany.edu/?q=DINAMelt/software'>download</a>)
  * ct2out (<a href='http://code.google.com/p/mirpara/downloads/detail?name=ct2out.tar.gz&can=2&q='>download</a>)
  * libsvm (http://www.csie.ntu.edu.tw/~cjlin/libsvm/)

Copy UNAFold.pl, ct2out, svm-predict to $PATH
```
#Install the packages from CPAN
#See the instruction at http://perl.about.com/od/packagesmodules/qt/perlcpan.htm

#Open a terminal under linux, type the following to install CPAN packages. (Ignore the "$")
$ perl -MCPAN -e shell
#for ubuntu root, type
$ sudo perl -MCPAN -e shell
#install CPAN packages
cpan> install Getopt::Long
cpan> install threads
cpan> install threads::shared
cpan> install Cwd
cpan> install File::chdir

#Install UNAFold
#See the INSTALL in the UNAFold package
#Check the success of installation by typing
$ UNAFold.pl
# if returns "UNAFold.pl: command not found", do the following under the UNAFold.pl folder and test again, you might need to be root to do so.
$ sudo cp UNAFold.pl /use/bin/

#Install ct2out
$ gfortran ct2out.f -o ct2out
#or
$g77 -o ct2out ct2out.f
#then copy ct2out to PATH, you might need to be root to do so
$ sudo cp ct2out /usr/bin/
#Test the ct2out by
$ ct2out
# Failed if returns "ct2out: command not found", success if returns nothing

#Install libsvm
#See the README in the libsvm package
#copy svm-predict to the path, you might need to be root to do so
$ sudo cp svm-predict /usr/bin/

```

_Version 5.3 or below_
  * Getopt::Long
  * Algorithm::SVM (<a href='http://code.google.com/p/mirpara/downloads/detail?name=algorithm_SVM_install_fix.pdf&can=2&q='>click here if you have install problem</a>)
  * Cwd
  * UNAFold (<a href='http://mfold.rna.albany.edu/?q=DINAMelt/software'>download</a>)
  * ct2out (<a href='http://code.google.com/p/mirpara/downloads/detail?name=ct2out.tar.gz&can=2&q='>download</a>)

Copy UNAFold.pl, ct2out to $PATH


---

# USAGE #

---


**Before Run**

A folder called "models" should be created and placed in the same directory as miRPara.pl, all the models should be copied into this folder. Or just created a link to the models folder and rename the link as "models". (It is not required if the models was provided by "-s" tag, only for version 6.0 or above)

miRPara.pl does not required to be copied to $PATH. If you did so, do copy the model folder links into $PATH.

The _mature.fa_ and _organisms.txt_ was required for the miRPara.pl 6.0 or above. miRPara.pl will check the $PATH and automatically download the files if they were not exist.

The selection of miRPara\_model\_trainer.pl have to match the miRPara.pl version:
|Model Trainer Version|miRPara Version|
|:--------------------|:--------------|
|3.0 or above         |6.0 or above   |
|2.3 or below         |5.3 or below   |
The previous models do not fit for miRPara.pl 6.0 or above!

**miRPara version 6.0 or above**

```
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
-t, --threads=<No. of threads> (defaults to 1)
--pmt, --Only calculate the parameters without further prediction

File (one of following):
*.fasta, --Only fasta format sequences were allowed
*.pmt, --Repredict from the parameter files

Report bugs to Yonggan Wu (weasteam@gmail.com) or Dr. Simon Rayner (raynere@wh.iov.cn)

Homepage: http://www.whiov.ac.cn/bioinformatics/mirpara
Google Project: http://code.google.com/p/mirpara/
Facebook: https://www.facebook.com/mirpara2009
```

_Examples:_

1) check the version and help message:
```
perl miRPara.pl -v or perl miRPara.pl --version
```

2) print the help message:
```
perl miRPara.pl -h or perl miRPara.pl --help
```

3) Basic usage:
```
perl miRPara.pl test.fasta
```

4) Advanced usage by providing species, cutoff, model level, prilengthlimit, overloplength
```
perl miRPara.pl -s animal -c 0.5 -l 12 -p 70 -o 100 test.fasta
```

5) Advanced usage by using species group model
```
perl miRPara.pl -s animal test.fasta
perl miRPara.pl -s Metazoa test.fasta
```

6) Advanced usage by using specific species model
```
perl miRPara.pl -s hsa test.fasta
```

7) Advanced usage by providing a model file
Instead of providing species short names, a trained model can be provided directly from "-s" tag. The level tag "-l" will be ignored:
```
perl miRPara.pl -s /home/testuser/mirpara/models/hsa_12.model test.fasta
```

8) Advanced usage by using a custom trained models
```
perl miRPara.pl -s customModelName test.fasta
perl miRPara.pl -s /home/testuser/mirpara/models/hocustomModelName_10.model test.fasta
```

9) Run with multiple cores
```
perl miRPara.pl -t 12 test.fasta
```

10) Only calculate the X.pmt files. skip the prediction process
```
perl miRPara.pl --pmt test.fasta
```

11) Repredict the result by change the cutoff or model level.
```
perl miRPara.pl -c 1 -l 20 test.pmt
```


**miRPara version 5.2 and 5.2 above**

```
Usage: miRPara.pl [Options ] file [File ]

Options:
-v, --version
-h, --help
-n, --name=<abbreviated name of the species>
-s, --species=species group name <overall, animal, plant, virus> (defaults as overall)
              or species short name <hsa, mmu, cel, ...>
              or path of a model file: /home/xxx/miRPara/model/hsa_10.model
-c, --cutoff=<the cutoff to svm prediction probabilities> (defaults to 0.8)
-l, --Levels=<1..20>(defaults to 1)
-p, --prilengthlimit=<limit to the pri-miRNA length> (defaults to 60)
-f, --foldlength=<The length to be split for folding> (defaults to 500)
-o, --overlaplength=<The overlap length of two nearby splited fragments> (defaults to 150)
-t, --cores=<No. of cores> (defaults to 1)
--pmt, --Only calculate the parameters without further prediction

File (one of following):
X.fasta, --Only fasta format sequences were allowed
X.pmt, --Repredict from the parameter files

Homepage: http://www.whiov.ac.cn/bioinformatics/mirpara
Google Project: http://code.google.com/p/mirpara/
Facebook: https://www.facebook.com/mirpara2009
```

_Examples:_

1) check the version and help message:
```
perl miRPara.pl -v or perl miRPara.pl --version
```

2) print the help message:
```
perl miRPara.pl -h or perl miRPara.pl --help
```

3) Basic usage:
```
perl miRPara.pl test.fasta
```

4) Advanced usage by providing the name, species, cutoff, model level, prilengthlimit, overloplength
```
perl miRPara.pl -n test -s animal -c 0.5 -l 12 -p 70 -o 100 test.fasta
```

5) Advanced usage by using species group model
```
perl miRPara.pl -s animal test.fasta
```

6) Advanced usage by using specific species model
```
perl miRPara.pl -s hsa test.fasta
```

7) Advanced usage by providing a model file directly
Note: If you have the miRPara model in a directory different to miRPara.pl program, you can redirect the models folder to miRPara.pl by giving trained models at any level.
-l option will be ignored if full path of a model provided
```
perl miRPara.pl -s /home/testuser/mirpara/models/hsa_12.model test.fasta
```

8) Advanced usage by using a custom trained models
```
perl miRPara.pl -s customModelName test.fasta
perl miRPara.pl -s /home/testuser/mirpara/models/hocustomModelName_10.model test.fasta
```

9) Run with multiple cores
Note: need BIG RAM (>16G, depend on the input data size), if your computer do not have big RAM, you can split the input fasta files into small files and run each of them in one core.
This function is still in beta and can be buggy.
```
perl miRPara.pl -t 12 test.fasta
```

10) Only calculate the X.pmt files. skip the prediction process<br>
<pre><code>perl miRPara.pl --pmt test.fasta<br>
</code></pre>

11) Repredict the result by change the cutoff or model level.<br>
<pre><code>perl miRPara.pl -c 1 -l 20 test.pmt<br>
</code></pre>


<b>miRPara version 5.2 or below</b>

<pre><code>Usage: miRPara.pl [Options ] file [File ]<br>
<br>
Options:<br>
-v, --version<br>
-h, --help<br>
-n, --name=&lt;abbreviated name of the species&gt;<br>
-s, --species=&lt;overall, animal, plant or virus&gt; (defaults as overall)<br>
-c, --cutoff=&lt;the cutoff to svm prediction probabilities&gt; (defaults to 0.8)<br>
-l, --Levels=&lt;1..20&gt;(defaults to 1)<br>
-p, --prilengthlimit=&lt;limit to the pri-miRNA length&gt; (defaults to 60)<br>
-f, --foldlength=&lt;The length to be split for folding&gt; (defaults to 500)<br>
-o, --overlaplength=&lt;The overlap length of two nearby splited fragments&gt; (defaults to 150)<br>
-t, --cores=&lt;No. of cores&gt; (defaults to 1)<br>
--pmt, --Only calculate the parameters without further prediction<br>
<br>
File (one of following):<br>
X.fasta, --Only fasta format sequences were allowed<br>
X.pmt, --Repredict from the parameter files<br>
</code></pre>

<i>Examples:</i>

1) check the version and help message:<br>
<pre><code>perl miRPara.pl -v or perl miRPara.pl --version<br>
</code></pre>

2) print the help message:<br>
<pre><code>perl miRPara.pl -h or perl miRPara.pl --help<br>
</code></pre>

3) Basic usage:<br>
<pre><code>perl miRPara.pl test.fasta<br>
</code></pre>

4) Advanced usage by providing the name, species, cutoff, model level, prilengthlimit, overloplength<br>
<pre><code>perl miRPara.pl -n test -s animal -c 0.5 -l 12 -p 70 -o 100 test.fasta<br>
</code></pre>

5) Run with multiple cores<br>
Note: need BIG RAM (>16G, depend on the input data size), if your computer do not have big RAM, you can split the input fasta files into small files and run each of them in one core.<br>
<pre><code>perl miRPara.pl -t 12 test.fasta<br>
</code></pre>

6) Only calculate the X.pmt files. skip the prediction process<br>
<pre><code>perl miRPara.pl --pmt test.fasta<br>
</code></pre>

7) Repredict the result by change the cutoff or model level.<br>
<pre><code>perl miRPara.pl -c 1 -l 20 test.pmt<br>
</code></pre>

<hr />
<h1>Result</h1>
<hr />

The predicted result named as "＊.out". It is a tab-delimited text file with following columns:<br>
<br>
<table><thead><th>Column Name</th><th>Example 1</th><th>Example 2</th><th>Description</th></thead><tbody>
<tr><td>priid      </td><td>lin:1-94 </td><td>lin:1-94 </td><td><fasta name>:<start position>-<end position></td></tr>
<tr><td>priseq     </td><td>augcuuccggcCUGUUCCCUGAGACCUCAAGU<br>gugaguguacuauugaugcuucacaccugggcucuccggg<br>uaccaggacgguuugagcagau</td><td>augcuuccggccuguUCCCUGAGACCUCAAGU<br>GUGaguguacuauugaugcuucacaccuggg<br>cucuccggguaccaggacgguuugagcagau</td><td>Upper case indicated the mature miRNA sequences</td></tr>
<tr><td>miid       </td><td>lin:1-94:12_21</td><td>lin:1-94:16_20</td><td>< priid >:<miRNA start>underline<miRNA length></td></tr>
<tr><td>miseq      </td><td>CUGUUCCCUGAGACCUCAAGU</td><td>UCCCUGAGACCUCAAGUGUG</td><td>miRNA sequence</td></tr>
<tr><td>strand     </td><td>5P       </td><td>5P       </td><td>location of the miRNA</td></tr>
<tr><td>SVM_probability</td><td>0.823885 </td><td>0.984632 </td><td>SVM Probability, the higher the better</td></tr>
<tr><td>miRBase    </td><td>NA       </td><td>lin-4-5p </td><td>Hit in miRBase</td></tr></tbody></table>

<hr />
<h1>Release</h1>
<hr />

<b>2013-03-13 miRPara 6.3</b>
<ul><li>fixed a bug in reading the fasta files in ./ folder<br>
</li><li>fixed a bug in outputing wrong format of pmt file when using multicore<br>
</li><li>added a feature to display how much time left to finish</li></ul>

2013-03-13 miRPara 6.2<br>
<ul><li>fix a bug in using the tag "-s animal", thanks to Anurag Chaturvedi for reporting the bug.</li></ul>

2012-11-24 miRPara 6.1<br>
<ul><li>optimize the multicore running, force to wait for unfinished cores.</li></ul>

2012-11-23 miRPara 6.0<br>
<ul><li>replace Algorithm::SVM with libsvm<br>
</li><li>remove "-n" tag<br>
</li><li>multi-core support<br>
</li><li>optimize the known miRNA data hash<br>
</li><li>the specific species support</li></ul>

2012-10-24 miRPara 5.3<br>
<ul><li>fix the bug to extract model level from a full path model file</li></ul>

2012-10-24 miRPara 5.2<br>
<ul><li>add support for custom models</li></ul>

2012-10-02 miRPara 5.1<br>
<ul><li>add facebook and google code page links</li></ul>

2012-09-10 miRPara 5.0<br>
<ul><li>multicore support (beta2)</li></ul>

2012-08-09 miRPara 4.2<br>
<ul><li>fix bug <a href='Use.md'>of uninitialized value $loop in lc at miRPara.pl line 813.</a>
</li><li>optimize the UNAFold.pl & ct2out check<br>
</li><li>multicore support (beta1)</li></ul>

2011-06-02 miRPara 4.1<br>
<ul><li>Correct share path problem</li></ul>

2010-11-15 miRPara 4.0<br>
<ul><li>return probability<br>
</li><li>Optimize the codes<br>
</li><li>Change the output format</li></ul>

2010-04-06 miRPara 3.0<br>
<ul><li>Changed the predicted models, replace species (animal, plant, virus) with miRNA families<br>
</li><li>new algorithm to calculate parameters<br>
</li><li>three different models (numerical, arphabetic and bulge)</li></ul>

2009-11-28 miRPara 2.0<br>
<ul><li>Retrain all the models<br>
</li><li>Different models been added, the miRNA could been predicted with specific models<br>
</li><li>Changed the SVM parameters for different models<br>
</li><li>Optimize the scripts</li></ul>

2009-9-16 miRPara 1.7<br>
<ul><li>Update the sub script of callunafold</li></ul>

2009-9-2 miRPara 1.6<br>
<ul><li>Retrained all the data with UNAFold.pl predicted structures.<br>
</li><li>Change the output alignment form, more easy to read the output results<br>
</li><li>Increse the level of models to 20<br>
</li><li>Added the enviroment, more easy to install</li></ul>

2009-5-12 miRPara 1.5<br>
<ul><li>Replaced the Bio::Seq and Bio::SeqIO module, again increased the speed<br>
</li><li>Change the input data read policy, allow multiple input sequences</li></ul>

2009-4-18 miRPara 1.4<br>
<ul><li>Retrain the models files, exclude the bad format files<br>
</li><li>Reorganize all miRNA parameters<br>
</li><li>Optimize output format</li></ul>

2009-4-17 miRPara 1.3<br>
<ul><li>Replace the XML module, increase the analysis speed.<br>
</li><li>Update the file output format of XXX.pmt</li></ul>

2009-3-9 miRPara 1.1<br>
<ul><li>Update the algorithm to extract UNAFold.pl exported data<br>
</li><li>Correct the candidates extraction algorithms</li></ul>

2009-2-19 miRPara 1.0 Beta<br>
<ul><li>All the script were rewrite with perl, and could be run under Linux</li></ul>

2008-9 miRPara 0.5<br>
<ul><li>miRPara Windows application software, write in Visual Basic</li></ul>

2007-12 miRPara 0.1<br>
<ul><li>Based on Microsoft Excel, the first miRPara version released<br>
</li><li>It calculate the parameters and filter based on artificial cutoffs<br>
</li><li>Write in VBA