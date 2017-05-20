#!/usr/bin/env perl -w

#Description: This perl script was meant to produce shell scripts for analysis the whole genome sequencing data based on provided parameters which call them use-help|? flag.
#
#WHOLE PIPELINE:raw data -> filter -> align -> correct bam -> call SNPs and indels -> call SVs 
#
#__Author__: Songjing ,__E-mail__: songjing@genomics.cn
#__Edit__: 14 Mar,2017
#__Modified__: -

use strict;
use File::Path;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use FindBin '$Bin';

sub usage {
	print STDERR << "USAGE";
usage: perl $0 [options]
options: 
	-c:*	FILE	configure file,including the path of tools used.(eg:HiseqWGS_2017.conf)(required)
	-f:*	FILE	sample fastq list(required)
	-R:	STR	read group header line such as "\@RG\\tID:id\\tSM:samplename\\tLB:lib"[Default: Auto]
	-m:	NUM	memory to use [Default:12 (Gb)]
	-p:	Boolens first fastq file consists of interleaved paired-end sequences[Default: False]
	-t:	INT	number of threads to use [default:16]
	-o:	STR	output dir prefix  [Defalut:'./outdir_\$date']
	-tmp:	STR	TMP dir prefix     [Default:'./TMP_\$date']
	
	-A:*	Boolens Call align(align+markdup+sort+index) [Default:False]	
	-V:*	Boolens Call var  (calling SNP,indels)       [Default:False]
	-S:*	Boolens Call sv   (calling sv,CNV)           [Default:False]
	
	-trio:	Boolens Necessary if your samples came from a family(eg:child,father and mother)[Defalut: False]
	-z:	Boolens Necessary if your samples are tumor/normal pairs[Default: False]	
	-muTec2 Boolens use muTect2(beta version) to call SNVs and indels in tumor/normal samples[Default:false. Will use muTect-1.1.7] 

	-help|?:	show this usage page
	-more:		more default information
 eg:
	1) For single/multiple sample(s):
	perl $0 -c HiseqWGS_2017.conf -f rawData.list -A -V -S
	2) For tumor/normal pair(s) samples:
	perl $0 -c HiseqWGS_2017.conf -f rawData.list -A -z -S
USAGE
	&details;
	exit 1;
} 

sub details {
	print STDOUT << "DETAIL";
Description: WGS Analysis pipeline: Filter -> Align -> correct bam -> Call variants ->Variants annotation
Author: Songjing ,E-mail: "songjing\@genomics.cn"
Date: 20170314
Modified: -

	How to config rawdata.list:

	1).if input sequencing data was single/multiple samples :(config like this)
	NA12878 lib1 NA12878 /home/songjing/HiseqWGS_2017/sample/NA12878_1.fq.gz,/home/songjing/HiseqWGS_2017/sample/NA12878_2.fq.gz 100:100 178 20465
	NA12878 lib2 NA12879 /home/songjing/HiseqWGS_2017/sample/NA12879_1.fq.gz,/home/songjing/HiseqWGS_2017/sample/NA12879_2.fq.gz 100:100 178 20465
	[Tips: column 1: project name
	       column 2: library name
	       column 3: sample name
	       column 4: data path of pair-end sequencing data
	       column 5: read length
	       column 6: insert fengzhi
	       column 7: total bases  
	please be sure that column 1 of each line was different. And column 3 = cloumn 1 was the prefix of sequencing data. 
	(eg:NA12878_1.fq.gz,prefix: NA12878)]
   Command: perl $0 -c HiseqWGS_2017.conf -f rawdata.list -A -V -S
	
	2).if input sequencing data was tumor/normal pair:(config like this)
	TCGA-B6-A0I6.tumor TCGA-B6-A0I6 TCGA-B6-A0I6.tumor /home/songjing/HiseqWGS_2017/TCGA-B6-A0I6.tumor_1.fq.gz,/home/songjing/HiseqWGS_2017/sample/TCGA-B6-A0I6.tumor_2.fq.gz 100:100 178 20465 
	TCGA-B6-A0I6.normal TCGA-B6-A0I6 TCGA-B6-A0I6.normal /home/songjing/HiseqWGS_2017/TCGA-B6-A0I6.normal_1.fq.gz,/home/songjing/HiseqWGS_2017/sample/TCGA-B6-A0I6.normal_2.fq.gz 100:100 178 20465
	[Tips: column 1~7 was same as below.please be sure that column 1 = column 3, column 2 = prefix of column 1.
	(eg:TCGA-B6-A0I6.tumor_1.fq.gz,column 1=TCGA-B6-A0I6.tumor,column 2=TCGA-B6-A0I6,column 3=TCGA-B6-A0I6.tumor)]
   Command: perl $0 -c HiseqWGS_2017.conf -f rawdata.list -A -z -S
	The config file must contain these information:
	'reference','include_bed','exclude_bed','seqType','soapnuke','speedseq','GATK','java','samtools','Genome','millsb37','dbsnpVCFb37','phaseb37','1000Gb37','omnib37','hapmapb37','leftAligned','java7','muTect'
DETAIL
	exit 0;
}

my ($config, $fqlist, $read_header, $memory, $help,$somatic,$pair_in_one,$threads,$outdir,$tmpdir,$speedseq_align,$speedseq_var,$use_sv,$trio,$mutect2,$helpmes,$default_settings);
GetOptions(
        "c=s" => \$config,
        "f=s" => \$fqlist,
        "R:s" => \$read_header,
	"m:i" => \$memory,
	"z+" => \$somatic,
	"p+" => \$pair_in_one,
	"t:i" => \$threads,
	"o:s" => \$outdir,
	"tmp:s" => \$tmpdir,
	"A+" => \$speedseq_align,
	"V+" => \$speedseq_var,
	"S+" => \$use_sv,
	"trio+" => \$trio,
	"muTec2+" => \$mutect2,
        "help|?" => \$helpmes,
	"more" => \$default_settings,
);
#die &usage() if (!defined $config || $helpmes);
#die &usage() if (!defined $fqlist || $helpmes);
if((!defined $config) || (!defined $fqlist)){
	if(!defined $default_settings){
		print "\nFatal error :Too less parameters.\n";
		&usage();
		die;
	}else{
		&details();
		die}
}

###setting Default###
#-z 0 -p 0 -t 16 -o outdir/ -T TMP_time -A 0 -V 0 -S 0 -trio 0 -muTec2 0
chomp(my $time=`date +%F-%H-%M`);
print "[INFO] Time: $time\n";
#if(!defined $read_header){
#	$read_header = "\@RG\tID:NA12877.Sample1\tSM:NA12877\tLB:lib1";}
if(!defined $memory){
	$memory = 12;}
if(!defined $somatic){
	$somatic = 0;}
if(!defined $pair_in_one){
	$pair_in_one = 0;}
if(!defined $threads){
	$threads = 16;}
if(!defined $outdir){
	$outdir = "./outdir";}
if(!defined $tmpdir){
	$tmpdir = "./TMP_$time";}
if(!defined $speedseq_align){
	$speedseq_align = 0;}
if(!defined $speedseq_var){
        $speedseq_var =0;}
if(!defined $use_sv){
        $use_sv =0;}
if(!defined $trio){
        $trio =0;}
if(!defined $mutect2){
	$mutect2 =0;}

# %info : path of tools.
# %hash : path of every fastq file in each  project(eg: patient1,patient2,...).
my (%info, %hash, @dependent);
# $cmd : shell running command.
my $cmd = "";

###---read configure file (-c)###
##---for path of corresponding tool###
#print $speedseq_align,$speedseq_var,$use_sv,$somatic,$pair_in_one;
#print $time ;
open IN,"$config" or die $!;
while (<IN>) {
	chomp;
	#ignore the blank line and lines including '#'
	next if ((/^\s*$/) || (/^\s*\#/));
	#delete space of line beginning and ending
	$_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
	if(/^(\w+) \s* = \s* (.*)$/x){
		next if ($2 =~ /^\s*$/);
		$info{$1} = $2;
		#delete space in keys
		$info{$1} =~ s/\s*//;
		next;
	}
}
print "[INFO] read config file....Done.\n";
close IN;

sub readlength{
        my $fqfile =$_[0];
        my $len = 0;
        if ($fqfile =~ /.gz$/){
                open RELE,"zcat $fqfile|" or die $!;
        }else{
                open RELE,"$fqfile" or die $!;
        }
        while(<RELE>){
                last if ($. > 5);
                if($. == 2){
                        chomp($_);
                        $len = length($_);
                }
        }
	close RELE;
        return $len;
}

sub total_base{
        my $fqfile =$_[0];
        my $reads_num = `cat $fqfile |wc -l`;
        $reads_num =$reads_num/4;
        my $basenum = &readlength($fqfile);
        $basenum =$basenum * $reads_num;
        return $basenum;
}

### pick only one or add one ###
### read rawData list ###
my (%sample, %header, %fqgz ,%finished, %unfinished, %report, %sbase, %dependent);
open FQ,"$fqlist" or die $!;
while (<FQ>) {
	chomp;
	next if ((/^\s*$/) || (/^\s*\#/));
	s/^\s*//;
        s/\s*$//;
	my @tmp = split /\s+/;
	my @fq = split /,/,$tmp[3];
	if ($tmp[0] =~ /^\s*\#/) {
                $tmp[0] =~ s/^\s*\#+//;
                $finished{$tmp[2]} = 1;
        }else {
                $unfinished{$tmp[0]} = 1;
        }
        print "[Warnning]: $tmp[0] format seems unfriendly!\n" if(@tmp != 7);
	$fqgz{$tmp[0]}->[0] = $fq[0];
	$fqgz{$tmp[0]}->[1] = $fq[1];
	# use in filter
	#path of fq1 of sample1
	$sample{$tmp[0]}{$tmp[1]}{$tmp[2]}->[0] = $fq[0];
	#path of fq2 of sample1
	$sample{$tmp[0]}{$tmp[1]}{$tmp[2]}->[1] = $fq[1];
	#my $pl ||= "ILLUMINA";
	#generate PU
	my $pu;
	if($time =~/(\d+)-(\d+)-(\d+)-\d+-\d+/){
		$pu = "$1$2$3_$tmp[0]_L1_$tmp[2]";}
	print "[INFO] PU(platform unit):$pu\n";
	#used only once in align
	$header{$tmp[0]} = "\"\@RG\\tID:$tmp[0]\\tPL:ILLUMINA\\tPU:$pu\\tSM:$tmp[2]\\tLB:$tmp[1]\"";
	#insert maximumn (unuseful now)
	$sample{$tmp[0]}{$tmp[1]}{$tmp[2]}->[2] = 300;
	# total bases (unuseful now)
        $sample{$tmp[0]}{$tmp[1]}{$tmp[2]}->[3] = &total_base($fq[0]);
	#read length 
        $report{$tmp[2]} = &readlength($fq[0]);
	my $eee = $report{$tmp[2]};
	print "[INFO] Read length:$eee\n";
        $sbase{$tmp[0]} = $sample{$tmp[0]}{$tmp[1]}{$tmp[2]}->[3];
	my $fff = $sbase{$tmp[0]};
	print "[INFO] total bases:$fff per fastq\n";
        $dependent{$tmp[2]} = $tmp[7] if(@tmp >= 8);
	
}
close FQ;
### default values ###
`mkdir -p $outdir`;
if (-e $outdir){
	print STDOUT "[INFO] creat $outdir successfully.\n";}
else{
	die "error: Could not creat directory $outdir [Permission denied].\n";}
$outdir = abs_path($outdir);
$info{'MScheck'} = abs_path($info{'MScheck'}) if(exists $info{'MScheck'});
my $shlog = "$outdir/shell.log";
open LOG, ">$shlog" or die $!;
print LOG "Pipeline process: WGS\n\n";
close LOG;

my $sh_dir = "$outdir/shell";
my $result_dir = "$outdir/result";
my $listdir = "$outdir/list";
my $processdir = "$outdir/process";
for($sh_dir, $processdir, $result_dir, $listdir){
        `mkdir -p $_` unless(-d $_);
}
`chmod -R 777 $outdir`;

### analysis content confirm###
if (($speedseq_align == 0) && ($speedseq_var == 0) && ($use_sv == 0) && ($somatic == 0)){
	print  "error: too less parameters ,it must contain at least one of '-A','-V','-S','-z'.\n";
	&usage;
	die;
}

######### starting analysis ###########
### configure PATH environment ###
my ($path1,$path2,$path3,$path4);
my ($rootpath,$speedseqpa);
$rootpath =`which root`;
die "[Error]:You haven't installed root yet!\n" if ($rootpath eq "");
$rootpath = dirname(dirname $rootpath);
$speedseqpa = dirname($info{'speedseq'}) if (exists $info{'speedseq'});
$rootpath =abs_path($rootpath);
open PA,">$outdir/set_path.sh" or die $!;
print PA "#!/bin/sh\n";
print PA "export SHELL=/bin/bash\n";
print PA "export HOME=./\n";
print PA "export TMPDIR=$outdir/process\n";
print PA "export TMP=$outdir\n";
$path1 = "export PATH=\"\$PATH:$speedseqpa:$rootpath/bin:$rootpath:$rootpath/lib\"\n";
$path2 ="ROOTSYS = \"export ROOTSYS=$rootpath\"\n";
$path3 ="LD_LIBRARY_PATH = \"export LD_LIBRARY_PATH=$rootpath/lib\"\n";

$path4 ="source ~/.bashrc\n";
print PA $path1,$path2,$path3,$path4;
close PA;

#print "$info{'soapnuke'}\n";
### filter ###
#$Bin .= "/opt/bin/soapnuk1.5.3"; 

#my $tile;
#my $fq1,$fq2;
#$fq1 = basename
#if((exists $info{'soapnuke'}) && (exists $info{'soapnuke_conf'})) {
	#$tile = `perl $info{'findNtile'} -fq1 $hash{$project[0]}->[0] -fq2 $hash{$project[0]}->[1]`;
#	$cmd = "$info{'soapnuke'} -c $info{'soapnuke_conf'} -st $info{'seqType'} -f $fqlist -s _all -o $outdir";
#	`echo "generate all filter shell command:\n$cmd\n\nClean Data List File: $fqlist\n" >>$shlog`;
#	`$cmd`;
#	push @dependent, "$listdir/${soft}_all_dependence.txt";
#}else{
#	die "ERROR:There is no filter tool in config file.\n";
#}

my $nofilter = 1;
if((exists $info{'soapnuke'}) && ($nofilter = 0)){
	$Bin = $info{'soapnuke'};
	my $OutDir = $outdir;
	my ($cut, $seqType, $propf, $suffix, $monitorOptions);
	$cut = 0;
	$suffix = '';
	$info{"soapnuke_para"} ||= "-n 0.1 -q 0.5 -i";
	$seqType = $info{'seqType'};
	`echo "SOAPnuke sample lane:" >>$shlog`;

	my %filterSamp;
	my $lanenum = 0;
	open OUT,">$listdir/cleanData$suffix.list" or die $!;
	foreach my $s (keys %sample) {
		`mkdir -p $result_dir/$s/clean_data`;
		foreach my $lib (keys %{$sample{$s}}) {
		     foreach my $lane (keys %{$sample{$s}{$lib}}) {
			`mkdir -p $sh_dir/$lane`;
			my $shell = "$sh_dir/$lane/filter_$lane.sh";
			open SHE,">$shell" or die $!;
			if (exists $finished{$lane}) {
				print OUT "\#$s\t$lib\t$lane\t$result_dir/$s/clean_data/$lane\_1.fq.gz,$result_dir/$s/clean_data/$lane\_2.fq.gz\t$report{$lane}\t$sample{$s}{$lib}{$lane}->[2]\t$shell:$memory\n";
				next;
			}else {
				my $content = "";
				my $lane_process_dir = "$processdir/$lane";
				`mkdir -p $lane_process_dir`;
				my ($fq1, $fq2, $insert) = ($sample{$s}{$lib}{$lane}->[0],$sample{$s}{$lib}{$lane}->[1],$sample{$s}{$lib}{$lane}->[2]);
				my $dir = dirname($fq1);
				my ($adapter1, $adapter2, $fqcheck1, $fqcheck2);

				opendir DIR,$dir or die "Cannot find the raw dir $dir: $!";
				while (my $name=readdir DIR)
				{
					$adapter1 = "$dir/$name" if($name =~ /1.+adapter.+list/);
					$adapter2 = "$dir/$name" if($name =~ /2.+adapter.+list/);
					$fqcheck1 = "$dir/$name" if($name =~ /1.fqcheck/);
					$fqcheck2 = "$dir/$name" if($name =~ /2.fqcheck/);
				}
				close DIR;

				print OUT "$s\t$lib\t$lane\t$result_dir/$s/clean_data/$lane\_1.fq.gz,$result_dir/$lane/clean_data/$lane\_2.fq.gz\t$report{$lane}\t$insert\t$shell:$memory\n";
				`echo "$s\t$lane" >>$shlog`;
				my $lbase = $sample{$s}{$lib}{$lane}->[3];
				my $lcutb = $cut * $lbase / $sbase{$s};
				my $lcutr = int($lcutb/(2*$report{$lane}*1024*1024))+1;
				
				if(defined $propf)
				{
					my ($depfile, $exdep) = ($hash{$s}->[0], $hash{$s}->[1]);
					$content .= "t_depth=\`grep Average_sequencing_depth: $depfile | awk '{print \$2}'\`\n";
					$content .= "t_depth2=\$(echo \"\$t_depth/1\"|bc)\n";
					$content .= "cutReads=\$[$lbase*$exdep/(\$t_depth2*2*$report{$lane}*1024*1024)+1]\n";
				}
                                
				$filterSamp{$s} = 1;
				$lanenum++;
				if((defined $seqType) && ($seqType == 0))
				{
					$content .= "tile=`perl $Bin/bin/findNtile.pl -fq1 $fq1 -fq2 $fq2` && \\\n";
					$content .= "echo \$tile && \\\n";
					$content .= "$Bin/bin/SOAPnuke1.5.2 filter ";
					$content .= $info{"soapnuke_para"};
					$content .= " -l 5 -G";
					$content .= " -1 $fq1 -2 $fq2 -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -M 2 \$tile -o $lane_process_dir -C $lane\_1.fq -D $lane\_2.fq";
					$content .= " -c $lcutr" if($cut > 0);
					$content .= " -c \$cutReads" if(defined $propf);
					$content .= " && \\\n$Bin/bin/fqcheck33 -r $lane_process_dir/$lane\_1.fq.gz -c $lane_process_dir/$lane\_1.fqcheck && \\\n";
					$content .= "$Bin/bin/fqcheck33 -r $lane_process_dir/$lane\_2.fq.gz -c $lane_process_dir/$lane\_2.fqcheck && \\\n";
				}
				else
				{
					$content .= "tile=`perl $Bin/bin/findNtile.pl -seqType 1 -fq1 $fq1 -fq2 $fq2` && \\\n";
					$content .= "echo \$tile && \\\n";
					$content .= "$Bin/bin/SOAPnuke1.5.2 filter ";
					$content .= $info{"soapnuke_para"};
					$content .= " -l 10 -Q 2 -5 1 -7 1 -G";
					$content .= " -1 $fq1 -2 $fq2 -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -M 2 \$tile -o $lane_process_dir -C $lane\_1.fq -D $lane\_2.fq";
					$content .= " -c $lcutr" if($cut > 0);
					$content .= " -c \$cutReads" if(defined $propf);
					$content .= " && \\\n$Bin/bin/fqcheck33 -r $lane_process_dir/$lane\_1.fq.gz -c $lane_process_dir/$lane\_1.fqcheck && \\\n";
					$content .= "$Bin/bin/fqcheck33 -r $lane_process_dir/$lane\_2.fq.gz -c $lane_process_dir/$lane\_2.fqcheck && \\\n";
				}
				$content .= "perl $Bin/bin/fqcheck_distribute.pl $lane_process_dir/$lane\_1.fqcheck $lane_process_dir/$lane\_2.fqcheck -o $lane_process_dir/$lane. && \\\n";
				$content .= "mv $lane_process_dir/$lane\_1.fq.gz $lane_process_dir/$lane\_2.fq.gz $lane_process_dir/*.png $result_dir/$lane/clean_data/";
				print SHE "#!/bin/sh\n";
				print SHE $content;
				close SHE;	
			}
		}
	}
}
close OUT;

	my $samplenum = keys %filterSamp;
	`echo "SOAPnuke $samplenum samples, $lanenum lanes\n" >>$shlog`;

	   #--- clean data statistics ---#
	foreach my $s (keys %sample) {
		my $memory = "${memory}G";
		next if(!exists $unfinished{$s});
		my $content = "";
		open LIST, ">$processdir/$s/$s.stat$suffix.list" or die $!;
		open SHET,">>$sh_dir/$s/filter_$s.sh" or die $!;
		foreach my $lib (keys %{$sample{$s}}) {
			foreach my $lane (keys %{$sample{$s}{$lib}}) {
				print LIST "$processdir/$lane/$lane.stat\n";
				$content .= "perl $Bin/bin/soapnuke_stat.pl $processdir/$lane/Basic_Statistics_of_Sequencing_Quality.txt $processdir/$lane/Statistics_of_Filtered_Reads.txt > $processdir/$lane/$lane.stat && \\\n";
			}	
		}
		$content .= "perl $Bin/bin/catRS.pl $processdir/$s/$s.stat$suffix.list $result_dir/$s/clean_data/$s.xls $s";
		print SHET "\n";
		print SHET $content;
		close LIST;
		close SHET;
	}
}
#else{
 #      die "ERROR:There is no filter tool in config file.\n";
#}


### align ###
#update the fastq filepath#

#foreach my $s2 (keys %sample){
#	foreach my $lib(keys %{$sample{$s2}}){
#		foreach my $lane(keys %{$sample{$s2}{$lib}}){
#			$sample{$s2}{$lib}{$lane}->[0] = "$result_dir/$lane/clean_data/$lane\_1.fq.gz";
#			$sample{$s2}{$lib}{$lane}->[1] = "$result_dir/$lane/clean_data/$lane\_2.fq.gz";
#		}
#	}
#}

if (exists $info{'speedseq'}) {
#	open LOG, ">$shlog" or die $!;
	if ($speedseq_align == 1) {
		foreach my $lane (keys %sample){
			foreach my $lib(keys %{$sample{$lane}}){
				foreach my $s(keys %{$sample{$lane}{$lib}}){
					`mkdir -p $sh_dir/$lane`;
					`mkdir -p $result_dir/$lane`;
					my $mkd .=$result_dir;
					$mkd .= "/";
					$mkd .= $s;
					`mkdir -p $mkd` unless (-d $mkd);
	                                my $read_header = $header{$lane};
				#	$read_header = '"'."$read_header".'"';
					print "[INFO] Header :\n[INFO] $read_header\n";
				#	$cmd = "$info{'speedseq'} align -R $read_header -M $memory -t $threads -o $outdir/$s/$s $info{'reference'} $sample{$lane}{$lib}{$s}->[0] $sample{$lane}{$lib}{$s}->[1]";
					$cmd .= "$info{'speedseq'} align ";
					$cmd .= "-R ";
					$cmd .= "$read_header ";
					$cmd .= "-M $memory -t $threads -o $result_dir/$s/$s $info{'reference'} $sample{$lane}{$lib}{$s}->[0] $sample{$lane}{$lib}{$s}->[1] \n";
					
                              		`echo "speedseq align $s command :\n$cmd\n" >>$shlog`;
                                	`echo ".................\n" >>$shlog`;
					open ALN,">$sh_dir/$lane/$lane.speedseq_align.sh" or die $!;
					print ALN "#\!/bin/sh\necho =======start at -----\`date\`-----=======\n$cmd\necho ========end at -----\`date\`------===========\n";
					close ALN;
				}
			}
		}
	}
	#`rm -rf $outdir/$s.*`;
}else{
	die "ERROR: please pick align software or there is no messages for align tool in config file.\n ";
}


### calling variants : snp SNVs###
#my $alnresult_
if (((keys %sample) == 1) && ($speedseq_var == 1) && ($somatic == 0)){
	foreach my $s (keys %sample){
	foreach my $lib (keys %{$sample{$s}}){
	foreach my $lane (keys %{$sample{$s}{$lib}}){
		my $allbam = "$result_dir/$lane/$lane.bam";
		if (exists $info{'speedseq'}){
        		if($speedseq_var == 1){
                		$cmd = "$info{'speedseq'} var -o $result_dir/$lane -t 8 -k -q 1 $info{'reference'} $allbam";
                		`echo "speedseq call $lane SNVs or snp command :\n$cmd\n" >>$shlog`;
                		`echo ".................\n">>$shlog`;
				open SVA,">$sh_dir/$lane/$lane.speedseq_var.sh" or die $!;
				print SVA "#\!/bin/sh\necho =============start at -----\`date\`-----============\n$cmd\necho ===========end at -----\`date\`------===========\n";
				close SVA;

        		}
		}else{
        		die "ERROR: please pick var software or there is no messages for var tool in config file.\n "}
	}
	}
	}
}elsif((keys %sample >= 1) && ($speedseq_var == 1) && ($somatic == 0)){
	my $allbam ="";
	foreach my $lane(keys %sample){
		foreach my $lib(keys %{$sample{$lane}}){
			foreach my $s(keys %{$sample{$lane}{$lib}}){
				$allbam .= "$result_dir/$s/$s.bam ";
				if (exists $info{'speedseq'}){
                        		if($speedseq_var == 1){
                                		$cmd = "$info{'speedseq'} var -o $result_dir/$s -t 8 -k -q 1 $info{'reference'} $allbam";
                                		`echo "speedseq call $s SNVs or snp command :\n$cmd\n" >>$shlog`;
                                		`echo ".................\n>>$shlog"`;
						open SSV,">$sh_dir/$lane/$lane.speedseq_var.sh" or die $!;
						print SSV "#\!/bin/sh\necho ===========start at -----\`date\`-----===========\n$cmd\necho ===========end at -----\`date\`------============\n";
						close SSV;

                        		}
                		}else{
                        		die "ERROR: please pick var software or there is no messages for var tool in config file.\n "}


			}
		}
	}
}
### call SV ###
if (((keys %sample) == 1) && ($use_sv == 1) && ($somatic == 0)){
	my ($allbam,$disbam,$splbam);
        foreach my $s(keys %sample){
        foreach my $lib(keys %{$sample{$s}}){
        foreach my $lane(keys %{$sample{$s}{$lib}}){
		$allbam .= "$result_dir/$lane/$lane.bam";
		$disbam .= "$result_dir/$lane/$lane.discordants.bam";
		$splbam .= "$result_dir/$lane/$lane.splitters.bam";
                if (exists $info{'speedseq'}){
                        if($use_sv == 1){
				$cmd ="";
                                $cmd .= "$info{'speedseq'} sv -g -d -t 8 -m 4 -k -o $result_dir/$lane -x $info{'exclude_bed'} -R $info{'reference'} -B $allbam -S $splbam -D $disbam";
                                `echo "speedseq call $lane SV command :\n$cmd\n" >>$shlog`;
                                `echo ".................\n" >>$shlog`;
				open SAV,">$sh_dir/$lane/$lane.speedseq_sv.sh" or die $!;
                                print SAV "#\!/bin/sh\necho ===========start at -----\`date\`-----===========\n$cmd\necho ===========end at -----\`date\`------============\n";
                                close SAV;
                        }
                }else{
                        die "ERROR: please pick SV software or there is no messages for SV tool in config file.\n "}
        }
        }
        }
}elsif((keys %sample >= 1) && ($use_sv == 1) && ($somatic == 0)){
        my ($allbam,$disbam,$splbam);
	my (@allbam,@disbam,@splbam);
	foreach my $lane(keys %sample){
                foreach my $lib(keys %{$sample{$lane}}){
                        foreach my $s(keys %{$sample{$lane}{$lib}}){
                                my @bamfile = `ls $result_dir/$s |grep "bam\$"`;
                                my %bam;
                                foreach (@bamfile) {
                                        if (/splitters.bam/){
                                                chomp;
                                                # sample1.splitter.bam
                                                $bam{$s}->[0]= "$result_dir/$s/$_";
                                                push @allbam,$bam{$s}->[0];
                                        }
                                        elsif(/discordants.bam/){
                                                chomp;
                                                # sample1.discordant.bam
                                                $bam{$s}->[1] = "$result_dir/$s/$_";
						push @disbam,$bam{$s}->[1];
                                        }
                                        else{
                                                chomp;
                                                # sample1.bam
                                                $bam{$s}->[2] = "$result_dir/$s/$_";
						push @splbam,$bam{$s}->[2];
                                        }
                                }
				@allbam = sort @allbam;
				@disbam = sort @disbam;
				@splbam = sort @splbam;
				$allbam = join ',',@allbam;
				$disbam = join ',',@disbam;
				$splbam = join ',',@splbam;
                                if (exists $info{'speedseq'}){
                                        if($use_sv == 1){
                                                $cmd = "$info{'speedseq'} sv -o $result_dir/$s -x $info{'exclude_bed'} -R $info{'reference'} -B $allbam -S $splbam -D $disbam";
                                                `echo "speedseq call $s SV command :\n$cmd\n" >>$shlog`;
						 open SCSV,">$sh_dir/$lane/$lane.speedseq_sv.sh" or die $!;
                                		print SCSV "#\!/bin/sh\necho ===========start at -----\`date\`-----===========\n$cmd\necho ===========end at -----\`date\`------============\n";
                                		close SCSV;
                                        }
                                }else{
                                        die "ERROR: please pick SV software or there is no messages for SV tool in config file.\n "}


                        }
                }
        }
}

### call SNVs and SVs if it's tumor/normal pairs ###
my %sampair;
my (@sams,@tum,@norm);
if (($somatic >= 1) && (exists $info{'speedseq'}) && ($speedseq_align == 1)){
	@sams = keys %sample;
	@sams = sort @sams;
	my $i = 1;
	my @tmpsams = @sams;
	foreach (@tmpsams){
		if($i % 2 == 0){
			@tum = shift @sams;}
		else{
			@norm = shift @sams;}
		$i++;
	}
	my @tmptum = @tum;
	my @tmpnorm = @norm;
	#jump to next if $flag % 2 = 0
	my $flag = 1;
#	print "@tmptum\n@tmpnorm\n";
	foreach my $s (keys %sample){
	foreach my $lib (keys %{$sample{$s}}){
	foreach my $lane (keys %{$sample{$s}{$lib}}){
		next if($flag % 2 == 0);
		my $tumor = shift @tmptum;
		my $normal = shift @tmpnorm;
		my $tum_splitter_bam = "$result_dir/$tumor/$tumor.splitters.bam";
		my $tum_discordant_bam = "$result_dir/$tumor/$tumor.discordants.bam";
		my $tum_allbam = "$result_dir/$tumor/$tumor.bam";
		my $norm_splitter_bam = "$result_dir/$normal/$normal.splitters.bam";
                my $norm_discordant_bam = "$result_dir/$normal/$normal.discordants.bam";
                my $norm_allbam = "$result_dir/$normal/$normal.bam";
		my @tumbam = `ls $result_dir/$tumor |grep "bam\$"`;
                my %tumbam;
                foreach (@tumbam) {
                        if (/splitters.bam/){
                                chomp;
                                $tumbam{'splitter'} = "$result_dir/$tumor/$_";
                        }
                        elsif(/discordants.bam/){
                                chomp;
                                $tumbam{'discordant'} = "$result_dir/$tumor/$_";
                        }
                        else{
                                chomp;
                                $tumbam{'allbam'} = "$result_dir/$tumor/$_";
                        }
                }
		my @normbam = `ls $result_dir/$normal |grep "bam\$"`;
                my %normbam;
                foreach (@normbam) {
                        if (/splitters.bam/){
                                chomp;
                                $normbam{'splitter'} = "$result_dir/$normal/$_";
                        }
                        elsif(/discordants.bam/){
                                chomp;
                                $normbam{'discordant'} = "$result_dir/$normal/$_";
                        }
                        else{
                                chomp;
                                $normbam{'allbam'} = "$result_dir/$normal/$_";
                        }
                }
		##call snv or snp##
		my $cmd ="";
		`mkdir -p $result_dir/$lib`;
		$cmd = "$info{'speedseq'} somatic -o $result_dir/$lib -w $info{'include_bed'} -F 0.05 -q 1 $info{'reference'} $tum_allbam $norm_allbam && \\\n";
		`echo "speedseq call $lib SNVs or snp command :\n$cmd\n" >>$shlog`;
		`echo ".................\n" >>$shlog`;
		open TAV,">$sh_dir/$tumor/$tumor.speedseq_somatic.sh" or die $!;
                print TAV "#\!/bin/sh\necho ===========start at -----\`date\`-----===========\n$cmd\necho ===========end at -----\`date\`------============\n";
                close TAV;
		## call SV ##
		$cmd = "$info{'speedseq'} sv -o $result_dir/$lib -x $info{'exclude_bed'} -R $info{'reference'} -B $tum_allbam,$norm_allbam -S $tum_splitter_bam,$norm_splitter_bam -D $tum_discordant_bam,$norm_discordant_bam && \\\n";
		`echo "speedseq call $lib SV command :\n$cmd\n" >>$shlog`;
                `echo ".................\n"`;
                open TSV,">$sh_dir/$tumor/$tumor.speedseq_sv.sh" or die $!;
                print TSV "#\!/bin/sh\necho ===========start at -----\`date\`-----===========\n$cmd\necho ===========end at -----\`date\`------============\n";
                close TSV;
		$flag++;
	}
	}
	}

}

### call indels and snp###
`echo "GATKrealn&recal analysis:\n" >>$shlog`;
my $method = 0;
if(($speedseq_align == 1) && ($speedseq_var == 1)){
	$method = 1;
}elsif(($speedseq_align == 1) && ($somatic == 1)){
	$method = 2;
}elsif(($speedseq_var == 1) && ($somatic == 1)){
	print "Error : parameter -A and -z can't be simultaneously occured!\n";
}else{
	$method = 0;
	print "Error : lack of parameter -A and -V/-z ,use -A -V or -A -z.";
	die;
}
`mkdir -p $processdir/java_tmp`;
if($method == 1){
	# calculate the number of project .
	my $number = 0;
	foreach my $s (keys %sample){
        foreach my $lib (keys %{$sample{$s}}){
        foreach my $lane (keys %{$sample{$s}{$lib}}){
		`mkdir -p $processdir/$lane`;
		my @bamfile = `ls $result_dir/$lane |grep "bam\$"`;
                my %bam;
		$bam{'allbam'} = "$result_dir/$lane/$lane.bam" unless (@bamfile >= 1);
		my $cmd = "";
		# GATK RealignerTargetCreator
		$cmd .= "$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T RealignerTargetCreator -R $info{'Genome'} -I $bam{'allbam'}";
		$cmd .= " -known $info{'millsb37'}";
		$cmd .= " -known $info{'phaseb37'}";
		$cmd .= " -o $result_dir/$lane/$lane.intervals && \\\n";
		# GATK IndelRealigner
		$cmd .= "$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T IndelRealigner -R $info{'Genome'} -I $bam{'allbam'}";
		$cmd .= " -known $info{'millsb37'}";
                $cmd .= " -known $info{'phaseb37'}";
		$cmd .= " --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4";
		$cmd .= " --filter_bases_not_stored";
		$cmd .= " -targetIntervals $result_dir/$lane/$lane.intervals -o $result_dir/$lane/$lane.realign.bam && \\\n";
		# GATK BaseRecalibrator
		$cmd .= "$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar -jar $info{'GATK'} -nct 5 -T BaseRecalibrator -R $info{'Genome'} -I $result_dir/$lane/$lane.realign.bam";
		$cmd .= " -knownSites $info{'dbsnpVCFb37'}";
		$cmd .= " -knownSites $info{'millsb37'}";
		$cmd .= " -knownSites $info{'phaseb37'}";
		$cmd .= " -o $result_dir/$lane/$lane.realign.recal.table && \\\n";
		# plots the recalibration result for comparing the before and after
		$cmd .= "$info{'java'} -jar $info{'GATK'} -T BaseRecalibrator -R $info{'Genome'} -I $result_dir/$lane/$lane.realign.bam";
		$cmd .= " -knownSites $info{'dbsnpVCFb37'}";
		$cmd .= " -knownSites $info{'phaseb37'}";
		$cmd .= " -BQSR $result_dir/$lane/$lane.realign.recal.table";
		$cmd .= " -o $result_dir/$lane/$lane.post_recal_data.table && \\\n";
		$cmd .= "$info{'java'} -jar $info{'GATK'} -T AnalyzeCovariates -R $info{'Genome'}";
		$cmd .= " -before $result_dir/$lane/$lane.realign.recal.table -after $result_dir/$lane/$lane.post_recal_data.table";
		$cmd .= " -plots $result_dir/$lane/$lane.recalibration_plots.pdf && \\\n";
		# GATK PrintReads
		$cmd .= "$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -nct 5 -T PrintReads -R $info{'Genome'} -I $result_dir/$lane/$lane.realign.bam -BQSR $result_dir/$lane/$lane.realign.recal.table -o $result_dir/$lane/$lane.realign.recal.bam";
		`echo "GATK realignRecal $lane command : \n\n" >>$shlog`;
		`mkdir -p $sh_dir/$lane`;
		open RIBP,">$sh_dir/$lane/$lane.GATK_RealignerTargetCreator_Indelrealigner_BaseRecalibator_PrintReads+plots.sh" or die $!;
                print RIBP "#\!/bin/sh\necho ===========start at -----\`date\`-----===========\n$cmd\necho ===========end at -----\`date\`------============\n";
                close RIBP;
		# Call GVCF (GATK)
		my @chrom =qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M/ unless(exists $info{'fai_index'});
		my $cat_cmd ="";
		$cat_cmd ="$info{'java'} -jar $info{'GATK'} -T CombineGVCFs -R $info{'Genome'}";
		foreach my $chr(@chrom){
			$chr ="chr$chr" unless ($info{'GenomeVersion'} eq "b37");
			my $cmd ="";
			$cmd .= "$info{'java'} -Xmx5g -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T HaplotypeCaller -R $info{'Genome'} -I $result_dir/$lane/$lane.realign.recal.bam -L $chr --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $processdir/$lane/$lane.$chr.g.vcf && \\\n";
			$cmd .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$result_dir/java_tmp -jar $info{'GATK'} -T GenotypeGVCFs -R $info{'Genome'} --variant $processdir/$lane/$lane.$chr.g.vcf -stand_call_conf 30 -stand_emit_conf 10 -allSites -o $processdir/$lane/$lane.$chr.vcf";
			`mkdir -p $sh_dir/$lane` unless(-d "$sh_dir/$lane");
			open SHS,">$sh_dir/$lane/callGVCF.$chr.sh" or die $!;
			print SHS "\#!/bin/sh\necho ========Start at:-----\`date\`-----========\n";
			print SHS $cmd;
			print SHS "\necho =========End at:-----\`date\` -----========\n";
			close SHS;
			$cat_cmd .= " --variant $processdir/$lane/$lane.$chr.vcf";
		}
		$cat_cmd .= " -o $processdir/$lane/$lane.vcf --assumeSorted";
		open CA,">$sh_dir/$lane/cat_vcf.sh" or die $!;
		print CA "\#!/bin/sh\necho =======Start at:-----\`date\`-----========\n";
		print CA $cat_cmd;
		print CA "\necho =========End at:-----\`date\` -----=======\n";
		close CA;
		$number++;
	}
	}
	}
	if($number > 1){
		#reGenotype and recat the vcf
		foreach my $s (keys %sample){
        	foreach my $lib (keys %{$sample{$s}}){
        	foreach my $lane (keys %{$sample{$s}{$lib}}){
			my $comb_outdir = "$processdir/combine/";
                	my $comb_shelldir = "$sh_dir/combine/";
                	`mkdir $comb_outdir` unless(-d $comb_outdir);
                	`mkdir $comb_shelldir` unless(-d $comb_shelldir);
			my $cat_cmd = "";
			my $cmd = "";
			$cat_cmd .= "$info{'java'} -jar $info{'GATK'} -T CombineGVCFs -R $info{'Genome'}";
			my @chrom =qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M/ unless(exists $info{'fai_index'});
			foreach my $chr(@chrom){
				$chr ="chr$chr" unless ($info{'GenomeVersion'} eq "b37");
				$cmd .= "$info{'java'} -Xmx4G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T GenotypeGVCFs -R $info{'Genome'}";
				foreach my $samp(%sample){
					$cmd .= " --variant $processdir/$samp/$samp.$chr.g.vcf";
				}
				$cmd .= "-o $comb_outdir/combine.$chr.vcf";
				open CO,">$comb_shelldir/callGVCF_GATK.$chr.sh" or die $!;
				print CO "\#/bin/sh\necho =======Start at :-----\`date\` -----=======\n";
				print CO $cmd;
				print CO "\necho ========End at:-----\`date\`-----========\n";
				close CO;
				$cat_cmd .= "--variant $comb_outdir/combine.$chr.vcf";
			}
			$cat_cmd .= " -o $comb_outdir/combine.vcf --assumeSorted";
			open CAC,">$comb_shelldir/cat_vcf.sh" or die $!;
			print CAC "\#/bin/sh\necho ========Start at :-----\`date\` -----=========\n";
			print CAC $cat_cmd;
			print CAC "\necho =========End at:-----\`date\`-----==========\n";
			close CAC;
				
		}	
		}
		}
	}	
	foreach my $s (keys %sample){
        foreach my $lib (keys %{$sample{$s}}){
        foreach my $lane (keys %{$sample{$s}{$lib}}){
		my $cmd ="";
		#Get perl5lib and rpath
		my $perl5lib = `which perl`;
		$perl5lib = dirname(dirname $perl5lib);
		$perl5lib = "$perl5lib/lib";
		my $rpath = `which R`;
		print "Warnning : There is no R language(>3.1.2) exists,would cause none plots out by GATK!!!\n" if($rpath eq "");
		$cmd .= "\#!bin/sh\n";
		$cmd .= "echo =======Start at : -----\`date\`-----=========\n";
		$cmd .= "export PERL5LIB=\"$perl5lib:\$PERL5LIB\"\n";
		$cmd .= "export RPATH=\"$rpath\"\n";
		# GATK SelectVariants
		$cmd .= "$info{'java'} -Xmx10G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T SelectVariants -R $info{'Genome'} -V $processdir/$lane/$lane.vcf -selectType SNP --excludeNonVariants -o $processdir/$lane/$lane.raw.snp.vcf.gz && \\\n";
		$cmd .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T VariantRecalibrator -R $info{'Genome'} -input $processdir/$lane/$lane.raw.snp.vcf.gz \\\n";
		$cmd .= "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $info{'hapmapb37'} \\\n";
        	$cmd .= "-resource:omni,known=false,training=true,truth=true,prior=12.0 $info{'omnib37'} \\\n";
        	$cmd .= "-resource:1000G,known=false,training=true,truth=false,prior=10.0 $info{'1000Gb37'} \\\n";
        	$cmd .= "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $info{'dbsnpVCFb37'} \\\n";
        	$cmd .= "-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\\n";
        	$cmd .= "-recalFile $processdir/$lane/$lane.recalibrate_SNP.recal -tranchesFile $processdir/$lane/$lane.recalibrate_SNP.tranches -rscriptFile $processdir/$lane/$lane.recalibrate_SNP_plots.R && \\\n";
		# GATK ApplyRecalibration
		$cmd .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T ApplyRecalibration -R $info{'Genome'} -input $processdir/$lane/$lane.raw.snp.vcf.gz -mode SNP --ts_filter_level 99.0 -recalFile $processdir/$lane/$lane.recalibrate_SNP.recal -tranchesFile $processdir/$lane/$lane.recalibrate_SNP.tranches -o $processdir/$lane/$lane.filtered_snp.vcf.gz && \\\n";
		# GATK SelectVariants
        	$cmd .= "$info{'java'} -Xmx10G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T SelectVariants -R $info{'Genome'} -V $processdir/$lane/$lane.filtered_snp.vcf.gz --excludeFiltered -o $processdir/$lane/$lane.exclude.filtered_snp.vcf.gz && \\\n";
		open SN,">$sh_dir/$lane/GATK_snp.sh" or die $!;
		print SN $cmd;
		print SN "echo =========End at :-----\`date\`-----=======\n";
		close SN;
		###################################
		#this place was for annotation snp#
		###################################
		#GATK call indels
		my $content ="";
		$content .= "\#!bin/sh\n";
		$content .= "echo =====Start at : -----\`date\`-----=======\n";
                $content .= "export PERL5LIB=\"$perl5lib:\$PERL5LIB\"\n";
                $content .= "export RPATH=\"$rpath\"\n";
		# GATK  SelectVariants
		$content .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T SelectVariants -R $info{'Genome'} -V $processdir/$lane/$lane.vcf -selectType INDEL --excludeNonVariants -o $processdir/$lane/$lane.raw.indel.vcf.gz && \\\n";
		# GATK VariantRecalibrator
	        $content .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T VariantRecalibrator -R $info{'Genome'} -input $processdir/$lane/$lane.raw.indel.vcf.gz \\\n";
        	$content .= "-resource:mills,known=true,training=true,truth=true,prior=12.0 $info{'millsb37'} \\\n";
        	$content .= "-an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL \\\n";
        	$content .= "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 \\\n";
        	$content .= "-recalFile $processdir/$lane/$lane.recalibrate_INDEL.recal -tranchesFile $processdir/$lane/$lane.recalibrate_INDEL.tranches -rscriptFile $processdir/$lane/$lane.recalibrate_INDEL_plots.R && \\\n";
		# GATK ApplyRecalibration
        	$content .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T ApplyRecalibration -R $info{'Genome'} -input $processdir/$lane/$lane.raw.indel.vcf.gz -mode INDEL --ts_filter_level 99.0 -recalFile $processdir/$lane/$lane.recalibrate_INDEL.recal -tranchesFile $processdir/$lane/$lane.recalibrate_INDEL.tranches -o $processdir/$lane/$lane.filtered_indel.vcf.gz && \\\n";
		# GATK SelectVariants
        	$content .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T SelectVariants -R $info{'Genome'} -V $processdir/$lane/$lane.filtered_indel.vcf.gz --excludeFiltered -o $processdir/$lane/$lane.exclude.filtered_indel.vcf.gz && \\\n";
		open IND,">$sh_dir/$lane/GATK_indel.sh" or die $!;
                print IND $content;
                print IND "echo =========End at :-----\`date\`-----=======\n";
                close IND;
		#####################################
                #this place was for annotation indel#
                #####################################
		#$processdir/combine/combine.vcf
	}
	}
	}
	if($number > 1){
                my $cmd .= "";
                #Get perl5lib and rpath
                my $perl5lib = `which perl`;
                $perl5lib = dirname(dirname $perl5lib);
                $perl5lib = "$perl5lib/lib";
                my $rpath = `which R`;
                $cmd .= "\#!bin/sh\n";
                $cmd .= "echo ========Start at : -----\`date\`-----========\n";
                $cmd .= "export PERL5LIB=\"$perl5lib:\$PERL5LIB\"\n";
                $cmd .= "export RPATH=\"$rpath\"\n";
                # GATK SelectVariants
                $cmd .= "$info{'java'} -Xmx10G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T SelectVariants -R $info{'Genome'} -V $processdir/combine/combine.vcf -selectType SNP --excludeNonVariants -o $processdir/combine/combine.raw.snp.vcf.gz && \\\n";
		$cmd .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T VariantRecalibrator -R $info{'Genome'} -input $processdir/combine/combine.raw.snp.vcf.gz \\\n";
                $cmd .= "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $info{'hapmapb37'} \\\n";
                $cmd .= "-resource:omni,known=false,training=true,truth=true,prior=12.0 $info{'omnib37'} \\\n";
                $cmd .= "-resource:1000G,known=false,training=true,truth=false,prior=10.0 $info{'1000Gb37'} \\\n";
                $cmd .= "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $info{'dbsnpVCFb37'} \\\n";
                $cmd .= "-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\\n";
                $cmd .= "-recalFile $processdir/combine/combine.recalibrate_SNP.recal -tranchesFile $processdir/combine/combine.recalibrate_SNP.tranches -rscriptFile $processdir/combine/combine.recalibrate_SNP_plots.R && \\\n";
                # GATK ApplyRecalibration
                $cmd .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T ApplyRecalibration -R $info{'Genome'} -input $processdir/combine/combine.raw.snp.vcf.gz -mode SNP --ts_filter_level 99.0 -recalFile $processdir/combine/combine.recalibrate_SNP.recal -tranchesFile $processdir/combine/combine.recalibrate_SNP.tranches -o $processdir/combine/combine.filtered_snp.vcf.gz && \\\n";
                # GATK SelectVariants
                $cmd .= "$info{'java'} -Xmx10G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T SelectVariants -R $info{'Genome'} -V $processdir/combine/combine.filtered_snp.vcf.gz --excludeFiltered -o $processdir/combine/combine.exclude.filtered_snp.vcf.gz && \\\n";
		open SN,">$sh_dir/combine/GATK_snp.sh" or die $!;
                print SN $cmd;
                print SN "echo =========End at :-----\`date\`-----========\n";
                close SN;
                ###################################
                #this place was for annotation snp#
                ###################################
                #GATK call indels
                my $content ="";
                $content .= "\#!bin/sh\n";
                $content .= "echo =======Start at : -----\`date\`-----======\n";
                $content .= "export PERL5LIB=\"$perl5lib:\$PERL5LIB\"\n";
                $content .= "export RPATH=\"$rpath\"\n";
                # GATK  SelectVariants
                $content .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T SelectVariants -R $info{'Genome'} -V $processdir/combine/combine.vcf -selectType INDEL --excludeNonVariants -o $processdir/combine/combine.raw.indel.vcf.gz && \\\n";
                # GATK VariantRecalibrator
                $content .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T VariantRecalibrator -R $info{'Genome'} -input $processdir/combine/combine.raw.indel.vcf.gz \\\n";
                $content .= "-resource:mills,known=true,training=true,truth=true,prior=12.0 $info{'millsb37'} \\\n";
                $content .= "-an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL \\\n";
                $content .= "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 \\\n";
                $content .= "-recalFile $processdir/combine/combine.recalibrate_INDEL.recal -tranchesFile $processdir/combine/combine.recalibrate_INDEL.tranches -rscriptFile $processdir/combine/combine.recalibrate_INDEL_plots.R && \\\n";
                # GATK ApplyRecalibration
                $content .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T ApplyRecalibration -R $info{'Genome'} -input $processdir/combine/combine.raw.indel.vcf.gz -mode INDEL --ts_filter_level 99.0 -recalFile $processdir/combine/combine.recalibrate_INDEL.recal -tranchesFile $processdir/combine/combine.recalibrate_INDEL.tranches -o $processdir/combine/combine.filtered_indel.vcf.gz && \\\n";
                # GATK SelectVariants
                $content .= "$info{'java'} -Xmx5G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T SelectVariants -R $info{'Genome'} -V $processdir/combine/combine.filtered_indel.vcf.gz --excludeFiltered -o $processdir/combine/combine.exclude.filtered_indel.vcf.gz && \\\n";
		open IND,">$sh_dir/combine/GATK_indel.sh" or die $!;
                print IND $content;
                print IND "echo ======End at :-----\`date\`-----==========\n";
                close IND;
                #####################################
                #this place was for annotation indel#
                #####################################
	
	}

}elsif($method == 2){
	if(exists $info{'muTect'}){
		foreach my $s (keys %sample){
                foreach my $lib (keys %{$sample{$s}}){
                foreach my $lane (keys %{$sample{$s}{$lib}}){
			my $cmd = "";
			$cmd .="#\!/bin/sh\n";
                        $cmd .= "echo =======Start at : -----\`date\`-----======\n";
                        # GATK RealignerTargetCreator
                        $cmd .= "$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T RealignerTargetCreator -R $info{'Genome'} -I $result_dir/$lane/$lane.bam";
                        $cmd .= " -known $info{'millsb37'}";
                        $cmd .= " -known $info{'phaseb37'}";
                        $cmd .= " -o $result_dir/$lane/$lane.intervals && \\\n";
                        # GATK IndelRealigner
                        $cmd .= "$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -T IndelRealigner -R $info{'Genome'} -I $result_dir/$lane/$lane.bam";
                        $cmd .= " -known $info{'millsb37'}";
                        $cmd .= " -known $info{'phaseb37'}";
			$cmd .= " --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4";
			$cmd .= " --filter_bases_not_stored";
                        $cmd .= " -targetIntervals $result_dir/$lane/$lane.intervals -o $result_dir/$lane/$lane.realign.bam && \\\n";
                        # GATK BaseRecalibrator
                        $cmd .= "$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar -jar $info{'GATK'} -nct 5 -T BaseRecalibrator -R $info{'Genome'} -I $result_dir/$lane/$lane.realign.bam";
                        $cmd .= " -knownSites $info{'dbsnpVCFb37'}";
                        $cmd .= " -knownSites $info{'millsb37'}";
                        $cmd .= " -knownSites $info{'phaseb37'}";
                        $cmd .= " -o $result_dir/$lane/$lane.realign.recal.table && \\\n";
                        # plots the recalibration result for comparing the before and after
                        $cmd .= "$info{'java'} -jar $info{'GATK'} -T BaseRecalibrator -R $info{'Genome'} -I $result_dir/$lane/$lane.realign.bam";
                        $cmd .= " -knownSites $info{'dbsnpVCFb37'}";
                        $cmd .= " -knownSites $info{'phaseb37'}";
                        $cmd .= " -BQSR $result_dir/$lane/$lane.realign.recal.table";
                        $cmd .= " -o $result_dir/$lane/$lane.post_recal_data.table && \\\n";
                        $cmd .= "$info{'java'} -jar $info{'GATK'} -T AnalyzeCovariates -R $info{'Genome'}";
                        $cmd .= " -before $result_dir/$lane/$lane.realign.recal.table -after $result_dir/$lane/$lane.post_recal_data.table";
                        $cmd .= " -plots $result_dir/$lane/$lane.recalibration_plots.pdf && \\\n";
                	# GATK PrintReads
                        $cmd .= "$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -nct 5 -T PrintReads -R $info{'Genome'} -I $result_dir/$lane/$lane.realign.bam -BQSR $result_dir/$lane/$lane.realign.recal.table -o $result_dir/$lane/$lane.realign.recal.bam";
#			$cmd .="$info{'java'} -Xmx4g -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'ContEst'} -T Contamination -I $result_dir/$lane/$lane.realign.recal.bam -R $info{'reference'} -B:pop,vcf $info{'pop_vcf'} -B:genotypes,vcf ";
                        `echo "GATK realignRecal $lane command : \n\n" >>$shlog`;
			
			open SOR,">$sh_dir/$lane/$lane.RealignerTargetCreator_IndelRealigner_BaseRecalibrator_PrintReads+plots.sh" or die $!;
			print SOR $cmd;
			print SOR "\necho End at :======-----\`date\`------=======\n";
			close SOR;
		}
		}
		}
		my (@sams,@tum,@norm);
	        @sams = keys %sample;
        	@sams = sort @sams;
		my $i = 1;
		my @tmpsams = @sams;
        	foreach (@tmpsams){
                	if($i % 2 == 0){
                        	@tum = shift @sams;}
                	else{
                        	@norm = shift @sams;}
			$i++;
        	}
		my @tmptum = @tum;
        	my @tmpnorm = @norm;
		#jump to next if $flag % 2 = 0
	        my $flag = 1;
        	foreach my $s (keys %sample){
        	foreach my $lib (keys %{$sample{$s}}){
        	foreach my $lane (keys %{$sample{$s}{$lib}}){
                	next if($flag % 2 == 0);
                	my $tumor = shift @tmptum;
                	my $normal = shift @tmpnorm;
			my $tumbam = "$result_dir/$tumor/$tumor.realign.recal.bam";
			my $normbam = "$result_dir/$normal/$normal.realign.recal.bam";
			my $cmd = "";
			#Estimate Contamination by Conpair
			$cmd .="#\!/bin/sh\n";
			$cmd .="echo Start at:=====-------\`date\`--------=======\n";
			##use conpair to estimate contamination( outdated:not use )
#			chomp(my $conpairpath =`find $info{'Bin'} -name run_gatk_pileup_for_sample.py`);
#
#                       die "[Error]:please install Conpair tool.[URL]https://github.com/nygenome/Conpair/archive/master.zip\n" if($conpairpath eq "");
#                      $conpairpath =abs_path($conpairpath);
#                        my $conpair_src =dirname($conpairpath);
#                        $cmd .="export GATK_JAR=$info{'GATK'}\n";
#                        $conpairpath =dirname(dirname $conpairpath);
#                        $cmd .="export CONPAIR_DIR=$conpairpath\n";
#                        chomp(my $pythonpath =`which python`);
#                        $pythonpath =abs_path($pythonpath);
#                        $cmd .="export PYTHONPATH=$pythonpath\n";
#                        $cmd .="export PYTHONPATH=\${PYTHONPATH}:$conpairpath/modules\n";
#			$cmd .="$conpair_src/run_gatk_pileup_for_sample.py -B $tumbam -O $result_dir/$tumor/$tumor._pileup --reference $info{'reference'} && \\\n";
#			$cmd .="$conpair_src/run_gatk_pileup_for_sample.py -B $normbam -O $result_dir/$normal/$normal._pileup --reference $info{'reference'} && \\\n";

			#$conpairpath =dirname($conpairpath);
#			$cmd .="$conpairpath/scripts/verify_concordance.py -T $result_dir/$tumor/$tumor._pileup -N $result_dir/$normal/$normal._pileup --min_cov 10 --outfile $tumor.pair.pileup.verify_concordance.txt && \\\n";
#			$cmd .="$conpairpath/scripts/estimate_tumor_normal_contamination.py -T $result_dir/$tumor/$tumor._pileup -N $result_dir/$normal/$normal._pileup --outfile $tumor.pair.estimate_contamination.txt && \\\n";
			
			##GATK ContEst to estimate contamination
			$cmd .="$info{'java'} -jar $info{'GATK'} -T ContEst -R $info{'Genome'} -I:eval $tumbam -I:genotype $normbam --popFile $info{'pop_vcf'} -L $info{'interval_list'} -isr INTERSECTION -o $result_dir/$tumor/ContEst_contamination_estimate.txt && \\\n";
			# muTect call point mutations
			$cmd .="$info{'java7'} -Xmx3G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'muTect'} --analysis_type MuTect --reference_sequence $info{'Genome'} --cosmic $info{'cosmic'} --dbsnp $info{'leftAligned'} --input_file:normal $normbam --input_file:tumor $tumbam --out $result_dir/$tumor/$tumor.call_stats.txt --coverage_file $result_dir/$tumor/$tumor.coverage.wig.txt && \n" if ($mutect2 == 0);
			$cmd .="$info{'java'} -Xmx3G -Djava.io.tmpdir=$processdir/java_tmp -jar $info{'GATK'} -MuTect2  -R $info{'Genome'} --cosmic $info{'cosmic'} --dbsnp $info{'leftAligned'} -I:normal $normbam -I:tumor $tumbam -o $result_dir/$tumor/$tumor.call_stats.txt --coverage_file $result_dir/$tumor/$tumor.coverage.wig.txt && \n" if ($mutect2 > 0);
			open MUO,">$sh_dir/$tumor/$tumor.somatic_snps_indels_muTect.sh" or die $!;
			print MUO $cmd;
			print MUO "echo End at :======-------\`data\`-------========\n";
			close MUO;
			$flag++;
		}
		}
		}
		foreach my $s (keys %sample){
                foreach my $lib (keys %{$sample{$s}}){
                foreach my $lane (keys %{$sample{$s}{$lib}}){
			my $cmdpre="";
			#CREST Call SVs
			chomp(my $crestpath =`find $info{'Bin'} -name CREST.pl`);
			print "[Error]:please install CREST  tool.\n" if($crestpath eq "");
                        $crestpath = abs_path($crestpath);
                        $crestpath = dirname($crestpath);
			my $genome2bit_dir = dirname($info{'Genome'});
			my $genome2bit_name = basename($info{'Genome'});
			my $cmd_2bit .="$crestpath/blatSuite-34/faToTwoBit $info{'Genome'} $genome2bit_dir/$genome2bit_name.2bit && \\\n";
			chomp(my $serverip =`/sbin/ifconfig eth0 |grep "inet addr" |cut -f2 -d":" |cut -d " " -f1`);
                        $cmdpre .="$crestpath/blatSuite-34/gfServer start $serverip 9001 $genome2bit_dir/$genome2bit_name.2bit && \\\n";
			$cmdpre .="export PATH=\"$crestpath:\$PATH\"\n";
			##my $perl5lib = `which perl`;
                	##$perl5lib = dirname(dirname $perl5lib);
                	##$perl5lib = "$perl5lib/lib";
			$cmdpre .="export PERL5LIB=\"$crestpath/lib/lib/perl5:$crestpath/lib/lib64/perl5:$crestpath:\$PERL5LIB\"\n";
	                my @chrom =qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M/ unless(exists $info{'crest_chr'});
			foreach my $chr(@chrom){
				$chr ="chr$chr";
				`mkdir -p $result_dir/$lane/crest`;
				my $hg19fa =dirname($info{'reference'});
				my $cmd = $cmd_2bit if($chr eq 'chr1');
				$cmd .= $cmdpre;
				#$cmd .="ln -s $hg19fa/human_g1k_v37.fasta $hg19fa/human_g1k_v37.fa &&\n";
				#$cmd .="ln -s $hg19fa/human_g1k_v37.fasta.bwt $hg19fa/human_g1k_v37.fa.bwt\n";
				$cmd .="perl $crestpath/extractSClip.pl -i $result_dir/$lane/$lane.realign.recal.bam --ref_genome $info{'reference'} -o $result_dir/$lane/crest -r $chr && \\\n";
				`mkdir -p $sh_dir/$lane/crest`;
				$cmd .="\necho end at :========-------\`date\`---------=========\n";		
				open SCLG,">$sh_dir/$lane/crest/$lane.$chr.SClip.sh";
				print SCLG "#\!/bin/bash\necho ======------\`date\`------=======\n";
				print SCLG $cmd;
				close SCLG;
				$cmd = $cmdpre;
				my $prefix = "$lane.realign.recal";
				$cmd .="perl $crestpath/CREST.pl -f $result_dir/$lane/$prefix.$chr.cover -d $result_dir/$lane/$lane.realign.recal.bam -p $prefix.$chr -r $chr ";
				$cmd .="--ref_genome $info{'reference'} -t $genome2bit_dir/$genome2bit_name.2bit --cap3 $crestpath/CAP3/cap3 --blatclient $crestpath/blatSuite-34/gfClient --blat $crestpath/blatSuite-34/blat --blatserver $serverip --blatport 9001 ";
				$cmd .="-o $result_dir/$lane/crest -l $report{$lane} && \\\n";
				$cmd .="echo ========-------\`date\`-------========\n";
				open CRES,">$sh_dir/$lane/crest/$lane.$chr.CREST.sh" or die $!;
				print CRES "#\!bin/sh\necho Start at :=========--------\`date\`----------=========\n";
				print CRES $cmd;
				close CRES;				
			}
			$cmd ="";
			$cmd ="#\!bin/sh\n";
			my $prefix = "$lane.realign.recal";
			$cmd .="cat $result_dir/$lane/crest/$prefix.*.cover > $result_dir/$lane/$prefix.cover && \\\n";
			$cmd .="cat $result_dir/$lane/crest/$prefix.*.sclip.cover > $result_dir/$lane/$prefix.sclip.cover && \\\n";
			$cmd .="cat  $result_dir/$lane/crest/$prefix.*.predSV.txt > $result_dir/$lane/$prefix.predSV.txt && \\\n";
			$cmd .="cat $crestpath/crest_head $result_dir/$lane/$prefix.predSV.txt >$result_dir/$lane/$prefix.cathead.predSV.txt && \\\n";
			$cmd .="\necho end at :=======--------\`date\`--------========\n";
			open CATH,">$sh_dir/$lane/crest/$prefix.cat.cover.sclip.predSV.head.sh" or die $!;
			print CATH $cmd;
			close CATH;
		}
		}
		}		
	}	
	
}



### run shell command ###

if($method == 1){
	foreach my $s (keys %sample){
	foreach my $lib (keys %{$sample{$s}}){
	foreach my $lane (keys %{$sample{$s}{$lib}}){
		if(keys %sample == 1){
			my $cmd ="";
			$cmd .="sh $outdir/set_path.sh && \\\n";
			$cmd .="sh $sh_dir/$lane/$lane.speedseq_align.sh && \\\n";
			$cmd .="sh $sh_dir/$lane/$lane.speedseq_var.sh && \\\n";
			$cmd .="sh $sh_dir/$lane/$lane.speedseq_sv.sh && \\\n";
			$cmd .="sh $sh_dir/$lane/$lane.GATK_RealignerTargetCreator_Indelrealigner_BaseRecalibator_PrintReads+plots.sh && \\\n";
			my @chrom =qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M/ unless(exists $info{'fai_index'});
			foreach my $chr(@chrom){
				$chr = "chr$chr" unless ($info{'GenomeVersion'} eq "b37");
				$cmd .="sh $sh_dir/$lane/callGVCF.$chr.sh \\\n";
			}
			$cmd .=" && sh $sh_dir/$lane/cat_vcf.sh && ";
			$cmd .="sh $sh_dir/$lane/GATK_snp.sh && ";
			$cmd .="sh $sh_dir/$lane/GATK_indel.sh";
			
		}else{
			my $cmd ="";
			$cmd .="sh $outdir/set_path.sh && \\\n";
                        $cmd .="sh $sh_dir/$lane/$lane.speedseq_align.sh && \\\n";
                        $cmd .="sh $sh_dir/$lane/$lane.speedseq_var.sh && \\\n";
                        $cmd .="sh $sh_dir/$lane/$lane.speedseq_sv.sh && \\\n";
                        $cmd .="sh $sh_dir/$lane/$lane.GATK_RealignerTargetCreator_Indelrealigner_BaseRecalibator_PrintReads+plots.sh && \\\n";			
			my @chrom =qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M/ unless(exists $info{'fai_index'});
			foreach my $chr(@chrom){
                                $chr = "chr$chr" unless ($info{'GenomeVersion'} eq "b37");
                                $cmd .="sh $sh_dir/$lane/callGVCF.$chr.sh \\\n";
                        }
			$cmd .=" && sh $sh_dir/$lane/cat_vcf.sh && ";
                        $cmd .="sh $sh_dir/$lane/GATK_snp.sh && ";
                        $cmd .="sh $sh_dir/$lane/GATK_indel.sh";

                        foreach my $chr(@chrom){
                                $chr = "chr$chr" unless ($info{'GenomeVersion'} eq "b37");
                                $cmd .="sh $sh_dir/combine/callGVCF.$chr.sh \\\n";
                        }
                        $cmd .=" && sh $sh_dir/combine/cat_vcf.sh && ";
                        $cmd .="sh $sh_dir/combine/GATK_snp.sh && ";
                        $cmd .="sh $sh_dir/combine/GATK_indel.sh";

		}	
	}
	}
	}

}elsif($method == 2){
	my $cmd ="";
	$cmd ="sh $outdir/set_path.sh && \\\n";
	my (@sams,@tum,@norm);
        @sams = keys %sample;
        @sams = sort @sams;
        my $i = 1;
        my @tmpsams = @sams;
       	foreach (@tmpsams){
        	if($i % 2 == 0){
                        @tum = shift @sams;}
                else{
                        @norm = shift @sams;}
                        $i++;
                }
        my @tmptum = @tum;
        my @tmpnorm = @norm;
        #jump to next if $flag % 2 = 0
        my $flag = 1;
	foreach my $s (keys %sample){
        foreach my $lib (keys %{$sample{$s}}){
        foreach my $lane (keys %{$sample{$s}{$lib}}){
		next if($flag % 2 == 0);
       	        my $tumor = shift @tmptum;
                my $normal = shift @tmpnorm;
		$cmd .="sh $sh_dir/$tumor/$tumor.speedseq_align.sh  \\\n";
		$cmd .="sh $sh_dir/$normal/$normal.speedseq_align.sh  && \\\n";
		$cmd .="sh $sh_dir/$tumor/$tumor.GATK_RealignerTargetCreator_Indelrealigner_BaseRecalibator_PrintReads+plots.sh \\\n";
		$cmd .="sh $sh_dir/$normal/$normal.GATK_RealignerTargetCreator_Indelrealigner_BaseRecalibator_PrintReads+plots.sh && \\\n";
		$cmd .="sh $sh_dir/$tumor/$tumor.somatic_snps_indels_muTect.sh && \\\n";
                my @chrom =qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M/ unless(exists $info{'fai_index'});
                foreach my $chr(@chrom){
              		$chr = "chr$chr" unless ($info{'GenomeVersion'} eq "b37");
                        $cmd .="sh $sh_dir/$tumor/crest/$tumor.$chr.SClip.sh \\\n";
                }
		$cmd .=" && ";
		foreach my $chr(@chrom){
                        $chr = "chr$chr" unless ($info{'GenomeVersion'} eq "b37");
                        $cmd .="sh $sh_dir/$normal/crest/$normal.$chr.SClip.sh \\\n";
                }
		$cmd .=" && ";
		 foreach my $chr(@chrom){
                        $chr = "chr$chr" unless ($info{'GenomeVersion'} eq "b37");
                        $cmd .="sh $sh_dir/$tumor/crest/$tumor.$chr.CREST.sh \\\n";
                }
		$cmd .=" && ";
		 foreach my $chr(@chrom){
                        $chr = "chr$chr" unless ($info{'GenomeVersion'} eq "b37");
                        $cmd .="sh $sh_dir/$normal/crest/$normal.$chr.CREST.sh \\\n";
                }
		$cmd .=" && ";
		$cmd .="sh $sh_dir/$tumor/crest/$tumor.realign.recal.cat.cover.sclip.predSV.head.sh \\\n";
		$cmd .="sh $sh_dir/$normal/crest/$normal.realign.recal.cat.cover.sclip.predSV.head.sh \\\n";
		$cmd .="sh $sh_dir/$tumor/crest/$tumor.speedseq_somatic.sh && \\\n";
		$cmd .="sh $sh_dir/$tumor/crest/$tumor.speedseq_sv.sh && \\\n";


		
		$flag++;
	}
	}
	}
}else{
	die "fatal error: not enough parameters specified.\n";
}


## END SCRIPT

#                             .       .
#                            / `.   .' \
#                    .---.  <    > <    >  .---.
#                    |    \  \ - ~ ~ - /  /    |
#                     ~-..-~             ~-..-~
#                 \~~~\.'                    `./~~~/
#       .-~~^-.    \__/                        \__/
#     .'  O    \     /               /       \  \
#    (_____,    `._.'               |         }  \/~~~/
#     `----.          /       }     |        /    \__/
#           `-.      |       /      |       /      `. ,~~|
#               ~-.__|      /_ - ~ ^|      /- _      `..-'   f: f:
#                    |     /        |     /     ~-.     `-. _||_||_
#                    |_____|        |_____|         ~ - . _ _ _ _ _>
#

























