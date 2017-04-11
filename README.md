# HiseqWGS_2017 (Test version)

A WGS pipeline constructed by multiple tools. Meant to  produce shell scripts for analysis the whole genome sequencing data based on provided parameters.<br>

Whole pipeline : raw data -> filter -> align -> correct bam -> call SNPs and indels -> call SVs -> call SNVs <br>

* __Author__: Songjing <br>
* __E-mail__: songjing@genomics.cn<br>
* __Edit__: 14 Mar,2017<br>

## Overview <br>
This pipeline contains germline variants calling as well as somatic variants calling.You can just use few parameters then it can procedure all shell scripts to run each step.

## Usage <br>
#### For single/multiple sample(s):
	 perl HiseqWGS_2017.pl -c HiseqWGS_2017.conf -f rawData.list -A -V -S
#### For tumor/normal pair(s) samples:
	 perl HiseqWGS_2017.pl -c HiseqWGS_2017.conf -f rawData.list -A -z -S

##### options: <br>
        -c:*    STR     configure file,including the path of tools used.(eg:HiseqWGS_2017.conf)(required)
        -f:*    STR     sample fastq list(required)
        -R:     STR     read group header line such as "\@RG\\tID:id\\tSM:samplename\\tLB:lib"[Default: Auto]
        -m:     STR	memory to use [Default:12 (Gb)]
        -p:     Boolens first fastq file consists of interleaved paired-end sequences[Default: False]
        -t:     INT     number of threads to use [default:16]
        -o:     STR     output dir prefix  [Defalut:'./outdir_\$date']
        -tmp:   STR     TMP dir prefix     [Default:'./TMP_\$date']

        -A:*    Boolens Call align(align+markdup+sort+index) [Default:False]
        -V:*    Boolens Call var  (calling SNP,indels)       [Default:False]
        -S:*    Boolens Call sv   (calling sv,CNV)           [Default:False]

        -trio:  Boolens Necessary if your samples came from a family(eg:child,father and mother)[Defalut: False]
        -z:     Boolens Necessary if your samples are tumor/normal pairs[Default: False]

        -help|?:        show this usage page
        -more:          more default information


## Prerequisites<br>
* GNU gcc =4.7 ~ 5.0 , (gcc 4.8.4 used here)
* g++
* python 2.7
   	* pysam 0.8.4 or newer(0.9.1 used here)
	*	numpy 
	*	scipy
* automake 
* make
* cmake 
* git 
*	libncurses5-dev
*	zlib1g-dev
*	g++ 
*	gcc 
*	gfortran
*	python2.7-dev 
*	python-pip 
*	libblas-doc 
*	libblas3
*	libblas-dev
*	liblapack3 
*	liblapack-dev
*	zlib-1.2.11
*	libX11-dev 
*	libxpm-dev
*	libxft-dev
*	libxext-dev 
*	dpkg-dev 
*	binutils
*	libssl-dev 
*	libpcre3-dev
*	xlibmesa-glu-dev
*	libglew1.5-dev 
*	libftgl-dev
*	libmysqlclient-dev
*	libfftw3-dev 
*	libcfitsio3-dev 
*	graphviz-dev 
*	libavahi-compat-libdnssd-dev
*	libldap2-dev 
*	python-dev
*	libxml2-dev
*	libkrb5-dev 
*	libgsl0-dev 
*	libqt4-dev
*	root(v5.34.20)
*	samtools(1.3.1)
*	bwa_mem(0.7.12)
*	java 8
*	picard-tools(1.54)
*	Speedseq(v0.1.2)
*	Conpair
*	CREST
*	muTect
*	GATK (3.6)
*	GASVPro
*	SOAPnuke(no need if can't get)
*	fastQC
*	Oncotator
*	VEP
*	GISTIC
## installation<br>
##### You can use existed script 'install_HiseqWGS_2017.pl' to install all the environment exclude reference data.<br>
	perl < path to install_HiseqWGS_2017.pl> /install_HiseqWGS_2017.pl [-o] < OUT_PATH >
	

