# HiseqWGS_2017

A WGS pipeline constructed by multiple tools.Meant to  produce shell scripts for analysis the whole genome sequencing data based on provided parameters.
WHOLE PIPELINE:raw data -> filter -> align -> correct bam -> call SNPs and indels -> call SVs

__Author__: Songjing ,__E-mail__: songjing@genomics.cn
__Edit__: 14 Mar,2017

#Overview
This pipeline contains germline variants calling as well as somatic variants calling.You can just use few parameters then it can procedure all shell scripts to run each step.

#Prerequisites
	GNU gcc =4.7 ~ 5.0 , (gcc 4.8.4 used here)<br>
	g++
	python 2.7
	automake 
	make
	cmake 
	git 
	libncurses5-dev
	zlib1g-dev
	g++ 
	gcc 
	gfortran
	python2.7-dev 
	python-pip 
	libblas-doc 
	libblas3
	libblas-dev
	liblapack3 
	liblapack-dev
	pysam 0.8.4 or newer(0.9.1 used here)
	numpy 
	scipy
	zlib-1.2.11
	libX11-dev 
	libxpm-dev
	libxft-dev
	libxext-dev 
	dpkg-dev 
	binutils
	libssl-dev 
	libpcre3-dev
	xlibmesa-glu-dev
	libglew1.5-dev 
	libftgl-dev
	libmysqlclient-dev
	libfftw3-dev 
	libcfitsio3-dev 
	graphviz-dev 
	libavahi-compat-libdnssd-dev
	libldap2-dev 
	python-dev
	libxml2-dev
	libkrb5-dev 
	libgsl0-dev 
	libqt4-dev
	root(v5.34.20)
	samtools(1.3.1)
	bwa_mem(0.7.12)
	java 8
	picard-tools(1.54)
	Speedseq(v0.1.2)
	Conpair
	CREST
	muTect
	GATK (3.6)
	GASVPro
	SOAPnuke(no need if can't get)
	fastQC
	Oncotator
	VEP
	GISTIC
	

