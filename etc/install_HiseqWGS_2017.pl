#!/usr/bin/env perl -w
# perl script #
#
# This script was meant to build the 'Running Environment' of WGS pipeline
#
# edit : 7 Apr,2017. 
# ©songjing   E-mail：kakaluote707@sina.com
# Version 1.2

# Limitation : GNU gcc =4.7 ~ 5.0 , (gcc 4.8.4 used here)
#              python 2.7 or newer , python modules(pysam > 0.8.1 )
#	       GATK 3.6 , java 8 , speedseq(svtyper 0.0.4 used here)
# 
 		
use strict;
use File::Path;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use FindBin '$Bin';

sub usage{
	print STDOUT << "USAGE";
OPTIONS:
	-o:		path of where you want to install [./]
	-help|h:	help messages
	[tips: this scripts was for fedora/Ubuntu 14.04]
USAGE
	exit 0;
}

my ($outdir,$help);

GetOptions(
	"o:s" => \$outdir,
	"help|h+" => \$help,
);
if($help >= 1){
	&usage;
	die;
}
$outdir = abs_path($Bin);
$outdir ||= "$outdir/HiseqWGS_2017";
print "install to :$outdir\n";
`mkdir -p $outdir/opt`;
`mkdir -p $outdir/opt/bin`;
`mkdir -p $outdir/opt/prerequisites` unless (-d "$outdir/opt/prerequisites");
`mkdir -p $outdir/opt/Database` unless (-d "$outdir/opt/Database");
`mkdir -p $outdir/tmp` unless(-d "$outdir/tmp");
print "creat directory: \n$outdir\n$outdir/opt\n$outdir/opt/bin\n$outdir/opt/prerequisites\n$outdir/opt/Database\n ";
print "creat tmp dir: $outdir/tmp\n";

print " \n * Inatalling HiseqWGS_2017 running environment.\n";
my $cmd = "";
#speedseq_v0.1.2 code
$cmd .="cd $outdir/tmp && ";
$cmd .="git clone --recursive https://github.com/cc2qe/speedseq.git && ";
$cmd .="mv $outdir/tmp/speedseq $outdir/opt/bin && ";
#prerequisites
$cmd .="apt-get install automake make cmake git libncurses5-dev zlib1g-dev g++ gcc gfortran && ";
$cmd .="apt-get install python2.7-dev python-pip libblas-doc libblas3 libblas-dev liblapack3 liblapack-dev && ";
#python2.7  modules
$cmd .="pip install --upgrade pip && ";
$cmd .="pip install pysam==0.9.1 && ";
$cmd .="pip install numpy && ";
$cmd .="pip install scipy && ";
#zlib_1.2.11 install
$cmd .="wget https://fossies.org/linux/misc/zlib-1.2.11.tar.gz && ";
$cmd .="tar -zxvf zlib-1.2.11* && ";
$cmd .="cd zlib* && ./configure --prefix $outdir/opt/bin/speedseq && make && make install && ";
$cmd .="rm zlib*.gz && mv zlib* $outdir/opt/prerequisites && ";
#root prerequisites
$cmd .="apt-get install libX11-dev libxpm-dev libxft-dev libxext-dev dpkg-dev g++ gcc binutils && ";
$cmd .="apt-get install git dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev && ";
$cmd .="apt-get install gfortran libssl-dev libpcre3-dev xlibmesa-glu-dev libglew1.5-dev libftgl-dev libmysqlclient-dev libfftw3-dev libcfitsio3-dev graphviz-dev libavahi-compat-libdnssd-dev libldap2-dev python-dev libxml2-dev libkrb5-dev libgsl0-dev libqt4-dev && ";
#root_v5.34.20 install
$cmd .="wget https://root.cern.ch/download/root_v5.34.20.source.tar.gz && ";
$cmd .="tar -vxf root_v5.34.20.source.tar.gz && ";
$cmd .="rm root_*gz && mv root_* $outdir/opt/prerequisites && ";
$cmd .="echo \"making root_v5.34.20,it may takes 1 hour~3 hours,..wait.......\"";
$cmd .=" && ";
$cmd .="cd $outdir/opt/prerequisites/root && ./configure && make -j4 && ";
$cmd .="echo \"source $outdir/opt/prerequisites/root/bin/thisroot.sh\" >> ~/.bashrc && source ~/.bashrc && ";
$cmd .="apt-get install ncftp && perl -MCPAN -e shell && ";
$cmd .="\`yes\` && \`yes\` && \`yes\` && \`yes\` && \`yes\` && echo -e \"install Archive::Extract\ninstall CGI\ninstall DBI\ninstall Time::HiRes\ninstall Archive::Tar\ninstall Archive::Zip\ninstall Module::Build\ninstall File::Copy::Recursive\n\" && echo -e \"exit\" && ";
$cmd .="cd $outdir/opt/bin/speedseq && make && make cnvnator-multi && ";

$cmd .="echo -e \"export PATH=$outdir/opt/bin/speedseq/bin:\$PATH\" >> ~/.bashrc && source ~/.bashrc && ";
$cmd .="cd $outdir/tmp && wget https://github.com/hall-lab/svtyper/archive/v0.0.4.tar.gz && tar -zxvf svtyper-0.0.4* && rm svtyper-0.0.4*gz && cp svtyper*/svtyper $outdir/opt/bin/speedseq/bin && ";
#bwamem 0.7.12,  samtools 1.3.1
$cmd .="wget https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2/download && tar -jxvf samtools-1.3.1.tar.bz2 && rm samtools*bz2 && cd samtools* && ./configure && make && cd ../ && cp -r samtools* $outdir/opt/bin/samtools && echo -e \"export PATH=$outdir/opt/bin/samtools:\$PATH\" >> ~/.bashrc && source ~/.bashrc && ";
$cmd .="wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz && tar -zxvf bwa* && rm bwa*gz && cd bwa* && ./configure && make && cp -r bwa* $outdir/opt/bin/bwa_mem && echo -e \"export PATH=$outdir/opt/bin/bwa_mem:\$PATH\" >> ~/.bashrc && source ~/.bashrc && ";
#configure path
$cmd .="echo \"export ROOTSYS=$outdir/opt/prerequisites/root\" >>~/.bashrc && echo \"export PATH=\$PATH:\$ROOTSYS/bin\" >>~/.bashrc && echo \"export LD_LIBRARY_PATH=/usr/lib:/usr/lib64:/root/lib:$outdir/opt/prerequisites/zlib-1.2.11:\$LD_LIBRARY_PATH:\$ROOTSYS/lib:\$LD_LIBRARY_PATH\" >>~/.bashrc && echo \"export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ROOTSYS/lib\" >>~/.bashrc && echo \"export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$outdir/opt/prerequisites/zlib-1.2.11\">>~/.bashrc && source ~/.bashrc && ";
#java running enrironmen 8(jdk1.8.0_121.tar.gz),GATK3.6
$cmd .="cd $outdir/tmp && wget http://download.oracle.com/otn-pub/java/jdk/8u121-b13/e9e7ea248e2c4826b92b3f075a80e441/jdk-8u121-linux-x64.tar.gz?AuthParam=1491548627_bb5d3a04d6d93aa6d3bf93322f75fba7 && tar -zxvf jdk-8u121-linux-x64.tar.gz && rm jdk*tar.gz && mv jdk1.8.0_121 $outdir/opt/bin && cd $outdir/opt/bin/jdk1.8.0_121 && echo \"export JAVA_HOME=$outdir/opt/bin/jdk1.8.0_121\" >>~/.bashrc && echo \"export JRE_HOME=\${JAVA_HOME}/jre\" >>~/.bashrc && echo \"export CLASSPATH=.:\${JAVA_HOME}/lib:\${JRE_HOME}/lib\" >>~/.bashrc echo \"export PATH=\${JAVA_HOME}/bin:\$PATH\" >>~/.bashrc && ";
# mutect1.1.5(required maven >3.0 java 7)
$cmd .="apt-get install maven && \\\n";
$cmd .="wget http://download.oracle.com/otn-pub/java/jdk/7u79-b15/jre-7u79-linux-x64.tar.gz && tar -zxvf jre-7u79-linux-x64.tar.gz && mv jre1.7.0* $outdir/opt/bin/ && \\\n";

#picard install 1.54
$cmd .="wget ";




##bioperl and Bio::DB::Sam package


print "			ALL Install Command :\n\n";
print $cmd;
my $clean;
$clean ="echo -e \"Clean....\" && cd $outdir && rm -r tmp && echo -e \"....Done.\"";
print $clean;

print "\n\n\nSuccessfully installed.\n"


