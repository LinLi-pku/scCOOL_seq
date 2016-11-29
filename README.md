# scCOOL_seq
single-cell COOL-Seq analysis pipline
October 19,2016
Contents
1 Preprocessing
 1.1 Installation of softwares
 1.2 Find the positions of WCG and GCH in the genome
2 From raw fastq to calculate WCG and GCH Methylation Level and define Nucleosome Depleted Regions (NDR) and Nucleosome Occupied Regions


1 Preprocessing
1.1 Installation of softwares
Firstly, we have to install softwares for our single-cell COOL-Seq analysis. 

mkdir software
cd software

## install python anaconda 4.1.1
wget https://repo.continuum.io/archive/ Anaconda2-4.1.1-Linux-x86_64.sh
bash Anaconda2-4.1.1-Linux-x86_64.sh

## install samtools 0.1.18
wget http://jaist.dl.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2
tar jxvf samtools-0.1.18.tar.bz2
cd samtools-0.1.18
make
cd ..

## install bedtools 2.17.0
wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
tar jxvf bedtools-2.26.0.tar.gz
cd bedtools-2.26.0
make
cd ..

## instal tabix 0.2.6
wget http://jaist.dl.sourceforge.net/project/samtools/tabix/tabix-0.2.6.tar.bz2
tar jxvf tabix-0.2.6.tar.bz2
cd tabix-0.2.6
make
cd ..

## install trim_galore 0.3.3
wget http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.3.3.zip
unzip trim_galore_v0.3.3.zip
cd ..

## install bismark 0.7.6
wget http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.7.6.tar.gz
tar zxvf bismark_v0.7.6.tar.gz
cd ..

## install bowtie 1.0.0
wget http://jaist.dl.sourceforge.net/project/bowtie-bio/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip
unzip bowtie-1.0.0-linux-x86_64.zip
cd ..

Next, put the module MethGC into the path for python packages. Then revise all the paths in  path_to_MethGC/setting/projpath.py to your own paths.

1.2 Find the positions of WCG and GCH in the genome
Firstly, we need to download a copy of genome sequence from UCSC and add lambda.fa to creat a reference genome, like mm9_lambda.fa tha we can align our single-cell COOL-Seq reads to it.

mkdir database_MethGC
cd database_MethGC
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
tar zxvf chromFa.tar.gz
rm chromFa.tar.gz
for I in `seq 1 19` X Y M;
do
	cat chr”$i”.fa >>mm9.fa
done

cat mm9.fa lambda.fa >mm9_lambda.fa
samtools faidx mm9_lambda.fa

Secondly, we can use get_enzymeSite to find out all the positions of WCG and GCH in the genome. 
path_to_MethGC/bin/get_enzymeSites/get_enzymeSite –s ACG.TCG path_to_ database_MethGC /mm9_lambda.fa
path_to_MethGC/bin/get_enzymeSites/get_enzymeSite –s GCA.GCC.GCT path_to_ database_MethGC /mm9_lambda.fa

In order to speed up the calculation of methylation level of each WCG and GCH, we divide the WCG and GCH position into each chromosome.
mkdir mm9_lambda.fa.ACG.TCG mm9_lambda.fa.GCA.GCC.GCT
each chromosome has 2 bed files:
 
 

2 From raw fastq to WCG and GCH Methylation Level and defining Nucleosome Depleted Regions (NDR) and Nucleosome Occupied Regions
Now, we can go to our analysis directory and use single-cell COOL-Seq analysis pipline(MethGC) to analysis our data.

mkdir mouse_scCOOL_Seq
cd mouse_scCOOL_Seq
cp path_to_MethGC/ run_meth.py ./

We have to put our single-cell COOL-Seq fastq data into 00.0.raw_data. Each sample has a single
