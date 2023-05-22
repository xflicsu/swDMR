
swDMR: A sliding window approach to identify differentially methylated regions from whole genome bisulfite sequencing.
https://www.ncbi.nlm.nih.gov/pubmed/26176536
https://sourceforge.net/projects/swdmr/
Software discription
 
   swDMR is a sliding window approach, which is used to identify differentially methylated regions (DMR) from bisulfite sequencing dataset at single base resolution. The dataset, like whole genome bisulfite sequencing (WGBS) or reduced representation bisulfite sequencing (RRBS) with the same coverage region of samples, are suitable for swDMR. This software integrated several useful statistics methods, satisfying to two or multiple samples test. In addition, swDMR provides genomic features annotation with BED or GFF file. It can also produce WIG format file to upload to UCSC genome browser for DMR visualization.

Platform and environment
    
 System: Linux or UNIX
 Software: perl >= v5.8.8, R >= 2.12.0 and BEDTools are required.
 
How to run swDMR

1. Install

>cd /path/of/swDMR----- Go to path of swDMR software installed
>tar xzvf swDMR-xxx.tar.gz
>cd swDMR-xxx
>sh install
### tips:
	If you see some warns, when you install the dependence packages of swDMR.
	It may be that the packages in our swDMR are too old for your system. You
	can download the latest packages for installation as the file 'install'.
	
	Sorry for the problem! We will solve this in the future.

2. Run the software

>cd swDMR-xxx
>./swDMR --help

Discription: This program is used to detect DMRs, which can also annotate DMRs through user's defination of BED format file.
Usage:
  Input methylation file format:
  must contain follow information (separated by table "\t")
  you can assign the corresponding information contained in your files
 
		information required: chromosome position ctype methyl_reads unmethy_reads
 
 Example:
 
 ./swDMRs [options] < ... >
 swDMR --samples ESFibro.chr1.cout.gz,Fibro.chr1.cout.gz --name ESFibro,Fibro --outdir outdir --statistics Fisher --cytosineType CG --window 1000 --stepSize 100 --length 100 --pvalue 0.01 --coverage 4 --fold 2 --fdr 0.05 --Rbin /usr/bin/R --chromosome 1 --position 2 --ctype 4 --methy 7 --unmethy 8 -a UTR_CDS_Intron_Upstream_Downstream.bed -g gene.bed --CGI cpgIslandExt.txt

 Basic swDMR options:
  -h | --help	print this help information on screen
  -sam | --samples	input samples' methylation files separated by	comma "," < a,b,c,...,g,h >
  -na | --name	input samples' names separated by comma "," < a,b,c,...,g,h >
  -chp | --chromosome	column of chromosome in your file, default <1>
  -pos | --position	column of position in your file, default <2>
  -cp | --ctype	column of cytosine type in your file, default <4>
  -mp | --methy	column of methy reads in your file, default <5>
  -up | --unmethy	column of unmethy reads in your file, default <6>
  -o | --outdir	swDMR result directory
  -co | --coverage	lowest coverage of cytosine reads to use, default <4>
  -s | --statistics	choose one method to detect DMRs
		< T_test||Kolmogorov||Fisher||ChiSquare||Wilcoxon||ANOVA||Kruskal >
        if input samples == 2 --statistics should choose
		< T_test||Kolmogorov||Fisher||ChiSquare||Wilcoxon >
        if input samples >= 3 -statistics could choose	
		< ANOVA||Kruskal >
  -CT | --cytosineType	cytosine type < C|CG||CHG||CHH >
  -w | --window	a sliding window to statistics, default <1000>
  -N | --points	lowest number of selected type of cytosine in the window, default <10>
  -z | --stepSize	step size of the sliding processes, default <100>
  -f | --fold 	max/min methylation level difference
  -d | --diff 	value of max-min methylation level
  -len | --length	lowest length to join two fragment into one, default <100>
  -V | --pvalue	p value to judge as a DMR, default <0.01>
  -R | --Rbin	R bin
  -fdr	fdr to adjust DMR p value, default <0.05>
  -pro | --processes	  parallel processes for DMR test and annotation, default <1>
 Annotation and DMR maps:
  -a | --annotation	your annotation file should be bed of GFF format < BED||GFF >
  -g | --Gene	BED file of gene coordinates
  --CGI	cpgIsland file:
	http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cpgIslandExt.txt.gz. 
	The other release or species cpgIsland file can also be downloaded from UCSC genome browser
  --left	left side length of DMR in a map, default <1000>
  --right	right side length of DMR in a map, default <1000>
  
  
  
## output fomart
https://github.com/xflicsu/swDMR/wiki/About-swDMR#output-of-swdmr

