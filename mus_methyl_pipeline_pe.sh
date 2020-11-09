#!/bin/sh

if [ $# -ne 4 ]
then
    echo 'Usage: sh xxx.sh <IN sample> <IN R1.fq.gz> <IN R2.fq.gz> <IN directory>'
    exit 0
fi

sample=$1
R1=$2
R2=$3
dir=$4

###################################################1-Trimming####################################################

echo "Methylome NGS data analysis pipeline"
echo $(date +"%Y-%m-%d %T")
echo

mkdir $dir
cd $dir
mkdir trimming
cd trimming
$trim_galore --path_to_cutadapt $cutadapt --quality 20 --phred33 --stringency 3 --length 36 --rrbs --paired --trim1 --output_dir $dir/trimming $R1 $R2 >trimming.log 2>&1

echo Finished trimming
echo $(date +"%Y-%m-%d %T")
echo

#############################################2-Calculate convertion##############################################

cd $dir
mkdir bismark_lambdaDNA
cd bismark_lambdaDNA
$bismark --path_to_bowtie $bowtie2 --quiet -o $dir/bismark_lambdaDNA --temp_dir $dir/bismark_lambdaDNA/tmp  --genome_fold $lambda_DNA -1 $dir/trimming/*_1.fq.gz -2 $dir/trimming/*_2.fq.gz >./bismark.log 2>&1

echo Finished calculating convertion
echo $(date +"%Y-%m-%d %T")
echo

#################################################3-Alignment (mouse)##############################################

cd $dir
mkdir bismark_mouse
cd bismark_mouse
$bismark --path_to_bowtie $bowtie2 --quiet --un -o $dir/bismark_mouse --temp_dir $dir/bismark_mouse/tmp  --genome_fold $mus_GRCm38 -1 $dir/trimming/*_1.fq.gz -2 $dir/trimming/*_2.fq.gz >./bismark.log 2>&1

$bismark --path_to_bowtie $bowtie2 --quiet -o $dir/bismark_mouse --temp_dir $dir/bismark_mouse/tmp-un-1 --genome_fold $mus_GRCm38 $dir/bismark_mouse/*_unmapped_reads_1.fq.gz >./bismark-un-1.log 2>&1

$bismark --path_to_bowtie $bowtie2 --quiet -o $dir/bismark_mouse --temp_dir $dir/bismark_mouse/tmp-un-2 --genome_fold $mus_GRCm38 $dir/bismark_mouse/*_unmapped_reads_2.fq.gz >./bismark-un-2.log 2>&1

echo Finished alignment
echo $(date +"%Y-%m-%d %T")
echo

################################################4-Sort and merge################################################

cd $dir/bismark_mouse
samtools sort $dir/bismark_mouse/*_bismark_bt2_pe.bam -o ${sample}.pair.sort.bam 1>samtools_sort_pair.log 2>&1
samtools sort $dir/bismark_mouse/*_unmapped_reads_1_bismark_bt2.bam -o ${sample}.unpair-1.sort.bam 1>samtools_sort_un_1.log 2>&1
samtools sort $dir/bismark_mouse/*_unmapped_reads_2_bismark_bt2.bam -o ${sample}.unpair-2.sort.bam 1>samtools_sort_un_2.log 2>&1
samtools merge ${sample}.bam ${sample}.pair.sort.bam ${sample}.unpair-1.sort.bam ${sample}.unpair-2.sort.bam 
samtools sort ${sample}.bam -o ${sample}.sort.bam 1>samtools_sort.log 2>&1

echo Finished merging and sorting
echo $(date +"%Y-%m-%d %T")
echo

##################################################5-MarkDuplicates################################################

cd $dir/bismark_mouse
java -jar $picard MarkDuplicates REMOVE_DUPLICATES=false INPUT=${sample}.sort.bam OUTPUT=${sample}.sort.dedup.bam METRICS_FILE=${sample}.sort.dedup.metrics

echo Finished removing duplicates
echo $(date +"%Y-%m-%d %T")
echo

###################################################6-statistics###################################################

cd $dir/bismark_mouse
samtools index ${sample}.sort.dedup.bam

bedtools genomecov -ibam ${sample}.sort.dedup.bam -bga -trackline 1>depth_output
perl $perl_cal_cov depth_output >covered.stat
rm depth_output

echo Finished coverage statistics
echo $(date +"%Y-%m-%d %T")
echo

##########################################7-Extract the methylation call##########################################

cd $dir
mkdir methylation_calls
cd methylation_calls
$bismark_methylation_extractor --bedGraph $dir/bismark_mouse/${sample}.sort.bam

echo Finished methylation calling
echo $(date +"%Y-%m-%d %T")
echo

###########################################8-Calculate methylation level##########################################

#cd $dir
#mkdir cal_met_level
#cd cal_met_level
#perl $perl_metlevel $mus_genome $dir/bismark_mouse/${sample}.pileup > ${sample}.SingleCmet

##################################################9-CpG-covered###################################################

export PATH=/datapool/usr/R-3.5.0/bin/:$PATH
cd $dir
mkdir methylKit
cd methylKit
mkdir 1-cov
cd 1-cov
Rscript $R_methylKit_1cov $sample $dir/methylKit/1-cov $dir/bismark_mouse/${sample}.sort.bam
cat *_CpG.txt|awk '{sum+=$6} END {print "Average=",sum/NR}' > ${sample}_CpG_methyl_level.txt
perl $perl_cpg_num ${sample}_CpG.txt cpg_count.txt
cd ..
mkdir 5-cov
cd 5-cov
Rscript $R_methylKit_5cov $sample $dir/methylKit/5-cov $dir/bismark_mouse/${sample}.sort.bam
cat *_CpG.txt|awk '{sum+=$6} END {print "Average=",sum/NR}' > ${sample}_CpG_methyl_level.txt
perl $perl_cpg_num ${sample}_CpG.txt cpg_count.txt

echo Finished CpG methylation calling
echo $(date +"%Y-%m-%d %T")
echo

##################################################10-CpG-annotation################################################

cd $dir
mkdir CpG-annotation
cd CpG-annotation
cut -f 2,3 ../methylKit/5-cov/${sample}_CpG.txt >${sample}-2-3
cut -f 3  ../methylKit/5-cov/${sample}_CpG.txt >${sample}-3
paste ${sample}-2-3 ${sample}-3 >${sample}.bed
sed -i -e '/base/d' ${sample}.bed
rm ${sample}-2-3
rm ${sample}-3

bedtools intersect -a ${sample}.bed -b $promoter_bed -wa |uniq >${sample}_promoter_cpg.txt
bedtools intersect -a ${sample}.bed -b $cpgisland_bed -wa |uniq >${sample}_cpgisland_cpg.txt
bedtools intersect -a ${sample}.bed -b $sine_bed -wa |uniq >${sample}_sine_cpg.txt
bedtools intersect -a ${sample}.bed -b $line_bed -wa |uniq >${sample}_line_cpg.txt
bedtools intersect -a ${sample}.bed -b $ltr_bed -wa |uniq >${sample}_ltr_cpg.txt
bedtools intersect -a ${sample}.bed -b $exon_bed -wa |uniq >${sample}_exon_cpg.txt
bedtools intersect -a ${sample}.bed -b $intron_bed -wa |uniq >${sample}_intron_cpg.txt
bedtools intersect -a ${sample}.bed -b $cpgisland_shore_bed -wa |uniq >${sample}_cpgisland_shore_cpg.txt
bedtools intersect -a ${sample}.bed -b $utr5_bed -wa |uniq >${sample}_5_utr_cpg.txt
bedtools intersect -a ${sample}.bed -b $utr3_bed -wa |uniq >${sample}_3_utr_cpg.txt

echo Finished detected CpG annotation
echo $(date +"%Y-%m-%d %T")
echo
