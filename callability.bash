#!/bin/bash

cd /humgen/atgu1/fs03/eminikel/050callab

# wgs.list is a list of BAM files for whole genome samples
# wes.list is a list of BAM files for whole exome samples

gatkjar=/humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar
b37ref=/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta

b37genome=/humgen/atgu1/fs03/eminikel/050callab/b37.genome
broadexome=/humgen/gsa-hpprojects/GATK/bundle/2.8/b37/Broad.human.exome.b37.interval_list
refseq=/humgen/atgu1/fs03/eminikel/050callab/b37.refseq.genes.name2.bed
musclegenes=/humgen/atgu1/fs03/eminikel/050callab/muscle-gene-symbols.txt

# prep the raw BED inputs for muscle disease genes and the Broad exome
\grep -w -f $musclegenes $refseq > muscle-genes-refseq.bed
\grep -v ^@ $broadexome | cut -f1,2,3  > broadexome.bed

# create 2 bed files, for muscle disease exons and a 10bp buffer around them
bedtools intersect -a muscle-genes-refseq.bed -b broadexome.bed > broad-muscle-targets-raw.bed
sortBed -i broad-muscle-targets-raw.bed > broad-muscle-targets-sorted.bed
bedtools merge -i broad-muscle-targets-sorted.bed > broad-muscle-targets.bed # merge any overlapping entries
bedtools slop -i broad-muscle-targets.bed -g $b37genome -b 10 > broad-muscle-targets-plus10.bed
bedtools merge -i broad-muscle-targets-plus10.bed > broad-muscle-targets-plus10-merged.bed
bedtools subtract -a broad-muscle-targets-plus10-merged.bed -b broad-muscle-targets.bed > broad-muscle-slop10only.bed

# create 2 bed files, for all exons and a 10bp buffer around them
bedtools merge -i broadexome.bed > broadexome-merged.bed
bedtools sort -i broadexome-merged.bed > broad-exome-targets.bed
bedtools slop -i broad-exome-targets.bed -g $b37genome -b 10 > broad-exome-targets-plus10.bed
bedtools merge -i broad-exome-targets-plus10.bed > broad-exome-targets-plus10-merged.bed
bedtools subtract -a broad-exome-targets-plus10-merged.bed -b broad-exome-targets.bed > broad-exome-slop10only.bed

cp broad-exome-targets.bed be.bed
cp broad-exome-slop10only.bed be10.bed
cp broad-muscle-targets.bed bm.bed
cp broad-muscle-slop10only.bed bm10.bed

mkdir jobtemp

for ilist in {be.bed,be10.bed,bm.bed,bm10.bed}
do
    for bamlist in {wgs.list,wes.list}
    do
        for minbq in {0,} #1,10,20}
        do
            for minmq in {0,} #1,10,20}
            do
                bsub -q bweek -W 20:00 -P $RANDOM -J callab -M 8000000 \
                    -o jobtemp/job.$ilist.$bamlist.out \
                    -e jobtemp/job.$ilist.$bamlist.err \
                    "java -Xmx8g -jar $gatkjar \
                         -R $b37ref \
                         -T DepthOfCoverage \
                         -o cov_${ilist}_${bamlist}_${minbq}_${minmq} \
                         -I $bamlist \
                         -L $ilist \
                         --omitDepthOutputAtEachBase \
                         --omitIntervalStatistics \
                         --omitPerSampleStats \
                         --minBaseQuality $minbq \
                         --minMappingQuality $minmq"
            done
        done
    done
done

# longest one took 29511 sec = 8.2h
# one failed. 
# $ cat jobtemp/job.be.bed.wgs.list.err | tail -2
# ##### ERROR MESSAGE: File /seq/picard_aggregation/G34053/68T_DR_1/current/68T_DR_1.bai is malformed: Premature end-of-file while reading BAM index file /seq/picard_aggregation/G34053/68T_DR_1/current/68T_DR_1.bai. It's likely that this file is truncated or corrupt -- Please try re-indexing the corresponding BAM file.
# ##### ERROR ------------------------------------------------------------------------------------------
# need to redo wgs be
bamlist=wgs.list
ilist=be.bed
minbq=0
minmq=0
bsub -q bweek -W 40:00 -P $RANDOM -J callab -M 8000000 \
     -o jobtemp/job.$ilist.$bamlist.$minbq.$minmq.out \
     -e jobtemp/job.$ilist.$bamlist.$minbq.$minmq.err \
    "java -Xmx8g -jar $gatkjar \
         -R $b37ref \
         -T DepthOfCoverage \
         -o cov_${ilist}_${bamlist}_${minbq}_${minmq} \
         -I $bamlist \
         -L $ilist \
         --omitDepthOutputAtEachBase \
         --omitIntervalStatistics \
         --omitPerSampleStats \
         --minBaseQuality $minbq \
         --minMappingQuality $minmq"



for ilist in {be.bed,be10.bed,bm.bed,bm10.bed}
do
    for bamlist in {wgs.list,wes.list}
    do
        for minbq in {1,10,20}
        do
            for minmq in {1,10,20}
            do
                bsub -q bweek -W 20:00 -P $RANDOM -J callab -M 8000000 \
                    -o jobtemp/job.$ilist.$bamlist.$minbq.out \
                    -e jobtemp/job.$ilist.$bamlist.$minmq.err \
                    "java -Xmx8g -jar $gatkjar \
                         -R $b37ref \
                         -T DepthOfCoverage \
                         -o cov_${ilist}_${bamlist}_${minbq}_${minmq} \
                         -I $bamlist \
                         -L $ilist \
                         --omitDepthOutputAtEachBase \
                         --omitIntervalStatistics \
                         --omitPerSampleStats \
                         --minBaseQuality $minbq \
                         --minMappingQuality $minmq"
            done
        done
    done
done



echo -en "ilist\tbamlist\tminbq\tminmq\tnsamp" > all_cc.txt
# now fuse everything into one big table
cat *coverage_counts | head -1 >> all_cc.txt
for ilist in {be.bed,be10.bed,bm.bed,bm10.bed}
do
    for bamlist in {wgs.list,wes.list}
    do
        for minbq in {0,}
        do
            for minmq in {0,}
            do
                cat cov_${ilist}_${bamlist}_${minbq}_${minmq}.sample_cumulative_coverage_counts | tail -n +2 | awk -v ilist="$ilist" -v bamlist="$bamlist" -v minbq="$minbq" -v minmq="$minmq" '{print ilist "\t" bamlist "\t" minbq "\t" minmq "\t" $0}' >> all_cc.txt
            done
        done
    done
done
