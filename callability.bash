#!/bin/bash

cd /humgen/atgu1/fs03/eminikel/050callab

# wgs.list is a list of BAM files for whole genome samples
# wes.list is a list of BAM files for whole exome samples

gatkjar=/humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar
b37ref=/humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta

b37genome=/humgen/atgu1/fs03/eminikel/050callab/b37.genome
broadexome=/humgen/gsa-hpprojects/GATK/bundle/2.8/b37/Broad.human.exome.b37.interval_list # 32950014 bp
refseq=/humgen/atgu1/fs03/eminikel/050callab/b37.refseq.genes.name2.bed
musclegenes=/humgen/atgu1/fs03/eminikel/050callab/muscle-gene-symbols.txt

agilentbaits=/humgen/atgu1/fs03/lek/resources/gatk/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.interval_list
icebaits=/seq/references/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list

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

# replace symlinks with hard links to underlying files
cat wgs.list | awk '{print "readlink -f "$1}' | bash > wgs.hard.list
cat wes.list | awk '{print "readlink -f "$1}' | bash > wes.hard.list

# clean up files from first attempt at this analysis, on May 6, 2014
mkdir old_20140506
mv cov* old_20140506

mkdir jobtemp

for ilist in {be.bed,be10.bed,bm.bed,bm10.bed}
do
    for bamlist in {wgs.hard.list,wes.hard.list}
    do
        for minbq in {0,1,10,20}
        do
            for minmq in {0,1,10,20}
            do
                bsub -q bweek -P $RANDOM -J callab -M 8000000 \
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

# for re-submitting individual jobs
# bamlist=wgs.list
# ilist=be.bed
# minbq=10
# minmq=0
# bsub -q bweek -W 60:00 -P $RANDOM -J callab -M 8000000 \
#      -o jobtemp/job.$ilist.$bamlist.$minbq.$minmq.out \
#      -e jobtemp/job.$ilist.$bamlist.$minbq.$minmq.err \
#     "java -Xmx8g -jar $gatkjar \
#          -R $b37ref \
#          -T DepthOfCoverage \
#          -o cov_${ilist}_${bamlist}_${minbq}_${minmq} \
#          -I $bamlist \
#          -L $ilist \
#          --omitDepthOutputAtEachBase \
#          --omitIntervalStatistics \
#          --omitPerSampleStats \
#          --minBaseQuality $minbq \
#          --minMappingQuality $minmq"

# try to get more reproducibility
mkdir -p 2
mkdir -p 2/jobtemp
for ilist in {be.bed,be10.bed,bm.bed,bm10.bed}
do
    for bamlist in {wgs.list,wes.list}
    do
        for minbq in {0,1,10,20}
        do
            for minmq in {0,1,10,20}
            do
                bsub -q bweek -W 60:00 -P $RANDOM -J callab -M 8000000 \
                    -o 2/jobtemp/job.$ilist.$bamlist.out \
                    -e 2/jobtemp/job.$ilist.$bamlist.err \
                    "java -Xmx8g -jar $gatkjar \
                         -R $b37ref \
                         -T DepthOfCoverage \
                         -o 2/cov_${ilist}_${bamlist}_${minbq}_${minmq} \
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



####
# merge output into single file

# cumulative coverage
echo -en "ilist\tbamlist\tminbq\tminmq\tnsamp" > all_cc.txt
# now fuse everything into one big table
cat *coverage_counts | head -1 >> all_cc.txt
for ilist in {be.bed,be10.bed,bm.bed,bm10.bed}
do
    for bamlist in {wgs.list,wes.list}
    do
        for minbq in {0,1,10,20}
        do
            for minmq in {0,1,10,20}
            do
                fname=cov_${ilist}_${bamlist}_${minbq}_${minmq}.sample_cumulative_coverage_counts
                if [ -f $fname ];
                    then
                        cat cov_${ilist}_${bamlist}_${minbq}_${minmq}.sample_cumulative_coverage_counts | tail -n +2 | awk -v ilist="$ilist" -v bamlist="$bamlist" -v minbq="$minbq" -v minmq="$minmq" '{print ilist "\t" bamlist "\t" minbq "\t" minmq "\t" $0}' >> all_cc.txt
                    else
                        echo $fname "does not exist"
                fi
            done
        done
    done
done


# cumulative proportions
echo -en "ilist\tbamlist\tminbq\tminmq\tsid" > all_cp.txt
# now fuse everything into one big table
cat *coverage_proportions | head -1 >> all_cp.txt
for ilist in {be.bed,be10.bed,bm.bed,bm10.bed}
do
    for bamlist in {wgs.list,wes.list}
    do
        for minbq in {0,1,10,20}
        do
            for minmq in {0,1,10,20}
            do
                cat cov_${ilist}_${bamlist}_${minbq}_${minmq}.sample_cumulative_coverage_proportions | tail -n +2 | awk -v ilist="$ilist" -v bamlist="$bamlist" -v minbq="$minbq" -v minmq="$minmq" '{print ilist "\t" bamlist "\t" minbq "\t" minmq "\t" $0}' >> all_cp.txt
            done
        done
    done
done



# ############
# find variants called in only exome or only genome
# SelectVariants for muscle exons only, then VariantsToTable to analyze with SQL/R

# generate VCF tables to get aggregate statistics for all exons (incl non-muscle)

java -Xmx8g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   -V $wesvcf \
   -L $broadexome \
   --excludeNonVariants \
   -sn 15E_DD_1 \
   -sn 16E_MD_1 \
   -sn 17E_PD_1 \
   -sn 23H_LM_1 \
   -sn 24H_CM_1 \
   -sn 25H_JM_1 \
   -sn 65T_CR_1 \
   -sn 66T_NG_1 \
   -sn 67T_SR_1 \
   -sn 68T_DR_1 \
   -sn 69T_GG_1 \
   -o wes.vcf
# duration: 395s

java -Xmx8g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   -V $wgsvcf \
   -L $broadexome \
   --excludeNonVariants \
   -sn 15E_DD_1 \
   -sn 16E_MD_1 \
   -sn 17E_PD_1 \
   -sn 23H_LM_1 \
   -sn 24H_CM_1 \
   -sn 25H_JM_1 \
   -sn 65T_CR_1 \
   -sn 66T_NG_1 \
   -sn 67T_SR_1 \
   -sn 68T_DR_1 \
   -sn 69T_GG_1 \
   -o wgs.vcf
# 214s

java -Xmx8g -jar $gatkjar \
    -R $b37ref \
    -T VariantsToTable \
    --splitMultiAllelic \
    -V wes.vcf \
    -o wes.table \
    --variant_index_type LINEAR \
    --allowMissingData \
    --showFiltered \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F DP \
    -GF GT -GF AD -GF DP -GF GQ -GF PL

java -Xmx8g -jar $gatkjar \
    -R $b37ref \
    -T VariantsToTable \
    --splitMultiAllelic \
    -V wgs.vcf \
    -o wgs.table \
    --variant_index_type LINEAR \
    --allowMissingData \
    --showFiltered \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F DP \
    -GF GT -GF AD -GF DP -GF GQ -GF PL


cat wgs.list | awk -F"/" '{print "WGS_"$5"\t"$0}' > samples.txt
cat wes.list | awk -F"/" '{print "WES_"$5"\t"$0}' >> samples.txt

# run R script to generate sample-locus file
callability.r

# generate screenshots of all discordant site-sample combinations in WGS and WES in muscle exons
mkdir -p igv3
/home/unix/mlek/bin/IGV_plotter/IGV_plotter -LOCUS_SAMPLES sample_locus.txt -SAMPLES_TO_BAMS samples.txt -SNAPSHOT_DIR igv3



# do equivalent for only the muscle targets
# didn't end up talking about this part in the blog post.

cd /humgen/atgu1/fs03/eminikel/050callab
java -Xmx8g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   -V $wesvcf \
   -L broad-muscle-targets.bed \
   --excludeNonVariants \
   -sn 15E_DD_1 \
   -sn 16E_MD_1 \
   -sn 17E_PD_1 \
   -sn 23H_LM_1 \
   -sn 24H_CM_1 \
   -sn 25H_JM_1 \
   -sn 65T_CR_1 \
   -sn 66T_NG_1 \
   -sn 67T_SR_1 \
   -sn 68T_DR_1 \
   -sn 69T_GG_1 \
   -o wesm.vcf

java -Xmx8g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   -V $wgsvcf \
   -L broad-muscle-targets.bed \
   --excludeNonVariants \
   -sn 15E_DD_1 \
   -sn 16E_MD_1 \
   -sn 17E_PD_1 \
   -sn 23H_LM_1 \
   -sn 24H_CM_1 \
   -sn 25H_JM_1 \
   -sn 65T_CR_1 \
   -sn 66T_NG_1 \
   -sn 67T_SR_1 \
   -sn 68T_DR_1 \
   -sn 69T_GG_1 \
   -o wgsm.vcf


java -Xmx8g -jar $gatkjar \
    -R $b37ref \
    -T VariantsToTable \
    --splitMultiAllelic \
    -V wesm.vcf \
    -o wesm.table \
    --variant_index_type LINEAR \
    --allowMissingData \
    --showFiltered \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F DP \
    -GF GT -GF AD -GF DP -GF GQ -GF PL

java -Xmx8g -jar $gatkjar \
    -R $b37ref \
    -T VariantsToTable \
    --splitMultiAllelic \
    -V wgsm.vcf \
    -o wgsm.table \
    --variant_index_type LINEAR \
    --allowMissingData \
    --showFiltered \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F DP \
    -GF GT -GF AD -GF DP -GF GQ -GF PL



# create a minimal example to prove non-monotonicity of GATK DepthOfCoverage's output over minmq
bamlist=wgs.list
ilist=bm10.bed
minbq=1
for minmq in {0..1}
do
  bsub -q bhour -W 00:45 -P $RANDOM -J callab -M 2000000 \
       -o nonmon3/job.$ilist.$bamlist.$minbq.$minmq.out \
       -e nonmon3/job.$ilist.$bamlist.$minbq.$minmq.err \
      "java -Xmx8g -jar $gatkjar \
           -R $b37ref \
           -T DepthOfCoverage \
           -o nonmon3/cov_${ilist}_${bamlist}_${minbq}_${minmq} \
           -I $bamlist \
           -L $ilist \
           --omitDepthOutputAtEachBase \
           --omitIntervalStatistics \
           --omitPerSampleStats \
           --minBaseQuality $minbq \
           --minMappingQuality $minmq"
done

# for nonmon and nonmon2, the -M was 8000000. for nonmon3 it was 2000000

