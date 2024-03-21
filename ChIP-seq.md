
# align seq data to genome
need to change directory respectively

cd $sample/

bowtie2 -p 32  -N 1 --local -x /lustre/home/fhliu/HiCPro_data/bowtie2-build/mm9_noYM/mm9/mm9_noYM -1 sample_raw_1.fq.gz -2 sample_raw_2.fq.gz --al-conc-gz sample_mm9al.fq.gz --un-conc-gz sample_mm9un.fq.gz -S sample_mm9p.sam

bowtie2 -p 32  -N 1 --local -x /lustre/home/fhliu/genome_reference/bowtie2_build/hg19/hg19 -1 sample_raw_1.fq.gz -2 sample_raw_2.fq.gz --al-conc-gz sample_hgal.fq.gz --un-conc-gz sample_hgun.fq.gz -S sample_hgp.sam

#mm9
samtools view -h -f 2 -q 30 -F 256 -bS sample_mm9p.sam > sample_mm9p_bam
samtools sort sample_mm9p_bam -o sample_mm9p_bam.tmp
mv sample_mm9p_bam.tmp sample_mm9p_bam
samtools index sample_mm9p_bam
samtools flagstat sample_mm9p_bam > sample_mm9p_bam_stat1
samtools idxstats sample_mm9p_bam > sample_mm9p_bam_stat2

java -jar /lustre/home/fhliu/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=sample_mm9p_bam O=sample_mm9p_bam_mkdup M=sample_mm9p_mkdup.txt
samtools flagstat sample_mm9p_bam_mkdup > sample_mm9p_bam_mkdup_stat1

#hg19
samtools view -h -f 2 -q 30 -F 256 -bS sample_hgp.sam > sample_hgp_bam
samtools sort sample_hgp_bam -o sample_hgp_bam.tmp
mv sample_hgp_bam.tmp sample_hgp_bam
samtools index sample_hgp_bam
samtools flagstat sample_hgp_bam > sample_hgp_bam_stat1
samtools idxstats sample_hgp_bam > sample_hgp_bam_stat2

java -jar /lustre/home/fhliu/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=sample_hgp_bam O=sample_hgp_bam_mkdup M=sample_hgp_mkdup.txt
samtools flagstat sample_hgp_bam_mkdup > sample_hgp_bam_mkdup_stat1





# convert bam to bed file
bedtools bamtobed -i sample_mm9p_bam_mkdup > sample_mm9p_bed
mv sample_mm9p_bed sample_mm9p_bed_final

# (optional) remove mitochondria reads
grep -v chrM sample_mm9p_bed_final > sample_mm9p_bed
mv sample_mm9p_bed sample_mm9p_bed_final

# (optional) remove reads in blacklist regions
bedtools intersect -v -a sample_mm9p_bed_final -b /lustre/home/fhliu/CHIP_seq/RUN_0004/mm9_blacklist.bed > sample_mm9p_bed
mv sample_mm9p_bed sample_mm9p_bed_final

# generate bigwig using cpm normalization
bedtools bedtobam -i sample_mm9p_bed_final -g /lustre/home/fhliu/HiCPro_data/genome_reference/mm9/mm9_noYM.sizes.txt > sample_mm9p_bam_final
samtools index sample_mm9p_bam_final
bamCoverage -bs 1 --normalizeUsing CPM -b sample_mm9p_bam_final -o sample_mm9p_cpm.bw
# generate bigwig using scf normalization
With the previous alignment results of spike in sequence(sample_hgp_bam_mkdup_stat1), a scaling factor can be calculated with the formula bellow: 
scf(scaling factor) = 2000000/spike in reads
Then run the following code to get the scf normalization bigwig.
"samtools index sample_mm9p_bam_final
bamCoverage -bs 1 --scaleFactor scf(the exact number generated from the formula) -b sample_mm9p_bam_final -o sample_mm9p_scf.bw
"



