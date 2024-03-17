# convert bam to bed file
cd /lustre/home/fhliu/fasta_files/WMZ/RUN_0167/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa/
bowtie2 -p 32  -N 1 --local -x /lustre/home/fhliu/HiCPro_data/bowtie2-build/mm9_noYM/mm9/mm9_noYM -1 wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_raw_1.fq.gz -2 wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_raw_2.fq.gz --al-conc-gz /lustre/home/fhliu/CHIP_seq/WMZ/RUN_0167/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9al.fq.gz --un-conc-gz /lustre/home/fhliu/CHIP_seq/WMZ/RUN_0167/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9un.fq.gz -S /lustre/home/fhliu/CHIP_seq/WMZ/RUN_0167/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p.sam
bowtie2 -p 32  -N 1 --local -x /lustre/home/fhliu/genome_reference/bowtie2_build/hg19/hg19 -1 wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_raw_1.fq.gz -2 wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_raw_2.fq.gz --al-conc-gz /lustre/home/fhliu/CHIP_seq/WMZ/RUN_0167/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgal.fq.gz --un-conc-gz /lustre/home/fhliu/CHIP_seq/WMZ/RUN_0167/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgun.fq.gz -S /lustre/home/fhliu/CHIP_seq/WMZ/RUN_0167/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp.sam

cd /lustre/home/fhliu/CHIP_seq/WMZ/RUN_0167/wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa/
#mm9
samtools view -h -f 2 -q 30 -F 256 -bS wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p.sam > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam
samtools sort wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam -o wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam.tmp
mv wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam.tmp wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam
samtools index wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam
samtools flagstat wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_stat1
samtools idxstats wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_stat2

java -jar /lustre/home/fhliu/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam O=wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_mkdup M=wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_mkdup.txt
samtools flagstat wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_mkdup > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_mkdup_stat1
#hg19
samtools view -h -f 2 -q 30 -F 256 -bS wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp.sam > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam
samtools sort wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam -o wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam.tmp
mv wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam.tmp wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam
samtools index wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam
samtools flagstat wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam_stat1
samtools idxstats wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam_stat2

java -jar /lustre/home/fhliu/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam O=wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam_mkdup M=wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_mkdup.txt
samtools flagstat wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam_mkdup > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_hgp_bam_mkdup_stat1





# convert bam to bed file
bedtools bamtobed -i wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_mkdup > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed
mv wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed_final

# (optional) remove mitochondria reads
grep -v chrM wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed_final > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed
mv wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed_final

# (optional) remove reads in blacklist regions
bedtools intersect -v -a wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed_final -b /lustre/home/fhliu/CHIP_seq/RUN_0004/mm9_blacklist.bed > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed
mv wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed_final

# revert bed back to bam
bedtools bedtobam -i wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bed_final -g /lustre/home/fhliu/HiCPro_data/genome_reference/mm9/mm9_noYM.sizes.txt > wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_final
samtools index wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_final
bamCoverage -bs 1 --normalizeUsing CPM -b wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_bam_final -o wmz-chipseq-S0493-stag1-B11-antipolii-rep1-4h_noa_mm9p_cpm.bw

