## CREATE AN ALLSITES VCF THAT INCLUDES FILTERING PERFORMED ON VARIANTS

# Everything should go through PrepForVariantCalling.sh script first...

# Script for variant calling in Samtools
# author: Bryce Robinson. Modified from K Carbeck scripts
# split by chromosome/contig

#### explore a bam file to figure out how to split up over 40 cores
samtools view -H h10.rescaled.bam_sorted_mark.bam | less

## get genome length
# grep lines that with @SQ
samtools view -H h10.rescaled.bam_sorted_mark.bam | less | grep "^@SQ" | cut -f 3 | cut -d: -f2 | awk '{sum+=$1}END{print sum}'

## get number of contigs
samtools view -H h10.rescaled.bam_sorted_mark.bam | grep "^@SQ" | cut -f 3 | cut -d: -f2 | wc -l

# make regions file
samtools view -H h10.rescaled.bam_sorted_mark.bam | grep "^@SQ" | cut -f 2,3  | awk '
BEGIN {
    cutoff=27000000
    FS=OFS="\t"
    num = 1
}
{
    split($1,ctg,":")
    split($2,ln,":")
    cumulLen += ln[2]
    print ctg[2] , "1" , ln[2] > "regions_" (num<10?"0":"") num ".txt"
    if (cumulLen > cutoff) {
        num++
        cumulLen = 0
    }
}'


# looking at mapping quality
samtools view --exclude-flags SECONDARY,SUPPLEMENTARY h10.rescaled.bam_sorted_mark.bam | cut -f5 | sort -nr | uniq -c | less


# http://www.htslib.org/doc/bcftools.html#mpileup
# Generate VCF or BCF containing genotype likelihoods for one or multiple alignment (BAM or CRAM) files. This is based on the original samtools mpileup command (with the -v or -g options) producing genotype likelihoods in VCF or BCF format, but not the textual pileup output. The mpileup command was transferred to bcftools in order to avoid errors resulting from use of incompatible versions of samtools and bcftools when using in the mpileup+bcftools call pipeline.
# Individuals are identified from the SM tags in the @RG header lines. Multiple individuals can be pooled in one alignment file, also one individual can be separated into multiple files. If sample identifiers are absent, each input file is regarded as one sample.

## mpileup
  # -O: output type. u = bcf
  # -C: Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50
  # -a: annotate. list of common INFO tags for vcf output
  # -f: Reference genome
## call
  # -multiallelic-caller: alternative model for multialleleic and rare-variant calling designed to overcome known limitations in -c calling mode (previous version)
  # -variants-only: output variant sites only
  # -O: output type. v=VCF
  # -o: output file name
  # -b = List of input alignment files, one file per line (samplename_sorted_RGmark.bam)


mkdir vcf2


## Need to remove variants only flag plus any -v

genome=rtha_genome.fa
bamfiles=*_sorted_mark.bam
for regions in regions_*.txt; do
   mpileupVCF=${regions%.txt}.raw.mpileup.vcf
    echo $mpileupVCF
    bcftools mpileup -a "FORMAT/AD,FORMAT/DP,INFO/AD,FORMAT/DV,FORMAT/DP4,FORMAT/DPR,INFO/DPR" --threads 2 -d 1000 --min-MQ 30 --regions-file $regions -Ou -f $genome $bamfiles | bcftools call --multiallelic-caller --threads 2 -m -Ov -o vcf2/$mpileupVCF 2>&1 | tee vcf2/${regions%.txt}.mpileup.$(date +%Y%m%d-%Hh%Mm%Ss).log &
    sleep 30
done


## Have to make sure the min alleles is set to 1, no MAF filter, etc. Really the appropriate step is to not run this, but instead filter variants after splitting variants and invariants

for mpileupVCF in *.raw.mpileup.vcf; do
    filteredVCF=${mpileupVCF%.raw.mpileup.vcf}.filtered.vcf
    echo $filteredVCF
    vcftools --vcf $mpileupVCF --max-missing 0.8 --min-meanDP 4 --max-meanDP 50 --minQ 35 --min-alleles 1 --max-alleles 2 --recode --out $filteredVCF 2>&1 | tee ${mpileupVCF%.raw.mpileup.vcf}.filter.$(date +%Y%m%d-%Hh%Mm%Ss).log &
    sleep 15
done

#################################################################################################
#### bcftools concatenate ####
# you can't just cat them together because each file has a header section. The bcftools command will handle that.

# first, each vcf file must be sorted prior to calling bcftools concat
for filteredVCF in *.filtered.vcf.recode.vcf; do
    sortedVCF=${filteredVCF%.filtered.vcf.recode.vcf}.sorted.filtered.vcf
    echo $filteredVCF
    bcftools sort $filteredVCF -Oz -o $sortedVCF &
done

# then, bcftools concatenate
bcftools concat regions_01.sorted.filtered.vcf regions_02.sorted.filtered.vcf regions_03.sorted.filtered.vcf regions_04.sorted.filtered.vcf regions_05.sorted.filtered.vcf regions_06.sorted.filtered.vcf regions_07.sorted.filtered.vcf regions_08.sorted.filtered.vcf regions_09.sorted.filtered.vcf regions_10.sorted.filtered.vcf regions_11.sorted.filtered.vcf regions_12.sorted.filtered.vcf regions_13.sorted.filtered.vcf regions_14.sorted.filtered.vcf regions_15.sorted.filtered.vcf regions_16.sorted.filtered.vcf regions_17.sorted.filtered.vcf regions_18.sorted.filtered.vcf regions_19.sorted.filtered.vcf regions_20.sorted.filtered.vcf regions_21.sorted.filtered.vcf regions_22.sorted.filtered.vcf regions_23.sorted.filtered.vcf regions_24.sorted.filtered.vcf regions_25.sorted.filtered.vcf regions_26.sorted.filtered.vcf regions_27.sorted.filtered.vcf regions_28.sorted.filtered.vcf regions_29.sorted.filtered.vcf regions_30.sorted.filtered.vcf -Oz -o filtered_rtha_total_samplesminQ35_ALLSITES.vcf.gz

# valadate vcf to make sure everything is ok

vcf-validator filtered_rtha_total_samplesminQ35_ALLSITES.vcf.gz

#15,767,149 variants



# SPLIT VCF, THEN FILTER VARIANTS, THEN CONCATENATE

#!/bin/bash
# requires bcftools/bgzip/tabix and vcftools

# create a filtered VCF containing only invariant sites
vcftools --gzvcf filtered_rtha_total_samplesminQ35_ALLSITES.vcf.gz \
--max-maf 0 \
--recode --stdout | bgzip -c > test_invariant.vcf.gz

# create a filtered VCF containing only variant sites
vcftools --gzvcf filtered_rtha_total_samplesminQ35_ALLSITES.vcf.gz  \
--remove-indels \
--maf 0.05 \
--min-alleles 2 --max-alleles 2 \
--recode --stdout | bgzip -c > filtered_rtha_total_samplesminQ35_variants_only.vcf.gz

# index both vcfs using tabix
tabix test_invariant.vcf.gz
tabix test_variant.vcf.gz

# combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
test_variant.vcf.gz test_invariant.vcf.gz \
-O z -o Allsites_filtered.vcf.gz

