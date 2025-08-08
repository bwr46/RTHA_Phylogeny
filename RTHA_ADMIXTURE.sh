#Admixture analysis for described Red-tailed Hawk subspecies
#Original script by Jen Walsh, March 2023, adapted by Bryce Robinson for RTHA Phylogeny
#Analyses - ADMIXTURE

#Remove outgroup individuals from VCF

vcftools --gzvcf filtered_rtha_total_samplesminQ35_variants_only.vcf.gz --remove-indv MR3 --remove-indv MR4 --recode --stdout | bgzip -c > filtered_rtha_total_samplesminQ35_variants_only_nooutgroup.vcf.gz

#Remove missing data and prune for linkage

vcftools --gzvcf filtered_rtha_total_samplesminQ35_variants_only_nooutgroup.vcf.gz --max-missing 1 --recode --stdout | bgzip > filtered_rtha_total_samplesminQ35_variants_only_nooutgroup_NoMissing.vcf.gz

vcftools --gzvcf test_variant_nooutgroup.vcf.gz --max-missing 1 --recode --stdout | bgzip > test_variant_nooutgroup_NoMissing.vcf.gz

#198,726 snps retained! Something might be off here...

#Prune for linkdate
bash ldPruning.sh nocontrol_NoMissing.vcf.gz

bash ldPruning.sh test_variant_nooutgroup_NoMissing.vcf.gz
#retained 86,273 SNPs

#Prepare input files for ADMIXTURE
vcftools --gzvcf test_variant_nooutgroup_NoMissing.LDpruned.vcf.gz --plink --out RTHA.noN.LDprunedPlinkformat
plink --file RTHA.noN.LDprunedPlinkformat --make-bed

# 200 iterations for k 1-7 (which is the number of sites/populations sampled)
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; \
	do admixture --cv=5 -C=200 plink.bed $K | tee RTHA_K${K}.out; done

# -C 100
#^flag to set the number of iterations to 100

echo "CV Error:"
grep -h CV RTHA_K*.out

