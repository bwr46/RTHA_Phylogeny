#IQTree for phylogenetic analysis - Bryce Robinson 

# Convert to phylip format

##Going to do a series of filtering steps to reduce the number of SNPs and individuals before running IQTree

##25% missing data, 3X coverage, linkage

vcftools --gzvcf filtered_rtha_total_samplesminQ35_variants_only.vcf.gz --max-missing 0.75 --min-meanDP 3 --max-meanDP 50 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --out IQTreeData
# kept XXX out of a possible 8495350 Sites

##Filtering data and making input files

## need to run this before tab file

export PERL5LIB=/programs/vcftools-v0.1.16/share/perl5/5.32:$PERL5LIB 

##Make tab file
cat IQTreeData.recode.vcf | vcf-to-tab > IQTreeData.tab

##Swap out indels for missing data
##Scripts are frorm L. Campagna. I am not including these on GitHub, as they aren't mine. But email me if you need perl scripts (jlw395@cornell.edu)


##Before running script below, open fix_tab_file.pl in a text editor and change the name/path for your input and output file names
perl fix_tab_file.pl

##Write Fasta File. 
perl vcf_tab_to_fasta_alignment.pl -i IQTreeData_fixed.tab > IQTreeData.fixed.fasta

##convert fasta file to phylip format.
perl fasta2relaxedPhylip.pl -f IQTreeData.fixed.fasta -o IQTreeData.fixed.phy

# Install program

/programs/iqtree-2.2.2.6-Linux/bin/iqtree2

#Run IQTree??. It will fail??, but it will print conflicting sites to nohup. Record how many invariable sites are removed
nohup iqtree s -s IQTreeData.fixed.phy -n make_exclude_file &

##Take make_exclude_file.pl and make a file that captures those sites
perl exclude_file_maker.pl

#this second run will make a new phy file without the SNPs in exclude.txt
/programs/RAxML-8.2.4/raxmlHPC-PTHREADS-SSE3 -T 24 -f a -x 12345 -p 12345 -# 10 -m ASC_GTRGAMMA --asc-corr=lewis -E exclude.txt -s IQTreeData.fixed.phy -n excluded
0 out of 6364524 columns

#So now you have a phy with only the important positions and you're ready to run IQTree.

#Run this first, which will result in a variable sites .phy that you use for a second run

/programs/iqtree-2.2.2.6-Linux/bin/iqtree2 -s IQTreeData.fixed.phy --prefix concat -st DNA -m GTR+ASC -B 1000 

# use concat.varsites.phy in -s

/programs/iqtree-2.2.2.6-Linux/bin/iqtree2 -s concat.varsites.phy --prefix concat -st DNA -m GTR+ASC -B 1000 -nt 40

# compute site concordance factor using likelihood with v2.2.2

/programs/iqtree-2.2.2.6-Linux/bin/iqtree2 -te concat.treefile -s concord2.varsites.phy --scfl 100 --prefix concord2 -st DNA -m GTR+ASC -nt 50








