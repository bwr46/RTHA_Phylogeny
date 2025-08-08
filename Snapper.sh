#Divergence time estimation in SNAPPER


#First need to filter and prep .vcf. We'll use the same vcf we filtered for IQTree


#IQTreeData.vcf snps are biallelic, now need to be thinned. 

vcftools --vcf IQTreeData.vcf --thin 1000 --out IQTreeData_thinned

#OR

#We could just use the fixed.phy from our IQTree run...

#From first run, we found out one of the fumosus specimens has a lot of missing data, which caused a misplacement for that group. Need to remove that individual (h4). Did so manually by opening the .phy file and deleting the sequences for h4, then changing the sample size in the header from 73 to 72. 


#Download ruby script for prep from Github

wget https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/snapp_prep.rb


#Need a populations file specifying which population each individuals belongs to (should already be created, but make sure headers are 'species' and 'individual')

#Also need a constraint file that sets the divergence time estimate (use the divergence time for RTHA and SWHA/RLHA, which from Timetree is estimated at 3.1 mya


#Put the following in a .txt file named constraints.txt. All taxa in the pop list need to be included here (i.e. include lagopus and swainsoni)

normal(0,3.1,0.74)	crown	abieticola, alascensis, borealis, calurus, costaricensis, fuertesi, fumosus, hadropus, harlani, jamaicensis, kriderii, kemsiesi, lagopus, socorroensis, solitudinus, suttoni, swainsoni, umbrinus, ventralis


#Time to run the prep script downloaded above

ruby snapp_prep.rb -a SNAPPER -p IQTreeData.fixed.NOh4.phy -t rthapop.txt -c constraints.txt -m 1000 -l 1000000


#try using a lot of threads (512 possible on 256 core large machine on biohpc)

/programs/beast-2.7.6/bin/beast -threads 500 snapper.xml

export BEAST_PACKAGE_PATH=/home/bwr46/.beast/2.7.6

#now annotate trees

/programs/beast-2.7.6/bin/treeannotator -burnin 10 -height mean snapper.trees snapper.tre










