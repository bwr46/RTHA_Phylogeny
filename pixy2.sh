# have to index vcf first

tabix Allsites_filtered.vcf.gz

# then run pixy

source /programs/miniconda3/bin/activate pixy

pixy --stats pi fst dxy --vcf Allsites_filtered.vcf.gz --populations rthapop.txt --window_size 10000 --n_cores 12