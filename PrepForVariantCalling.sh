

##Going to use sambamba to sort
#https://www.basepairtech.com/blog/sorting-bam-files-samtools-vs-sambamba/
# loop for Sambamba sort
[ -f /workdir/bwr46/rescaled/SambambaSortCommands.txt ] && rm /workdir/bwr46/rescaled/SambambaSortCommands.txt

cd /workdir/bwr46/rescaled/

for BAM in *.bam; do
  OUT=/workdir/bwr46/rescaled/${BAM}_sorted.bam

  #create sort file
  echo "/programs/sambamba-0.7.1/sambamba sort -t 30 -m 45G -o $OUT $BAM" >> /workdir/bwr46/rescaled/SambambaSortCommands.txt

done

parallel -j 4 < /workdir/bwr46/rescaled/SambambaSortCommands.txt


# Script for prepping mapped reads for variant caller
# author: Jen Walsh. Adapted from katherine carbeck

# A new Bam file will be produced for each step. In the end, you only need to keep samplename_sortedRGmark.bam

# Reserved a medium machine (24 cores)
# Bam files from mapping script above are indexed and sorted already
# copy into workdir:
\ls *.bam | parallel -j 20 cp -a {} /workdir/bwr46 &



#### Prepare the genome: index it (fai and dict files)
# R:Refernece
# O:Output
java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar CreateSequenceDictionary R=rtha_genome.fa O=rtha_genome.fa.dict
samtools faidx rtha_genome.fa 

#Traditionally, I would add readgroup info, but I have SM info in the bam files, which is used for labeling individuals in VCF and am fine with that. I have no other use for RG info for this project

###### mark duplicates - identifies duplicate reads: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
# Metrics File: File name to write duplicate metrics to
# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: maximum number of file handles to keep open when spilling read ends to a desk. keep this set at 1000
# uses a lot of memory (could probably run more than 8 at a time if I decreased -Xmx to ~10g?)

for RGBAM in *_sorted.bam; do

  #string replacement command
  OUT=${RGBAM/%_sorted.bam/_sorted_mark.bam}
  METRICS=${RGBAM/%_sorted.bam/.metrics.txt}

  #create index file
  echo "java -Xmx15g -jar /programs/picard-tools-2.8.2/picard.jar MarkDuplicates INPUT=$RGBAM OUTPUT=$OUT METRICS_FILE=$METRICS MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000" >> /workdir/bwr46/rescaled/markDuplicatesCommands.txt

done

parallel -j 8 < /workdir/bwr46/rescaled/markDuplicatesCommands.txt



#### Validate files
# since running HaplotypeCaller, don't need to realign or fix mates unless there is an error in ValidateSamFile

[ -f /workdir/bwr46/rescaled/validateCommands.txt ] && rm /workdir/bwr46/rescaled/validateCommands.txt

cd /workdir/bwr46/rescaled

for MARKBAM in *mark.bam; do

  OUTFILE=${MARKBAM/%.bam/_validate}

  echo "java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar ValidateSamFile I=$MARKBAM OUTPUT=$OUTFILE MODE=SUMMARY" >> /workdir/bwr46/rescaled/ValidateCommands.txt

done

parallel -j 22 < /workdir/bwr46/rescaled/ValidateCommands.txt

# concatenate all validate output into one text file to see if there are any errors
cat *validate > summary_validate.txt
#No errors


#### index .bam files for HaplotypeCaller
#can run multiple at a time
[ -f /workdir/bwr46/rescaled/IndexCommands.txt ] && rm /workdir/bwr46/rescaled/IndexCommands.txt

cd /workdir/bwr46/rescaled

for MARKBAM in *mark.bam; do

  #create index file
  echo "java -Xmx10g -jar /programs/picard-tools-2.8.2/picard.jar BuildBamIndex I=$MARKBAM" >> /workdir/bwr46/rescaled/IndexCommands.txt

done

parallel -j 22 < /workdir/bwr46/rescaled/IndexCommands.txt

