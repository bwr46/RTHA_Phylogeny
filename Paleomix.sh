## Script for running paleomix on server

1. Create a config file

# make directory

mkdir $HOME/.paleomix/

# add config file (bam_pipeline.ini) to above directory

# In your working directory (workdir/USER NAME) - store the ref and fastq.gz files as well as the yaml files on a ROCKY 9 server

# Run the following before running Paleomix script

export PALEOMIX="singularity run --bind /workdir/$USER --pwd /workdir/$USER /programs/paleomix-1.3.6/paleomix.sif"

#Now run the paleomix command on the indivdual .yaml file in a tmux session. Max out the machine, but be sure to leave room for variation in how many cores particular parts of the .yaml file takes.

$PALEOMIX paleomix bam_pipeline run example.yaml &




