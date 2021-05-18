#!/bin/bash -l

set -e
# use set -ex for debugging instead

proj=u2015018
mail=carolin.seyfferth@umu.se
in=/mnt/picea/projects/aspseq/jfelten/ERF2018
out=/mnt/picea/projects/aspseq/jfelten/ERF2018/multiqc

if [ ! -d $out ]; then
	mkdir -p $out
fi

module load bioinfo-tools multiqc

sbatch --mail-user=$mail -o $in/multiqc.out -e $in/multiqc.err \
-A $proj $UPSCb/pipeline/runMultiQC.sh $in $out
