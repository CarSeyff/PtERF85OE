#!/bin/bash -l

set -ex

proj=u2015018
mail=carolin.seyfferth@umu.se
in=/mnt/picea/projects/aspseq/jfelten/ERF2018/H201SC19030195/raw_data
out=/mnt/picea/projects/aspseq/jfelten/ERF2018
genome=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/STAR/2.5.2b/Potra01
gff3=/mnt/picea/storage/reference/Populus-tremula/v1.1/gff3/Potra01-gene-synthetic-transcripts-wo-intron.gff3
gtf=/mnt/picea/storage/reference/Populus-tremula/v1.1/gtf/Potra01-gene-mRNA-wo-intron.gtf
kallisto_fasta=/mnt/picea/storage/reference/Populus-tremula/v1.1/fasta/Potra01-mRNA.fa
kallisto_index=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/kallisto/Potra01-mRNA.fa.inx
start=2
end=9

module load bioinfo-tools FastQC Trimmomatic sortmerna star samtools htseq multiqc kallisto

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz"`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
  bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -g $genome -m 128 \
  -G $gtf -H $gff3 -t -f $kallisto_fasta -K $kallisto_index $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out
done
