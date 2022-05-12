#!/bin/bash
#SBATCH --job-name=counts.htseq.all
#SBATCH --mail-user=lihe.liu@wisc.edu
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2990mb
#SBATCH --time=72:00:00

module load htseq
module load samtools

##
#cd /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day14/ 
#samtools view -Su $1.sam | samtools sort -o $1".sorted.bam"
#samtools index $1".sorted.bam"

python -m HTSeq.scripts.count --format=bam --order=name --stranded=no --minaqual=0 --type=exon --idattr=gene_id --mode=intersection-nonempty --quiet /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day14/$1.sorted.bam /blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/Bos_taurus.ARS-UCD1.2.101.clean.gtf >> /blue/mateescu/lihe.liu/CoolHeat0708/htseqCount/Day14/$1'.counts.raw.txt'

##
#cd /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day42/       
#samtools view -Su $1.sam | samtools sort -o $1".sorted.bam"
#samtools index $1".sorted.bam"

python -m HTSeq.scripts.count --format=bam --order=name --stranded=no --minaqual=0 --type=exon --idattr=gene_id --mode=intersection-nonempty --quiet /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day42/$1.sorted.bam /blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/Bos_taurus.ARS-UCD1.2.101.clean.gtf >> /blue/mateescu/lihe.liu/CoolHeat0708/htseqCount/Day42/$1'.counts.raw.txt'

##
#cd /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day84/
#samtools view -Su $1.sam | samtools sort -o $1".sorted.bam"
#samtools index $1".sorted.bam"

python -m HTSeq.scripts.count --format=bam --order=name --stranded=no --minaqual=0 --type=exon --idattr=gene_id --mode=intersection-nonempty --quiet /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day84/$1.sorted.bam /blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/annotation/Bos_taurus.ARS-UCD1.2.101.clean.gtf >> /blue/mateescu/lihe.liu/CoolHeat0708/htseqCount/Day84/$1'.counts.raw.txt'








