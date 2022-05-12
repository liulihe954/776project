#!/bin/bash
#SBATCH --job-name=hisat.all
#SBATCH --mail-user=lliu356@wisc.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3990mb
#SBATCH --time=240:00:00

module load hisat2 

hisat2 -p 8 --dta-cufflinks --novel-splicesite-outfile /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day14/splicesite/$1'.bed' -x /blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/index/index-hisat2-clean/Bos_taurus.ARS-UCD1.2 -1 /blue/mateescu/lihe.liu/CoolHeat0708/editData/Day14/$1/$1'-14_R1_val_1.fq' -2 /blue/mateescu/lihe.liu/CoolHeat0708/editData/Day14/$1/$1'-14_R2_val_2.fq' -S /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day14/$1'.sam' &> /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day14/$1'.hisat2.log'


hisat2 -p 8 --dta-cufflinks --novel-splicesite-outfile /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day42/splicesite/$1'.bed' -x /blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/index/index-hisat2-clean/Bos_taurus.ARS-UCD1.2 -1 /blue/mateescu/lihe.liu/CoolHeat0708/editData/Day42/$1/$1'-42_R1_val_1.fq' -2 /blue/mateescu/lihe.liu/CoolHeat0708/editData/Day42/$1/$1'-42_R2_val_2.fq' -S /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day42/$1'.sam' &> /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day42/$1'.hisat2.log'


hisat2 -p 8 --dta-cufflinks --novel-splicesite-outfile /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day84/splicesite/$1'.bed' -x /blue/mateescu/lihe.liu/AltSplicing/ARS-UCD1.2/index/index-hisat2-clean/Bos_taurus.ARS-UCD1.2 -1 /blue/mateescu/lihe.liu/CoolHeat0708/editData/Day84/$1/$1'-84_R1_val_1.fq' -2 /blue/mateescu/lihe.liu/CoolHeat0708/editData/Day84/$1/$1'-84_R2_val_2.fq' -S /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day84/$1'.sam' &> /blue/mateescu/lihe.liu/CoolHeat0708/hisat2Align/Day84/$1'.hisat2.log'





