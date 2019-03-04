#!/bin/bash

# job name
#
#PBS -N assemblyFromContigsAndNanopores

# when and who to mail
#PBS -m bea
#PBS -M maja.fabijanic@gmail.com

# queue:
#
#PBS -q MASTER
# request resources
#
#PBS -l select=1:ncpus=12:mem=400G

cd /common/WORK/mfabijanic/Sponges/ESU_v0.03

INFILE=/common/WORK/mfabijanic/Sponges/ESU_v0.03/redundans/scaffolds.reduced.fa
OUTDIRPREFIX=pilonMacOnredundans
OUTFILE=pilonMacrogenOnredundansScafolds
FRAGSBAM=MacrogenOnredundansbwa.bam

java -Xmx350G -jar /common/WORK/mfabijanic/programs/pilon-1.22.jar --genome $INFILE \
--frags $FRAGSBAM \
--output $OUTFILE \
--changes --vcf --tracks \
--threads 10 --verbose \
--outdir $OUTDIRPREFIX
