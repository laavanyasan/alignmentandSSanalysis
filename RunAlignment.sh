#!/bin/bash
#
#SBATCH --mem=250000
#SBATCH -J align
#SBATCH -o Runalign.%J.output
#SBATCH -e Runalign.%J.error
#SBATCH --array=1-3
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=ls331@duke.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -p new,all

#LOAD SOFTWARE
export PATH="/data/reddylab/gjohnson/preseq/:$PATH"

module load bowtie2
module load samtools
module load bedtools2
module load R/3.2.0-gcb01
module load gcc

source /data/reddylab/software/miniconda2/bin/activate alex


#SET VARIABLES
R1=$(head -${SLURM_ARRAY_TASK_ID} /data/reddylab/Laavanya/pcosStarrSeq/experiments/20210410-DENND1APOPSS-OPlib-v2/inputs/inputfile | tail -1 | awk '{print$1}')
R2=$(head -${SLURM_ARRAY_TASK_ID} /data/reddylab/Laavanya/pcosStarrSeq/experiments/20210410-DENND1APOPSS-OPlib-v2/inputs/inputfile | tail -1 | awk '{print$2}')
SAMPLE=$(head -${SLURM_ARRAY_TASK_ID} /data/reddylab/Laavanya/pcosStarrSeq/experiments/20210410-DENND1APOPSS-OPlib-v2/inputs/inputfile  | tail -1 | awk '{print$3}')

echo ${R1}
echo ${R2}
echo ${SAMPLE}

#cd /data/reddylab/Laavanya/pcosStarrSeq/experiments/20191126-SS_H295R_v1/
#RUN FASTQC BEFORE TRIMMING
/data/reddylab/software/FastQC/fastqc ${R1} ${R2} \
-t 31 \
-o /data/reddylab/Laavanya/pcosStarrSeq/experiments/20210410-DENND1APOPSS-OPlib-v2/fastqc

#MAKE SAMPLE DIRECTORY
mkdir alignments/${SAMPLE}


#ALIGN SAMPLES
bowtie2 \
-x /data/reddylab/gjohnson/reference_data/hg38 \
-p 31 \
-X 2000 \
--sensitive \
--un alignments/${SAMPLE}/${SAMPLE}.unmapped.sam \
-1 ${R1} \
-2 ${R2} \
| samtools view -@ 31 -Shu -L /data/reddylab/gjohnson/reference_data/hg38.genome.bed -f 3 -q 10 - \
| samtools sort -@ 31 - -o alignments/${SAMPLE}/${SAMPLE}.f3q10.bam


#INDEX BAM FILE
samtools index alignments/${SAMPLE}/${SAMPLE}.f3q10.bam

#CREATE BEDPE FILE FOR SAMPLE
samtools view -@ 31 -Shu alignments/${SAMPLE}/${SAMPLE}.f3q10.bam \
| samtools sort -@ 31 -n - \
| bedtools bamtobed -i stdin -bedpe \
| cut -f 1,2,6 \
| sortBed > alignments/${SAMPLE}/${SAMPLE}.bedpe


#MAKE NORMALIZED BIGWIGS
bamCoverage \
--bam alignments/${SAMPLE}/${SAMPLE}.f3q10.bam \
--binSize 1 \
--normalizeUsingRPKM \
--numberOfProcessors max \
--extendReads \
--outFileName bigwigs/${SAMPLE}.bw


#GET FASTQ READ COUNT FOR EACH SAMPLE
less ${R1} \
| wc -l \
| awk -v var="${SAMPLE}" '{print $0/4, var}' >> counts/fastq.counts


#GET f3q10 FRAGMENT COUNT FOR EACH SAMPLE
samtools view -@ 31 alignments/${SAMPLE}/${SAMPLE}.f3q10.bam \
| wc -l \
| awk -v var="${SAMPLE}" '{print $1/2, var}' >> counts/f3q10.frag.counts


#GET f3q10 RMDUP FRAGMENT COUNT FOR EACH SAMPLE
samtools rmdup alignments/${SAMPLE}/${SAMPLE}.f3q10.bam alignments/${SAMPLE}/dups.bam

samtools view alignments/${SAMPLE}/dups.bam \
| wc -l \
| awk -v var="${SAMPLE}" '{print $1/2, var}' >> counts/rmdup.f3q10.frag.counts

rm alignments/${SAMPLE}/dups.bam


#GET FRAGMENT LENGTHS FOR THE 1ST 0.5M FRAGMENTS FOR EACH SAMPLE
samtools view -q 10 -f 3 -F 284 alignments/${SAMPLE}/${SAMPLE}.f3q10.bam \
| head -500000 \
| awk '{print $9}'  > lengths/${SAMPLE}.lengths


#GET PRESEQ ESTIMATES FOR EACH SAMPLE
preseq c_curve -v -P -B alignments/${SAMPLE}/${SAMPLE}.f3q10.bam -o preseq.complexity/${SAMPLE}.ccurve.txt 2> preseq.complexity/${SAMPLE}.ccurve.verbose.txt

preseq lc_extrap -v -P -B alignments/${SAMPLE}/${SAMPLE}.f3q10.bam -o preseq.complexity/${SAMPLE}.lc.extrap.txt  2> preseq.complexity/${SAMPLE}.lc.extrap.verbose.txt

less preseq.complexity/${SAMPLE}.lc.extrap.txt \
| cut -f 2 \
| tail -n +2 \
| sort -k1,1nr \
| head -1 \
| awk -v OFS="\t" -v var="${SAMPLE}" '{print $1, var}'  >> preseq.maxs

