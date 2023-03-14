#!/bin/bash


## this is RNA-seq PE processing pipeline for GBM collaboration
## including umi_adapter de-duplicated
## example

##install bcl2fastq.sif

#singularity exec --bind /mnt/nas/yh/10.RNA_seq_poly_A/PQ111A06/1./:/data1/yh/ bcl2fastq.sif bash

#bcl2fastq --output-dir ./fastq/  --minimum-trimmed-read-length 8 --mask-short-adapter-reads 8

##cat 

    #for f in $(ls -h ./ | tr "_" "\t"| cut -f 1 | uniq | grep "RQG")

        #do

        #cat ../0.origin/${f}*R1_001.fastq.gz > ../1.fastq/${f}_R1.fq.gz
        #cat ../0.origin/${f}*R3_001.fastq.gz > ../1.fastq/${f}_R2.fq.gz
        #cat ../0.origin/${f}*R1_001.fastq.gz > ../1.fastq/${f}_R1.fq.gz

    #done &

export R1=${1}
export R2=${2}
export UMI=${3}
export FASTPPEADAPTERSFASTA=./pipeline/fastp.pe.adapters.fasta
export HISAT2REF=./pipeline/ref/hisat2/alignref/align/chm13v2.0_ERCC
export REFREF=./pipeline/ref/picard/GCF_009914755.1_T2T-CHM13v2.0_genomic_ERCC92_refFlat_2.txt
export REFSEQ=./pipeline/ref/hisat2/alignref/chm13v2.0_ERCC.fa
export HTSEQGTF=./pipeline/ref/htseq/GCF_009914755.1_T2T-CHM13v2.0_genomic_ERCC92.gtf
export CPU="12"

# check for empty parameters!
[ -z "${R1}" ] && { echo "R1 is empty" ; exit 1; }
[ -z "${R2}" ] && { echo "R2 is empty" ; exit 1; }
[ -z "${UMI}" ] && { echo "UMI is empty" ; exit 1; }

# owner and group access; block all others
umask 0007

echo $(date) RNA PE STARTED

# report parameters
echo R1=${R1}
echo R2=${R2}
echo UMI=${UMI}
echo HISAT2REF=${HISAT2REF}
echo REFREF=${REFREF}
echo REFSEQ=${REFSEQ}
echo HTSEQGTF=${HTSEQGTF}

# check parameters
[ ! -f "${R1}" ] && { echo "${R1} does not exist" ; exit 2 ; }
[ ! -f "${R2}" ] && { echo "${R2} does not exist" ; exit 2 ; }
[ ! -f "${UMI}" ] && { echo "${UMI} does not exist" ; exit 2 ; }
[ ! -f "${HISAT2REF}.1.ht2" ] && { echo "${HISAT2REF} does not exist" ; exit 2 ; }
[ ! -f "${REFREF}" ] && { echo "${REFREF} does not exist" ; exit 2 ; }
[ ! -f "${REFSEQ}" ] && { echo "${REFSEQ} does not exist" ; exit 2 ; }
[ ! -f "${HTSEQGTF}" ] && { echo "${HTSEQGTF} does not exist" ; exit 2 ; }

# create output directory
export OUTPUTDIR="${R1}-output"
[ ! -d "${OUTPUTDIR}" ] && mkdir -p "${OUTPUTDIR}"
export OUTPREFIXR1="${OUTPUTDIR}/$(basename ${R1})"
export OUTPREFIXR2="${OUTPUTDIR}/$(basename ${R2})"

# set up the command
export FASTQCCMD="./pipeline/fastqc_v0.11.9.sif fastqc"
export FASTPCMD="./pipeline/fastp_v0.23.2.sif fastp"
export HISAT2CMD="./pipeline/hisat2.v2.1.0_2_deb_cv1.sif hisat2"
export UMITOOLCMD="./pipeline/umi_tools_v1.0.1.sif umi_tools "
export PICARDCMD="./pipeline/picard-slim_v2.27.4.sif picard" 
export STDJAVAMEM="-Xms4g -Xmx4g"
export HIGHJAVAMEM="-Xms4g -Xmx8g"
export SVRJAVAOPTS="-XX:ParallelGCThreads=2 -XX:ConcGCThreads=2"
export SAMTOOLSCMD="./pipeline/samtools_v1.15.1.sif samtools"
export HTSEQ="./pipeline/htseq.v0.11.2-1_deb_py3_cv1.sif htseq-count"


# add UMI_hadder
umi_tools extract \
--bc-pattern=NNNNNNNNN \
--stdin=${UMI} \
--read2-in=${R1} \
--stdout=${OUTPREFIXR1}_UMI_R1.gq.gz \
--read2-stdout


#umi_tools extract \
--bc-pattern=NNNNNNNNN \
--stdin=${UMI} \
--read2-in=${R2} \
--stdout=${OUTPREFIXR1}_UMI_R2.gq.gz \
--read2-stdout
echo $(date) umi_tools extract COMPLETED


# fastqc 
 ${FASTQCCMD} \
 -o ${OUTPUTDIR} \
 --noextract \
 --nogroup \
 ${OUTPREFIXR1}_UMI_R1.gq.gz &

 ${FASTQCCMD} \
 -o ${OUTPUTDIR} \
 --noextract \
 --nogroup \
 ${OUTPREFIXR1}_UMI_R2.gq.gz &
 echo $(date) FASTQC COMPLETED

#trim
 fastp \
 -w 4 \
 -i ${OUTPREFIXR1}_UMI_R1.gq.gz \
 -I ${OUTPREFIXR1}_UMI_R2.gq.gz \
 -o ${OUTPREFIXR1}_UMI_R1_trimmed.fq.gz \
 -O ${OUTPREFIXR2}_UMI_R2_trimmed.fq.gz \
 --unpaired1 ${OUTPREFIXR1}_UMI_read2_failed.gz  \
 --unpaired2 ${OUTPREFIXR1}_UMI_read1_failed.gz \
 --failed_out ${OUTPREFIXR1}_UMI_failed.gz \
 --adapter_fasta ${FASTPPEADAPTERSFASTA} \
 -j ${OUTPREFIXR1}_UMI_fastp.gz.json \
 -h ${OUTPREFIXR1}_UMI_fastp.gz.html \
  -r -W 4 -M 20 -l 36 -5 \
 -5 --cut_front_window_size 1 \
 --cut_front_mean_quality 3 \
 -3 --cut_tail_window_size 1 \
 --cut_tail_mean_quality 3 \
  2>fastp.log
 echo $(date) fastp COMPLETED

#fastqc 
 ${FASTQCCMD} \
 -o ${OUTPUTDIR} \
 --noextract \
 --nogroup \
 ${OUTPREFIXR1}_UMI_R1_trimmed.fq.gz &

 ${FASTQCCMD} \
 -o ${OUTPUTDIR} \
 --noextract \
 --nogroup \
 ${OUTPREFIXR1}_UMI_R2_trimmed.fq.gz &
 echo $(date) FASTQC COMPLETED


#align
 ${HISAT2CMD} --dta --seed 1680 --rna-strandness FR \
    -p ${CPU} -x ${HISAT2REF} \
    -1 ${OUTPREFIXR1}_UMI_R1_trimmed.fq.gz  \
    -2 ${OUTPREFIXR2}_UMI_R2_trimmed.fq.gz \
    2>${OUTPREFIXR2}_align.log \
    | ${SAMTOOLSCMD} view -Sbu - \
    | ${SAMTOOLSCMD} sort -@ ${CPU} -o ${OUTPREFIXR1}.bam
 echo $(date) HISAT2 COMPLETED

 ${PICARDCMD} ${STDJAVAMEM} ${SVRJAVAOPTS} CollectRnaSeqMetrics \
 --REF_FLAT ${REFREF} \
 --INPUT ${OUTPREFIXR1}.bam \
 --OUTPUT ${OUTPREFIXR1}.bam_UMI_rna_metrics \
 --STRAND SECOND_READ_TRANSCRIPTION_STRAND \
 echo $(date) CollectRnaseqMetrics COMPLETED

 ${PICARDCMD} ${STDJAVAMEM} ${SVRJAVAOPTS} CollectAlignmentSummaryMetrics \
 --REFERENCE_SEQUENCE ${REFSEQ} \
 --INPUT ${OUTPREFIXR1}.bam \
 --OUTPUT ${OUTPREFIXR1}.bam_UMI_alignment_umi_metrics
 echo $(date) CollectAlignmentSummaryMetrics COMPLETED

 AnnotateBamWithUmis
 fgbio AnnotateBamWithUmis \
 -i RQG11108A01_prealin.bam \
 -f RQG11108A01_UMI.fq -o  \
 RQG11108A01_anno.bam \
 -t RX -r +M  -Xms4g -Xmx8g 

#Read grouping
 umi_tools group \
 -I ${OUTPREFIXR1}_align.bam \
 --paired --group-out=groups.tsv \
 --output-bam \
 -S ${OUTPREFIXR1}_group.bam \
 -L ${OUTPREFIXR1}_group.log

 samtools sort \
 -@ 12 -m 3gb \
 ${OUTPREFIXR1}_group.bam  > RQG11108A01_sort.bam

#index
samtools index \
-@ 12 ${OUTPREFIXR1}.bam

#Deduplication
${UMITOOLCMD}  dedup \
--read-length 71  \
-I ${OUTPREFIXR1}.bam \
--paired \
--output-stats=${OUTPREFIXR1}_deduplicated-stats.txt \
-S ${OUTPREFIXR1}_UMI_dedup.bam -L ${OUTPREFIXR1}_UMI_dedup.log
    


${SAMTOOLSCMD} view -@ ${CPU} \
-F 524 -f 2 \
${OUTPREFIXR1}_UMI_dedup.bam \
| ${SAMTOOLSCMD} sort -@ ${CPU} -u -o ${OUTPREFIXR1}_UMI_filter.bam
echo $(date) samtools de-duplication COMPLETED


# count
${HTSEQ} -f bam --stranded=yes -a 30 \
--secondary-alignments=ignore \
--supplementary-alignments=ignore \
--additional-attr=gene_name \
${OUTPREFIXR1}_UMI_filter.bam \
${HTSEQGTF}  > ${OUTPREFIXR1}_UMI_count.txt
echo $(date) htseq-count COMPLETED



#
echo $(date) RNA-seq with PE ENDED





 #umi_tools docker : https://hub.docker.com/r/fredhutch/umi_tools/tags
    #singularity build  umi_tools_v1.0.1.sif docker://fredhutch/umi_tools:1.0.1

