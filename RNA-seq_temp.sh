	#!/bin/bash

##
## this is RNA-seq PE processing pipeline for CGMH GBM collaboration
##
## example


##   export SINGULARITY_BINDPATH=/mnt/nas/GBM
##   cd /mnt/nas/GBM
##   ./pipeline/wgs.pe.sh testdata/testdata_S1_L001_R1_001.fastq.gz testdata/testdata_S1_L001_R2_001.fastq.gz 2>&1 | tee testdata/testdata_S1_L001_R1_001.fastq.gz.run.log
##


export R1=${1}
export R2=${2}
export HISAT2REF=./pipeline/ref/hisat2/alignref/align/chm13v2.0_ERCC
export REFSEQ=./pipeline/ref/hisat2/alignref/chm13v2.0_ERCC.fa
export FASTPPEADAPTERSFASTA=./pipeline/fastp.pe.adapters.fasta
export REFREF=./pipeline/ref/picard/GCF_009914755.1_T2T-CHM13v2.0_genomic_ERCC92_refFlat_2.txt
export HTSEQGTF=./pipeline/ref/htseq/GCF_009914755.1_T2T-CHM13v2.0_genomic_ERCC92.gtf
export CPU="12"

# check for empty parameters!
[ -z "${R1}" ] && { echo "R1 is empty" ; exit 1; }
[ -z "${R2}" ] && { echo "R2 is empty" ; exit 1; }

# owner and group access; block all others
umask 0007

echo $(date) RNA PE STARTED
# report parameters
echo R1=${R1}
echo R2=${R2}
echo HISAT2REF=${HISAT2REF}
echo REFSEQ=${REFSEQ}
echo HTSEQGTF=${HTSEQGTF}
echo REFREF=${REFREF}
# check parameters
[ ! -f "${R1}" ] && { echo "${R1} does not exist" ; exit 2 ; }
[ ! -f "${R2}" ] && { echo "${R2} does not exist" ; exit 2 ; }
[ ! -f "${HISAT2REF}.1.ht2" ] && { echo "${HISAT2REF} does not exist" ; exit 2 ; }
[ ! -f "${REFSEQ}" ] && { echo "${REFSEQ} does not exist" ; exit 2 ; }
[ ! -f "${HTSEQGTF}" ] && { echo "${HTSEQGTF} does not exist" ; exit 2 ; }
[ ! -f "${REFREF}" ] && { echo "${REFREF} does not exist" ; exit 2 ; }

# create output directory
export OUTPUTDIR="${R1}-output"
[ ! -d "${OUTPUTDIR}" ] && mkdir -p "${OUTPUTDIR}"
export OUTPREFIXR1="${OUTPUTDIR}/$(basename ${R1})"
export OUTPREFIXR2="${OUTPUTDIR}/$(basename ${R2})"

# set up the command
export FASTQCCMD="./pipeline/fastqc_v0.11.9.sif fastqc"
export FASTPCMD="./pipeline/fastp_v0.23.2.sif fastp"
export HISAT2="/mnt/nas/GBM/pipeline/hisat2.v2.1.0_2_deb_cv1.sif hisat2"
export HTSEQ="/mnt/nas/GBM/pipeline/htseq.v0.11.2-1_deb_py3_cv1.sif htseq-count"
export SAMTOOLSCMD="./pipeline/samtools_v1.15.1.sif samtools"
export PICARDCMD="./pipeline/picard-slim_v2.27.4.sif picard" 
export STDJAVAMEM="-Xms4g -Xmx4g"
export HIGHJAVAMEM="-Xms4g -Xmx8g"
export SVRJAVAOPTS="-XX:ParallelGCThreads=2 -XX:ConcGCThreads=2"
[ ! -f "$(echo $FASTQCCMD | cut -f 1 -d\  )" ] && { echo "$(echo $FASTQCCMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $FASTPCMD | cut -f 1 -d\  )" ] && { echo "$(echo $FASTPCMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $HISAT2 | cut -f 1 -d\  )" ] && { echo "$(echo $HISAT2 | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $HTSEQ | cut -f 1 -d\  )" ] && { echo "$(echo $HTSEQ | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $SAMTOOLSCMD | cut -f 1 -d\  )" ] && { echo "$(echo $SAMTOOLSCMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $PICARDCMD | cut -f 1 -d\  )" ] && { echo "$(echo $PICARDCMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }


# FASTQC of original reads
${FASTQCCMD} -o  ${OUTPUTDIR} --noextract --nogroup ${R1} &
${FASTQCCMD} -o  ${OUTPUTDIR} --noextract --nogroup ${R2} &
wait
echo $(date) FastQC COMPLETED


# adapter and quality trimming
${FASTPCMD} -w ${CPU} \
-i ${R1} -I ${R2} \
-o ${OUTPREFIXR1}_trimmed.fq.gz -O ${OUTPREFIXR2}_trimmed.fq.gz \
--unpaired1 ${OUTPREFIXR1}_read2_failed.gz --unpaired2 ${OUTPREFIXR2}_read1_failed.gz \
--failed_out ${OUTPREFIXR1}_failed.gz \
--adapter_fasta ${FASTPPEADAPTERSFASTA} \
-j ${OUTPREFIXR1}_fastp.gz.json \
-h ${OUTPREFIXR1}_fastp.gz.html \
-r -W 4 -M 20 -l 36 \
-5 --cut_front_window_size 1 --cut_front_mean_quality 3 \
-3 --cut_tail_window_size 1 --cut_tail_mean_quality 3
echo $(date) Adapter and quality trimming COMPLETED


# FASTQC of trimmed reads
${FASTQCCMD} -o  ${OUTPUTDIR} --noextract --nogroup ${OUTPREFIXR1}_trimmed.fq.gz &
${FASTQCCMD} -o  ${OUTPUTDIR} --noextract --nogroup ${OUTPREFIXR2}_trimmed.fq.gz &
wait
echo $(date) FastQC COMPLETED


# genome mapping
${HISAT2} --dta --seed 1680 --rna-strandness FR \
	-p ${CPU} -x ${HISAT2REF} \
	-1 ${OUTPREFIXR1}_trimmed.fq.gz \
	-2 ${OUTPREFIXR2}_trimmed.fq.gz \
	2>${OUTPREFIXR2}_align.log \
	| ${SAMTOOLSCMD} view -Sbu - \
	| ${SAMTOOLSCMD} sort -@ ${CPU} -o ${OUTPREFIXR1}.bam

echo $(date) Genome mapping COMPLETED


# mapping metrics
${PICARDCMD} ${STDJAVAMEM} ${SVRJAVAOPTS} CollectRnaSeqMetrics \
--REF_FLAT ${REFREF} \
--INPUT ${OUTPREFIXR1}.bam \
--OUTPUT ${OUTPREFIXR1}.bam_rna_metrics \
--STRAND SECOND_READ_TRANSCRIPTION_STRAND 
echo $(date) CollectWgsMetrics COMPLETED

${PICARDCMD} ${STDJAVAMEM} ${SVRJAVAOPTS} CollectAlignmentSummaryMetrics \
--REFERENCE_SEQUENCE ${REFSEQ} \
--INPUT ${OUTPREFIXR1}.bam \
--OUTPUT ${OUTPREFIXR1}.bam_alignment_metrics
echo $(date) CollectAlignmentSummaryMetrics COMPLETED

# filter Q<30 ,unmapped reads
${SAMTOOLSCMD} view -@ 12 -b -q 30 \
-u ${OUTPREFIXR1}.bam \
-o ${OUTPREFIXR1}_HQ.bam
echo $(date) Extract high-quality mapped reads COMPLETED

${SAMTOOLSCMD} view -@ ${CPU} \
-F 524 -f 2 \
-u ${OUTPREFIXR1}_HQ.bam \
| ${SAMTOOLSCMD} sort -@ ${CPU} -n -o ${OUTPREFIXR1}_filter.bam
echo $(date) samtools de-duplication COMPLETED

${PICARDCMD} ${STDJAVAMEM} ${SVRJAVAOPTS} CollectAlignmentSummaryMetrics \
--REFERENCE_SEQUENCE ${REFSEQ} \
--INPUT ${OUTPREFIXR1}_filter.bam \
--OUTPUT ${OUTPREFIXR1}_filter_alignment_metrics
echo $(date) CollectAlignmentSummaryMetrics COMPLETED

# de-duplication
${PICARDCMD} ${HIGHJAVAMEM} ${SVRJAVAOPTS} MarkDuplicates \
--INPUT ${OUTPREFIXR1}_filter.bam \
--OUTPUT ${OUTPREFIXR1}_dedup.bam \
--METRICS_FILE ${OUTPREFIXR1}_dedup_metrics \
--CREATE_INDEX true \
--REMOVE_DUPLICATES true
echo $(date) MarkDuplicates COMPLETED

# count
${HTSEQ} -f bam --stranded=yes -a 30 \
--secondary-alignments=ignore \
--supplementary-alignments=ignore \
--additional-attr=gene_name \
${OUTPREFIXR1}_dedup.bam \
${HTSEQGTF}  > ${OUTPREFIXR1}_count.txt
echo $(date) htseq-count COMPLETED

# bam file metrics
${SAMTOOLSCMD} flagstat ${OUTPREFIXR1}.bam > ${OUTPREFIXR1}.bam_flagstats_metrics
echo $(date) BAM flags metrics COMPLETED
# bam file metrics
${SAMTOOLSCMD} index ${OUTPREFIXR1}.bam
echo $(date) BAM file indexing COMPLETED
# bam file metrics
${SAMTOOLSCMD} idxstats ${OUTPREFIXR1}.bam > ${OUTPREFIXR1}.bam_idxstats_metrics
#echo $(date) BAM index metrics COMPLETED


#
echo $(date) RNA-seq PE ENDED
