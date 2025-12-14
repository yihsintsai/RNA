#!/bin/bash
set -euo pipefail

# ------------------------ Usage & Args ------------------------
usage() {
  cat <<'USAGE'
bash sbatch_cellranger.sh -p PROJECT_ID -g ORGANISM -s SC_CLASS -c CUSTOM_REF

Required:
  -p	PROJECT_ID		e.g. HS...
  -s    SC_CLASS		e.g. NORMAL / OCM / PLEX ..
Optional:
  -g    ORGANISM                species（human or mouse...)
  -c    CUSTOM_REF              path for CUSTOM_REF
  -w    WORK_PATH		path to run PROJECT_ID (default:/mnt/NFS/leo/yihsin/scRNA)
  -h				show this help
USAGE
}

PROJECT_ID=""
SC_CLASS=""
ORGANISM=""
CUSTOM_REF=""
WORK_PATH="/mnt/scRNA"

while getopts ":p:s:g:c:w:h" opt; do
    case "$opt" in
          p) PROJECT_ID="$OPTARG" ;;
          s) SC_CLASS="$OPTARG" ;;
          g) ORGANISM="$OPTARG" ;;
          c) CUSTOM_REF="$OPTARG" ;;
          w) WORK_PATH="$OPTARG" ;;
          h) usage; exit 0 ;;
          :)  echo "Option -$OPTARG requires an argument." ; usage; exit 1 ;;
          \?) echo "Unknown option: -$OPTARG"              ; usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

# --------- check parameter ------------------------------------------------
echo "CHECK PARAMETER"
if [[ -z ${PROJECT_ID} || -z ${SC_CLASS} ]] ; then
        echo "❗ Missing : Check -s SC_CLASS and -p PROJECT_ID" 
        exit 1
fi

if [[ -n "${CUSTOM_REF}" && -d "${CUSTOM_REF}" ]]; then
        ref_path=${CUSTOM_REF}
elif [[ "${ORGANISM}" == "human" ]]; then
        ref_path=""
	      ref_probe=""
elif [[ "${ORGANISM}" == "mouse" ]]; then
        ref_path=""
	      ref_probe=""
else
        echo "❗ --- ORGANISM must be 'human' or 'mouse', or CUSTOM_REF must exist ---" >&2
        exit 1
fi

## ------- CHECK PROJECT_DIR ----------------------------------------------
echo "CHECK PROJECT DIR"
PROJECT_DIR="${WORK_PATH}/${PROJECT_ID}"
if [[ -d ${PROJECT_DIR} ]]; then
  	  echo "PROJECT PATH : ${PROJECT_DIR}"
elif [[ -d "/mnt/NFS/pisces/yihsin/scRNA/${PROJECT_ID}" ]]; then
  	  echo "PROJECT PATH : /mnt/${PROJECT_ID}"	
  	  PROJECT_DIR="/mnt/${PROJECT_ID}"
else
  	  echo "Create project dir: ${PROJECT_DIR}"
  	  mkdir -p ${PROJECT_DIR}
fi

if [[ -d "${PROJECT_DIR}/data" ]]; then
  	  echo "PROJECT data path already exists: ${PROJECT_DIR}/data"
else
  	  echo "Create project data dir: ${PROJECT_DIR}/data"
  	  mkdir -p "${PROJECT_DIR}/data"
fi


## -------- buld soft file in data file -----------------------------------
echo "BUILD SOFT LINK IN PROJECT DATA FILE"
cd ${PROJECT_DIR}/data || {
	    echo "FAIL TO cd TO ${PROJECT_DIR}/data" >&2
      exit 1
}

shopt -s nullglob
check_fastq=( *.fastq.gz)
if ((${#check_fastq[@]} > 0)) ; then
	    if [[ "${check_fastq[*]}" =~ _S[0-9]_L[0-9]{3}_R[0-9]_[0-9]{3}\.fastq\.gz ]]; then
		      echo "FASTQ FILE already exists with correct format: ${check_fastq[*]}, PLEASE CHECK"
	    else
		      echo "RENAME DATA FORMATE"
		      bash /mnt/replace_name.sh
	    fi
else
	    src_fastq=(/mnt/NFS/pisces/projects/${PROJECT_ID}/FASTQ/concat/rawData/*gz)
	    if ((${#src_fastq[@]} > 0)) ; then
		        ln -s /mnt/${PROJECT_ID}/rawData/*gz ./
		        echo "RENAME DATA FORMATE"
		        bash /mnt/replace_name.sh
	    else 
		        echo "DIR NOT EXIST : /mnt/${PROJECT_ID}/rawData/ , PLEASE CHECK"
		        exit 1
      fi
fi
shopt -u nullglob

## -------- check single cell class ---------------------------
### cellranger count
if [[ ${SC_CLASS} == "NORMAL" ]] ; then
	echo "START CELLRANGER COUNT ANALYSIS"
	sample_id=$(ls *fastq.gz | cut -f 1 -d "_" | sort |uniq)
	for s in ${sample_id} 
		do
			echo "sbatch ${PROJECT_ID} : ${s} file"
				
			sed -e "s/SC_PROJECT/${PROJECT_ID}/g" \
			    -e "s|REF_PATH|${ref_path}|g" \
			    -e "s/SAMPLE_ID/${s}/g" \
			    cellranger_gex.sh \
			    > ${PROJECT_DIR}/cellranger_gex_${s}.sh

			sbatch ${PROJECT_DIR}/cellranger_gex_${s}.sh 
	done

### pip-seq
elif [[ ${SC_CLASS} == "PIP" ]]; then
	echo "START PIP-SEQ ANALYSIS"
	sbatch pip-seq.sh ${PROJECT_ID} ${ORGANISM}


### cellranger multi-PLEX
elif [[ ${SC_CLASS} == "PLEX" ]]; then

	echo "CHECK GEM PLEX INPUT FILE"
        shopt -s nullglob
        src_fastq=( *R1_001.fastq.gz)
        for fq in "${src_fastq[@]}"
                do
                        f=${fq%%_*}
                        echo ${f}               
                        GEX_FLEX_FILE=${PROJECT_DIR}/GEX_PLEX_${f}.csv
                        sed -e "s/SC_PROJECT/${f}/g" \
                            -e "s|INPUT_FILE|${GEX_FLEX_FILE}|g" \
                            cellranger_plex.sh \
                            > ${PROJECT_DIR}/cellranger_plex_${f}.sh
                        sed -e "s|REF_PATH|${ref_path}|g" \
                            -e "s|REF_PROBE|${ref_probe}|g" \
                            -e "s/SC_PROJECT/${f}/g" \
                            -e "s/PROJECT_ID/${PROJECT_ID}/g" \
                           GEX_PLEX_file.csv \
                            > ${PROJECT_DIR}/GEX_PLEX_${f}.csv

                echo "add GEX PLEX sample name and barcode id to ${PROJECT_DIR}/GEX_PLEX_${f}.csv"
        done
        shopt -u nullglob

### cellranger multi-OCM
elif [[ ${SC_CLASS} == "OCM" ]]; then
	echo "CHECK GEM-X OCM INPUT FILE"
	shopt -s nullglob
	src_fastq=( *R1_001.fastq.gz)
	for fq in "${src_fastq[@]}" 
		do
			f=${fq%%_*}
			echo ${f}		
			GEX_OCM_FILE=${PROJECT_DIR}/GEX_OCM_${f}.csv
			sed -e "s/SC_PROJECT/${f}/g" \
			    -e "s|INPUT_FILE|${GEX_OCM_FILE}|g" \
			    cellranger_flex.sh \
		    	    > ${PROJECT_DIR}/cellranger_flex_${f}.sh
			sed -e "s|REF_PATH|${ref_path}|g" \
			    -e "s|REF_PROBE|${ref_probe}|g" \
			    -e "s/SC_PROJECT/${f}/g" \
			    -e "s/PROJECT_ID/${PROJECT_ID}/g" \
			    GEX_OCM_file.csv \
			    > ${PROJECT_DIR}/GEX_OCM_${f}.csv

		echo "add GEX OCM sample name and barcode id to ${PROJECT_DIR}/GEX_OCM_${f}.csv"
	done
	shopt -u nullglob
else	
	echo "ERROR INPUT FOR SC_CLASS PARAMETER : ${SC_CLASS}, PLEASE CHECK"

fi
		
