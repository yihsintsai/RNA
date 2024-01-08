#!/bin/bash

ucscref="/mnt/nas2/yh/10.RNA_seq_poly_A/reference/htseq/chm13.ucsc_v2.0.gene_annotation.sorted.filtered.filter_ky_ERCC92.gff3"
ncbiref="/mnt/nas2/yh/10.RNA_seq_poly_A/reference/htseq/GCF_009914755.1_ncbi.T2T-CHM13v2.0_genomic_ERCC92_filter_20230314.gtf"
bam=${1}
dir=${2}
echo $ucscref
echo $ncbiref
echo $bam

echo check $dir
if [[ ! -d $dir ]] ; then
    echo $dir not exit
    mkdir $dir
 fi

#interscet
# echo `date` $bam "intersect reference"
#     if [[ -f $bam && $(echo $bam | grep -o "bam") ]] ; then
          sample=$(echo ${bam%_R1.*} | sed -e 's/.*\///g')
          echo "sample="$sample
#         for f in $ucscref $ncbiref
#             do
#                 if [[ -f $f &&  $(echo $f | grep -o "ucsc" ) ]] ; then
#                     echo $f
#                     echo ucsc
#                     name="ucsc"
#                     elif [[ -f $f &&  $(echo $f | grep -o "ncbi" ) ]] ; then
#                         echo $f
#                         echo ncbi
#                         name="ncbi"
#                             else
#                                 echo $f not exit && exit 0
#                     fi

#                 bedtools intersect -wa \
#                 -a $bam \
#                 -b ${f} \
#                 -s -f 0.1 \
#                 > $dir/${bam}_${name}.bed
#             done
#         else
#             echo $bam not exit
#         fi
#         wait
# echo `date` $bam "intersect reference completed"

 echo `date` $bam "grep ambiguous read"
    #parameter setting
        sam=$(echo ${bam%_R1.*}.sam)
        chr="chr2"
        readfile=$(echo ${sample}.amb.read)
        echo $sam $chr $readfile

    #make  $dir/$readfile
         cat $sam | awk '{if($3=="chr2") {print $0}}' | grep "ambiguous" >> $dir/$readfile
         wait
    #make splite bed file
        if [[ -f $dir/$readfile ]] ; then
        #check parameter
            file=$(echo $readfile)
            thread_num="12"
            echo $file $thread_num

        echo splite file
            split -d -l 1000 $dir/${file} $dir/${sample}_split

        echo check splite files
            readrow=$(wc -l $file)
            filenum=$(ls $dir/${sample}_split* | wc -l)
            if [[ $filenum==$(echo $($readrow/1000+1)) ]]; then
                echo splite file number correct
                list=$(echo $(ls $dir/${sample}_split*))
                else
                echo "splite wrong" && exit 0
                fi


        echo create pipe
            [ ! -p tmp ] && mkfifo tmp
            exec 9<>tmp
        # fill in the pipe
            for ((i=0;i<thread_num;i++)); do
                echo >&9
                done

        echo mkde splite bed file
            for s in $list;
                do
                {
                # remove one line per thread
                    read -u 9
                    {
                        IFS=$'\n'
                        for line in $(cat < $s)
                            do
                                if [[ $(echo $line | grep "ambiguous") ]]; then
                                    position=$(echo $line | cut -f 3,4,8)
                                    gene_name=$(echo $line | cut -f 1)
                                    score=$(echo $line | cut -f 5)
                                    strand=$(echo $line | cut -f 9)
                                    #check for/rev
                                        if [[ $(echo $strand | sed "s/[0-9]//g") ]]; then
                                        # reverse
                                        gene_id=$(echo $gene_name/2)
                                        strand="-"
                                        chr=$(echo $position | cut -f 1)
                                        str=$(echo $position | cut -f 2)
                                        ends=$(echo $position | cut -f 3)
                                        else
                                        # forward
                                            gene_id=$(echo $gene_name/1)
                                            strand="+"
                                            chr=$(echo $position | cut -f 1)
                                            str=$(echo $position | cut -f 3)
                                            ends=$(echo $position | cut -f 2)
                                        fi
                                        echo -e $chr'\t'$str'\t'$ends'\t'$gene_id'\t'$score'\t'$strand >> $dir/${s}.bed
                                        else
                                            echo $line
                                            echo $line error >> chm13_test_err && exit 0
                                    fi
                            done
                            echo >&9
                    } &
                }
                done
                wait

                echo  close pipe
                    exec 9>&-
                    rm tmp
            else
                echo $dir/$readfile "not exit"
             fi
    #cat splite file
        echo `date` "read pooling"
            # check splite bed file
            n=$(ls $dir/${sample}_split*bed|wc -l)
            if [[ $n==$filenum ]] ; then
                for f in $(ls $dir/${sample}_split*bed)
                        do
                        cat ${f} >> $dir/${sample}.amb.chr2.bed
                done &
                else
                echo "splite file number wrong" && exit 0
                fi

        echo `date` "read pooling completed"



echo `date` "make bed from bam"
    cat $sam | grep "ambiguous" | cut -f 1| sort |uniq > $dir/${sample}_ambiguous.name
    if [[ -f $dir/${sample}_ambiguous.name ]] ;then

        while read i ;
            do
            cat $dir/${sample}_R1.fq.gz.HQ_filter.bam_ucsc.bed | grep $i >> $dir/${sample}.ucsc.chr2.bed
            cat $dir/${sample}_R1.fq.gz.HQ_filter.bam_ncbi.bed | grep $i >> $dir/${sample}.ncbi.chr2.bed
            wait
            done < $dir/${sample}_ambiguous.name
        else
            echo $dir/${sample}_ambiguous.name "not exit" && exit 0
        fi
echo `date` "make bed from bam done"
