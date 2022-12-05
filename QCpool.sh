#!/bin/bash

#copy data
dir="/mnt/nas/GBM"
cpdir="/mnt/nas/yh/10.RNA_seq_poly_A/R/1.counts"
for f in PQ111A06 PQ111A08 PQ111A19 PQ111A10 

do

list=$(ls -h $dir/$f/|tr "_" "\t"| cut -f 1 | uniq | grep "RQG")

	for i in $list
	do
	mkdir $cpdir/$f/5.QC/$i/
	html=$(ls -h $dir/$f/${i}_R1.fq.gz-output/ | grep "html")	
	metrics=$(ls -h $dir/$f/${i}_R1.fq.gz-output/ | grep "metrics")
	json=$(ls -h $dir/$f/${i}_R1.fq.gz-output | grep "json")
	zip=$(ls -h $dir/$f/${i}_R1.fq.gz-output | grep "zip")
	txt=$(ls -h $dir/$f/${i}_R1.fq.gz-output | grep "txt")
	count=$(ls -h $dir/$f/${i}_R1.fq.gz-output | grep "count_2")
		for n in $count
		do 
		cp $dir/$f/${i}_R1.fq.gz-output/$n $cpdir/
		#cp -r $dir/$f/${i}_R1.fq.gz-output/$n $cpdir/$f/5.QC/$i/
		done &
	done 
done 


#cat QCdata 
dir="/mnt/nas/yh/10.RNA_seq_poly_A"

#catch QC data

plID=$(ls -h $dir/ |cut -f 1 | uniq | grep "PQ")

for p in $plID 

	do

	list=$(ls -h $dir/$p/5.QC/ | cut -f 1 | uniq | grep "RQG")	

		for l in $list

			do

			wdir="/mnt/nas/yh/10.RNA_seq_poly_A/$p/5.QC/$l"
			rm $wdir/${l}_QC.txt
			rm $wdir/${l}_nonfilter.txt
			rm $wdir/${l}_MAPQ30.txt

			echo $(date) rm COMPLETE
			wait


			echo  Sample_ID"\t"${l} > $wdir/${l}_QC.txt
			

			##Total pair reads
			grep -A 2 "\"before_filtering"\" \
			$wdir/${l}_R1.fq.gz_fastp.gz.json \
			| grep "\"total_reads"\" \
			| echo $(awk '{gsub(/^\s+|,|"/, "");print}') \
			| tr ":" "\t" \
			|cut -f 2 \
			| echo  Total_pair_reads"\t"$(awk '{print $1/2}') \
			>>$wdir/${l}_QC.txt

			##Total pair reads passed q20
			grep -A 2 "\"after_filtering"\" \
			$wdir/${l}_R1.fq.gz_fastp.gz.json \
			| grep "\"total_reads"\" \
			| echo $(awk '{gsub(/^\s+|,|"/, "");print}') \
                        | tr ":" "\t" \
                        | cut -f 2 \
			| echo  Total_pair_reads_passed_q20"\t"$(awk '{print $1/2}') \
			>>$wdir/${l}_QC.txt

			##Total mapped reads
			grep "^PAIR" \
			$wdir/${l}_R1.fq.gz.bam_alignment_metrics \
			| echo  Total_mapped_reads"\t"$(awk '{print $6}') \
			>>$wdir/${l}_QC.txt

			##Total mapped reads_HQ
			grep "total" \
			$wdir/${l}_R1.fq.gz.HQ.bam_flagstats_metrics\
			| tr "+" "\t" \
			| cut -f 1 \
			echo Total_mapped_reads_HQ"\t"$(awk '{print}') \
			>>$wdir/${l}_QC.txt
			
			echo $(date) QC.txt COMPLETED


######cat htseq QC data########

			HTSEQ count origin
			echo $l > $wdir/${l}_nonfilter.txt
			echo "No quality filter" >> $wdir/${l}_nonfilter.txt
			grep "^__" \
                        $wdir/${l}_R1.fq.gz_count.txt \
                        |cut -f 3 \
                        >>$wdir/${l}_nonfilt.txt
			
			echo $(date) nonfilt.txt COMPLETED


#####cat HQ htseq QC data###########
			MAPQ >=30
			echo $l > $wdir/${l}_MAPQ30.txt
			echo "Q30(MAPQ >= 30)" >> $wdir/${l}_MAPQ30.txt
			grep "^__" \
                        $wdir/${l}_R1.fq.gz_count_2.txt \
                        |cut -f 3 \
                        >>$wdir/${l}_MAPQ30.txt
                  	
			echo $(date) MAPQ30.txt COMPLETED
			wait
			done&
		
		done 

wait


