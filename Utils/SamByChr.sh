#!/bin/bash
#Record current date

job_name=SamByChr
input=$1
dir=$2
out=$3


function bamsort {
    echo `date` "starting job $i" >> ${job_name}.log
	echo $i
	sambamba view -F "ref_name=='${i}'" -h -f bam $input  >$dir/${out}.${i}.bam
	perl /data/rnomics5/Zhanghena/project/m5C/new_m5C_map/ref_rDNA/test/Pile_mc_lev.pl -i $dir/${out}.${i}.bam -t $out -r /data/rnomics5/Zhanghena/project/m5C/new_m5C_map/hg38_ref/hg38_primary_raw/${i}.fa --rlength 101 --minBQ 30 --overhang 6 --reads 1 --cRatio 0 --variants 1 --depth 8000 >$dir/${out}.${i}.sites
	rm $dir/${out}.${i}.bam
}
echo
		
max_jobs=$4
FIFO_FILE="./$$.fifo"
mkfifo $FIFO_FILE
exec 6<>$FIFO_FILE
rm -f $FIFO_FILE

chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrM chrX chrY)
#chr=(chr10)

for((i=1; i<=$max_jobs; i++));do
	echo
done >&6
			
for i in ${chr[@]};do
	read -u6
	{
		bamsort && {
			echo `date` "finish job $i" >> ${job_name}.log
			sleep 1
		} || {
			echo job $i error >> ${job_name}.log
		}
		echo >&6
	}& 
done
																							
# wati for all jobs to complete
wait
exec 6>&-
