
#! /bin/sh
###
 # @author Citu
 # @email Citu.Citu@uth.tmc.edu
 # @create date 2024-07-05 23:16:53
 # @modify date 2024-07-05 23:16:53
 # @desc Usage: sh rte_processing_pipeline.sh srr_indidual_ids.txt
###
meta_file=$1
pw="$(pwd)"
mkdir -p quantification
quantification="$(pwd)/quantification"
transcriptome="/data2/ccitu/software/cellranger-7.1.0/refdata-gex-GRCh38-2020-A"
index="/data2/ccitu/software/scTE/bin/hg38.exclusive.idx"
cell_ranger="/data2/ccitu/software/cellranger-7.1.0/bin/cellranger count"
scte="/data2/ccitu/software/scTE/bin/scTE"
reference_data="$(pwd)/Reference"
cut -f2 ${meta_file}|sort -u|while read indi; do
	echo $indi
	grep -w "$indi" $meta_file|cut -f1|sort -u|awk 'BEGIN{col=1}{print $1,col; col+=1}' OFS='\t'|while read srrd; do
		echo $srrd
		srr=$(echo "$srrd"|cut -f1)
		echo $srr
		labe=$(echo "$srrd"|cut -f2)
		echo $labe
		fasterq-dump -p -S --include-technical $srr -O $indi/
		if [[ `ls $indi/${srr}*|wc -l` == 2 ]];then
			FILE1="$indi/${srr}_1.fastq"
			echo $indi/$FILE1
			if [[ -f "$FILE1" ]]; then
    				echo "$FILE1 exists."
    				mv "$FILE1" $indi/${indi}_S1_L00${labe}_R1_001.fastq
			fi
			FILE2="$indi/${srr}_2.fastq"
			if [[ -f "$FILE2" ]]; then
    				echo "$FILE2 exists."
    				mv "$FILE2" $indi/${indi}_S1_L00${labe}_R2_001.fastq
			fi
		fi

		if [[ `ls $indi/${srr}*|wc -l` == 3 ]];then
			FILE1="$indi/${srr}_1.fastq"
			echo $indi/$FILE1
			if [[ -f "$FILE1" ]]; then
				 echo "$FILE1 exists." 
				mv "$FILE1" $indi/${indi}_S1_L00${labe}_I1_001.fastq
			fi
			FILE2="$indi/${srr}_2.fastq"
			if [[ -f "$FILE2" ]]; then
    				echo "$FILE2 exists."
    				mv "$FILE2" $indi/${indi}_S1_L00${labe}_R1_001.fastq
			fi
			FILE3="$indi/${srr}_3.fastq"
			if [[ -f "$FILE3" ]]; then
    				echo "$FILE3 exists."
   				 mv "$FILE3" $indi/${indi}_S1_L00${labe}_R2_001.fastq
			fi
		fi

	done

	fastq_path="$pw/$indi"
	bam="${quantification}/${indi}/outs/possorted_genome_bam.bam"
	cd $quantification
	if [[ ! -d "$indi/outs" ]];then	
	$cell_ranger --id=${indi} --transcriptome=${transcriptome} --fastqs=${fastq_path} --sample=${indi} --localcores=16 --localmem=64
	fi

	source ~/scte_env/bin/activate
	$scte -i ${bam} -o ${indi} -x $index -UMI UB --min_genes 200 --thread 64

	cd ../
done
conda activate R
Rscript rscript_final.R "$quantification"
