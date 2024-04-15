#!/bin/bash
echo "script is started"
# load variables
source variables.txt
mkdir -p ${log_dir}
mkdir -p ${ref_dir}
mkdir -p ${out_dir}
# make indexes
bowtie2-build ${ref} ${indexes} &>${log_dir}/bowtie2-build.log
echo "indexes are made"

# map reads
bowtie2 -p 10 -x ${indexes} -1 ${in_F} -2 ${in_R} -S ${sam} 2>${log_dir}/bowtie2.log
echo "reads are maped"
# convert file
echo "convertation ${sam} to ${bam}" 1>> ${log_dir}/samtools.log
samtools sort -o ${bam} ${sam} 2>>${log_dir}/samtools.log

# index bam
echo "indexing of ${bam}" 1>> ${log_dir}/samtools.log
samtools index ${bam} 2>${log_dir}/bam_index.log
echo "${bam} is indexed"

# get reads mapped to chromosome
echo "getting reads mapped to ${chr}" 1>> ${log_dir}/samtools.log
samtools view -h -bS ${bam} ${chr} 1> ${bam_on_chr} 2>>${log_dir}/samtools.log
# index bam
echo "indexing of ${bam_on_chr}" 1>> ${log_dir}/samtools.log
samtools index ${bam_on_chr}
samtools mpileup ${bam_on_chr} -A --output-extra FLAG,POS,RNEXT,PNEXT -a -o ${out_dir}/all_cover.txt

# Delly
echo "discovering of SVs with Delly"
delly call -o ${delly_out_bcf} -g ${ref} ${bam_on_chr} &>>${log_dir}/delly.log
bcftools view ${delly_out_bcf} > ${delly_out_vcf} 2>>${log_dir}/delly.log
echo "Delly is done"

# statistics with picard
echo "CollectInsertSizeMetrics"
picard CollectInsertSizeMetrics -I ${bam_on_chr} -O ${TLEN_txt} -H ${TLEN_pdf} -M 0.5 2>> ${log_dir}/picard.log
grep -v ^# ${TLEN_txt} | sed -n "2,3p" | datamash transpose --output-delimiter="=" 1> ${TLEN_var}

source ${TLEN_var}
echo "getting reads mapped to ${chr} with large TLEN" 1>> ${log_dir}/samtools.log
samtools view --input-fmt-option filter="tlen>${MEDIAN_INSERT_SIZE}+3*${MEDIAN_ABSOLUTE_DEVIATION} || tlen<-${MEDIAN_INSERT_SIZE}-3*${MEDIAN_ABSOLUTE_DEVIATION}" \
-bS ${bam_on_chr} 1> ${bam_filtered} 2>>${log_dir}/samtools.log

samtools mpileup ${bam_filtered} -A --output-extra FLAG,POS,RNEXT,PNEXT -a -o ${out_dir}/filtered_cover.txt

# python
python ./cls_del/main.py --chr ${chr} --chrcov ${out_dir}/filtered_cover.txt --chrcov-total ${out_dir}/all_cover.txt \
--median ${MEDIAN_INSERT_SIZE} --sd ${MEDIAN_ABSOLUTE_DEVIATION} \
--out ./python_out_new --save-temp-files --debug  2>>${log_dir}/python.log
