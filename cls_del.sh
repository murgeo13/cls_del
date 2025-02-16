#!/bin/bash
echo "script has been started"
# load variables
source variables.txt
mkdir -p ${log_dir}
mkdir -p ${ref_dir}
mkdir -p ${out_dir}
cat <(echo -en "script has been started\t") <(date) 1>>${log_dir}/cls_del.log
# make indexes
echo "making of indexes"
bowtie2-build ${ref} ${indexes} &>>${log_dir}/cls_del.log
echo "indexes have been made"
cat <(echo -en "indexes have been made\t") <(date) 1>>${log_dir}/cls_del.log

# mapping reads
echo "mapping reads"
if [ -v reads ]
then
  bowtie2 -p 10 -x ${indexes} --interleaved ${reads} -S ${sam} 2>${log_dir}/cls_del.log
  echo "reads have been mapped"
  cat <(echo -en "reads have been mapped\t") <(date) 1>>${log_dir}/cls_del.log
else
  if [ -v in_F ] && [ -v in_R ]
  then
    bowtie2 -p 10 -x ${indexes} -1 ${in_F} -2 ${in_R} -S ${sam} 2>${log_dir}/cls_del.log
  echo "reads have been mapped"
  cat <(echo -en "reads have been mapped\t") <(date) 1>>${log_dir}/cls_del.log
  fi
fi

# convert file
echo "converting ${sam} to ${bam}"
samtools sort -o ${bam} ${sam} 2>>${log_dir}/cls_del.log
echo "${sam} has been converted to ${bam}"
cat <(echo -en "${sam} has been converted to ${bam}\t") <(date) 1>>${log_dir}/cls_del.log

# index bam
echo "indexing ${bam}"
samtools index ${bam} 2>${log_dir}/cls_del.log
echo "${bam} has been indexed"
cat <(echo -en "${bam} has been indexed\t") <(date) 1>>${log_dir}/cls_del.log

if [ -v cntrchr ]
then
  # get reads mapped to control chromosome
  echo "getting reads mapped to ${cntrchr}"
  samtools view -h -bS ${bam} ${cntrchr} 1> ${bam_on_cntrchr} 2>>${log_dir}/cls_del.log
  echo "reads mapped to ${cntrchr} are in ${bam_on_cntrchr}"
  cat <(echo -en "reads mapped to ${cntrchr} are in ${bam_on_cntrchr}\t") <(date) 1>>${log_dir}/cls_del.log
  # index bam
  echo "indexing ${bam_on_cntrchr}"
  samtools index ${bam_on_cntrchr} 2>${log_dir}/cls_del.log
  echo "${bam_on_cntrchr} has been indexed"
  cat <(echo -en "${bam_on_cntrchr} has been indexed\t") <(date) 1>>${log_dir}/cls_del.log
  echo "making control_cover.txt"
  samtools mpileup ${bam_on_cntrchr} -A --output-extra FLAG,POS,RNEXT,PNEXT -a -o ${out_dir}/control_cover.txt
  echo "control_cover.txt has been made"
  cat <(echo -en "control_cover.txt has been made\t") <(date) 1>>${log_dir}/cls_del.log
fi

# get reads mapped to mitochromosome
echo "getting reads mapped to ${mitochr}"
samtools view -h -bS ${bam} ${mitochr} 1> ${bam_on_mitochr} 2>>${log_dir}/cls_del.log
echo "reads mapped to ${mitochr} are in ${bam_on_mitochr}"
cat <(echo -en "reads mapped to ${mitochr} are in ${bam_on_mitochr}\t") <(date) 1>>${log_dir}/cls_del.log
# index bam
echo "indexing ${bam_on_mitochr}"
samtools index ${bam_on_mitochr} 2>${log_dir}/cls_del.log
echo "${bam_on_mitochr} has been indexed"
cat <(echo -en "${bam_on_mitochr} has been indexed\t") <(date) 1>>${log_dir}/cls_del.log
echo "making all_mito_cover.txt"
samtools mpileup ${bam_on_mitochr} -A --output-extra FLAG,POS,RNEXT,PNEXT -a -o ${out_dir}/all_mito_cover.txt
echo "all_mito_cover.txt has been made"
cat <(echo -en "all_mito_cover.txt has been made\t") <(date) 1>>${log_dir}/cls_del.log

# Delly
echo "discovering of SVs with Delly"
delly call -o ${delly_out_bcf} -g ${ref} ${bam_on_mitochr} &>>${log_dir}/cls_del.log
bcftools view ${delly_out_bcf} > ${delly_out_vcf} 2>>${log_dir}/cls_del.log
echo "Delly has been done"
cat <(echo -en "Delly has been done\t") <(date) 1>>${log_dir}/cls_del.log

# statistics with picard
echo "CollectInsertSizeMetrics"
picard CollectInsertSizeMetrics -I ${bam_on_mitochr} -O ${TLEN_txt} -H ${TLEN_pdf} -M 0.5 2>> ${log_dir}/cls_del.log
grep -v ^# ${TLEN_txt} | sed -n "2,3p" | datamash transpose --output-delimiter="=" 1> ${TLEN_var}

source ${TLEN_var}
if [ TLEN_trh = "default" ]
then
  TLEN_trh=${MEDIAN_INSERT_SIZE}+3*${MEDIAN_ABSOLUTE_DEVIATION}
fi
echo "getting reads mapped to ${mitochr} with large TLEN"
samtools view --input-fmt-option filter="tlen>${TLEN_trh} || tlen<-${TLEN_trh}" \
-bS ${bam_on_mitochr} 1> ${bam_filtered} 2>>${log_dir}/cls_del.log
echo "reads mapped to ${mitochr} with large TLEN are in ${bam_on_mitochr}"
cat <(echo -en "reads mapped to ${mitochr} with large TLEN are in ${bam_on_mitochr}\t") <(date) 1>>${log_dir}/cls_del.log

echo "making filtered_mito_cover.txt"
samtools mpileup ${bam_filtered} -A --output-extra FLAG,POS,RNEXT,PNEXT -a -o ${out_dir}/filtered_mito_cover.txt
echo "filtered_mito_cover.txt has been made"
cat <(echo -en "filtered_mito_cover.txt has been made\t") <(date) 1>>${log_dir}/cls_del.log


# python
if [ $(( mode & 1 )) = 1 ]
then
  python ./cover_plot.py --strain ${strain} --chr ${mitochr} --cntrchr ${cntrchr} --window 2*3*${MEDIAN_ABSOLUTE_DEVIATION} \
  --chrcov ${out_dir}/all_mito_cover.txt --cntrcov ${out_dir}/control_cover.txt --lang EN
  --out ./python_out 2>>${log_dir}/python.log
fi
mode=$(( mode >> 1 ))
if [ $(( mode & 1 )) = 1 ]
then
  python ./cls_del/main.py --chr ${mitochr} --chrcov ${out_dir}/filtered_mito_cover.txt --chrcov-total ${out_dir}/all_mito_cover.txt \
  --median ${MEDIAN_INSERT_SIZE} --sd ${MEDIAN_ABSOLUTE_DEVIATION} \
  --out ./python_out --save-temp-files --debug  2>>${log_dir}/python.log
fi
mode=$(( mode >> 1 ))
if [ $(( mode & 1 )) = 1 ]
then
  python ./cls_sel/for_delly.py --delly-vcf ${delly_out_vcf} --chrcov-total ${out_dir}/all_mito_cover.txt \
  --chrname ${mitochr} --window 2*3*${MEDIAN_ABSOLUTE_DEVIATION} \
  --out ./python_out --save-temp-files --debug
fi
