#!/bin/bash
echo "script is started"
# берём переменные из конфиг файла
source variables.txt
mkdir -p ${log_dir}
mkdir -p ${ref_dir}
mkdir -p ${out_dir}
# make indexes
bowtie2-build ${ref} ${indexes} &>${log_dir}/bowtie2-build.log
echo "indexes are made"

# картируем чтения
bowtie2 -p 10 -x ${indexes} -1 ${in_F} -2 ${in_R} -S ${sam} 2>${log_dir}/bowtie2.log
echo "reads are mapped"
# конвертация файла
echo "конвертация файла" 1>> ${log_dir}/samtools.log
samtools sort -o ${bam} ${sam} 2>>${log_dir}/samtools.log

# индексируем bam файл
echo "индексируем bam файл" 1>> ${log_dir}/samtools.log
samtools index ${bam} 2>${log_dir}/bam_index.log
echo "bam is indexed"

# Получение чтений, картированных на вашу хромосому
echo "Получение чтений, картированных на хромосому "${chr} 1>> ${log_dir}/samtools.log
samtools view -h -bS ${bam} ${chr} 1> ${bam_on_chr} 2>>${log_dir}/samtools.log
# индексируем правильный bam файл
echo "индексируем правильный bam файл" 1>> ${log_dir}/samtools.log
samtools index ${bam_on_chr}
samtools mpileup ${bam_on_chr} -A --output-extra FLAG,POS,RNEXT,PNEXT -a -o all_cover.txt

#используем delly
echo "Поиск перестроек с помощью Delly"
delly call -o ${delly_out_bcf} -g ${ref} ${bam_on_chr} &>>${log_dir}/delly.log
bcftools view ${delly_out_bcf} > ${delly_out_vcf} 2>>${log_dir}/delly.log
echo "Delly is done"

#статистика picard
echo "Собрать статистику по размеру вставки между чтениями"
picard CollectInsertSizeMetrics -I ${bam_on_chr} -O ${TLEN_txt} -H ${TLEN_pdf} -M 0.5 2>> ${log_dir}/picard.log
grep -v ^# ${TLEN_txt} | sed -n "2,3p" | datamash transpose --output-delimiter="=" 1> ${TLEN_var}

source ${TLEN_var}
echo "Получение чтений, картированных на хромосому ${chr} c большим TLEN" 1>> ${log_dir}/samtools.log
samtools view --input-fmt-option filter="tlen>${MEDIAN_INSERT_SIZE}+3*${MEDIAN_ABSOLUTE_DEVIATION} || tlen<-${MEDIAN_INSERT_SIZE}-3*${MEDIAN_ABSOLUTE_DEVIATION}" \
-bS ${bam_on_chr} 1> ${bam_filtered} 2>>${log_dir}/samtools.log

samtools mpileup ${bam_filtered} -A --output-extra FLAG,POS,RNEXT,PNEXT -a -o filtered_cover.txt

#python
echo "python is started"
python ./cls_del/main.py --chr ${chr} --chrcov filtered_cover.txt --chrcov-total all_cover.txt--median ${MEDIAN_INSERT_SIZE} --sd ${MEDIAN_ABSOLUTE_DEVIATION} \
--out ./python_out_new --save-temp-files --debug
