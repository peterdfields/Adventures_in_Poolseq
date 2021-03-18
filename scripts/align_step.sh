# This bash script does:
# 	mapping reads against the reference
#	claculating simple mapping stats
#	converting to bam file and sort bam
#	removing dublicates
#	indelrealignment
#	calling variants (SNP) with unified genotyper (GATK)

# version 1
# 18.03.2021
# author: meret j. halter with modifications by peter d. fields

#!/bin/sh

sample=$1
describer=$(echo ${sample} | sed 's/.qc.fq.gz//')
REF=ref.fasta

# mapping reads to reference
bwa mem -t 8 -p -M  "$REF" ${describer}.qc.fq.gz | gzip -3 > ${describer}.sam.gz

# print mapping summaries
samtools flagstat ${describer}.sam.gz > ${describer}_ref_flagstat.txt

# convert to bam
samtools view ${describer}.sam.gz -b -@ 10 -o ${describer}.bam

# sort bam
samtools sort ${describer}.bam -o ${describer}.sort.bam

# add readgroup id
java -jar picard.jar AddOrReplaceReadGroups \
I=${describer}.sort.bam \
O=${describer}.bam \
RGID=${describer}  \
RGLB=${describer}  \
RGPL=illumina  \
RGPU=${describer}  \
RGSM=${describer}

# index readgrouped bam
samtools index ${describer}.bam

# remove duplicates
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE \ 
I=${describer}.bam O=${describer}_rdub.bam M=${describer}_ref_dup_metrics.txt

# index duplicate removed bam
samtools index ${describer}_rdub.bam

#remove intermediate files
rm ${describer}.sam.gz
rm ${describer}.bam
rm ${describer}.sort.bam

# call indels for doing realignment
java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R "$REF" \
-I ${describer}_rdub.bam \
-glm INDEL \
-o ${describer}_INDEL_ref.vcf

# create target intervals list using RealignerTargetCreator
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R "$REF"
-known ${describer}_INDEL_ref.vcf
-I ${describer}_rdub.bam
-o ${describer}_realignertargetcreator.intervals

# Realign reads using IndelRealigner
## be careful with memory allowance here, as well as tmp directory assignment
java -Xmx32G -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar \
-T IndelRealigner \
-R "$REF" \
-targetIntervals ${describer}_realignertargetcreator.intervals \
-known ${describer}_INDEL_Xinb3.vcf \
-I ${describer}_rdub.bam \
-o ${describer}_rdub_indelrealigner.bam

# removing unmapped reads
samtools view -b -F 4 ${describer}_rdub_indelrealigner.bam > \ 
${describer}_mapped_rdub_indelrealigner.bam

# index bam file
samtools index ${describer}_mapped_rdub_indelrealigner.bam

# remove intermediate files
rm ${describer}_rdub_indelrealigner.bam
rm ${describer}_rdub.bam
rm ${describer}_rdub.bam.bai

# move to joint variant calling step...