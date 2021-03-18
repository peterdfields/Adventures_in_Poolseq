#!/bin/sh
# The present script can be run to do UnifiedGenotyper on all bams created by
# align_step.sh, so should have *_mapped_rdub_indelrealigner.bam suffix
# there are multithreading options not included here, and there could also be
# modifications to allow for a large number of files (eg memory allocations)

REF=some.fasta

java -jar -T GenomeAnalysisTK.jar -T UnifiedGenotyper -mbq 20 -R "$REF" \
$(for file in *_mapped_rdub_indelrealigner.bam; do echo "-I $file "; done) \
-o out_UG.vcf 