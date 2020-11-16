#!/bin/bash -login
#PBS -l walltime=24:00:00,nodes=1:ppn=1,mem=8gb
#PBS -M annat.haber@jax.org
#PBS -m abe
#PBS -j oe

cd ${PBS_O_WORKDIR} 

module load bcftools/1.3.1
#module load vcftools/0.1.12a

InPath=/sdata/carter-lab/carter/AMPAD/WGS_PrelimProcessing/output/Joint/MAF05
inputvcf=NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_8.recalibrated_variants.vcf.gz

bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' -i 'POS=27219987' $InPath/$inputvcf >> \
                $InPath/ADvariants/allADVariantsApril2020allSamples_chr8.txt

bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' -i 'POS=145154222' $InPath/$inputvcf >> \
                $InPath/ADvariants/allADVariantsApril2020allSamples_chr8.txt
