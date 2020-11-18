#!/bin/bash -login
#PBS -l walltime=12:00:00,nodes=1:ppn=1,mem=4gb
#PBS -M annat.haber@jax.org
#PBS -m abe
#PBS -j oe

cd ${PBS_O_WORKDIR} 

module load bcftools/1.3.1
#module load vcftools/0.1.12a

InPath=/sdata/carter-lab/carter/AMPAD/WGS_PrelimProcessing/output/Joint/MAF05

#vcftools --gzvcf $InPath/Joint_allSNPjointMAF05.vcf.gz \
#		--snps allADVariantsApril2020.txt --recode --recode-INFO-all --out ADvariants

#| \
#                  bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' > \
#                  $OutPath/allADVariantsApril2020allSamples.txt

 #bcftools view -i'ID=@allADVariantsApril2020.txt' $InPath/Joint_allSNPjointMAF05.vcf.gz > $InPath/ADvariants_Joint_allSNPjointMAF05.vcf
 
 bcftools view -i'ID=@allADVariantsApril2020.txt' -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' \
 		$InPath/Joint_allSNPjointMAF05.vcf.gz > $InPath/allADVariantsApril2020allSamples.txt