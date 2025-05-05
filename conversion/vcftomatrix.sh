#!/bin/bash

#SBATCH -J vmat
#SBATCH -o VCFTOMAT.%j.out
#SBATCH -e VCFTOMAT.%j.error
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch
#SBATCH --qos=240c-1h_batch
#SBATCH --mail-user=r.d.pasco@irri.org
#SBATCH --mail-type=ALL
#SBATCH --requeue

invcf=169.vcf.gz

pref=169

mkdir -p $pref

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' $invcf |  tr "|" "/"  | \
   tr -d "/" | \
   sed "s:chr0\?::"  | \
   awk -v OFS="\t" '{$1 = $1 + 0; print $0}'  > ${pref}/mat_vcf.txt

cut -f1-4 $pref/mat_vcf.txt > $pref/pos.txt

# The following fails for some reason (grep CHR part)

bcftools view -o -   -h $invcf  |  grep CHR | tr "\t" "\n"  | tail -n+10 > $pref/sample_list.txt
