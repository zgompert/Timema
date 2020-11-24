# Timema Polygenic Selection

```
module load gatk

cd /uufs/chpc.utah.edu/common/home/u6000989/data/timema/combind_wgs_dovetailV3/variants_tcr

java -Xmx460g -jar ~/bin/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /uufs/chpc.utah.edu/common/home/u6000989/data/timema/tcrDovetail/version3/mod_map_timema_06Jun2016_RvNkF702.fasta -hets 0.001 -nt 32 --variant combinded_1.g.vcf --variant combinded_2.g.vcf --variant combinded_3.g.vcf --variant combinded_4.g.vcf --variant re_uniqe_merge_timemaHVC_8021008.g.vcf --variant re_uniqe_merge_timemaHVC_8021018.g.vcf -o tcr_wgs_variants_x.vcf

```
