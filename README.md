# Timema Polygenic Selection

1. DNA sequence alignment

2. Variant calling and filtering

I used GATK (version 3.5) for variant calling. GatkFork.pl was used to control the alignment with the HaplotypeCaller.

```
module load gatk

cd /uufs/chpc.utah.edu/common/home/u6000989/data/timema/combind_wgs_dovetailV3/alignments_tcr

perl GatkFork.pl un*bam
```
Here is GatkFork.pl.
```
#!/usr/bin/perl
#
# make g.vcf files 
#

use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

my $idir = '/uufs/chpc.utah.edu/common/home/u6000989/data/timema/combind_wgs_dovetailV3/alignments_tcr';
my $odir = '/uufs/chpc.utah.edu/common/home/u6000989/data/timema/combind_wgs_dovetailV3/variants_tcr';
my $tdir = '/scratch/general/lustre/tcrGatk';
my $genome = "/uufs/chpc.utah.edu/common/home/u6000989/data/timema/tcrDovetail/version3/mod_map_timema_06Jun2016_RvNkF702.fasta";

FILES:
foreach $bam (@ARGV){
        $pm->start and next FILES; ## fork
        $out = $bam;
        $out =~ s/bam/g.vcf/ or die "failed here: $out\n";
        $out = "re_$out";
        system "java -Xmx48g -jar ~/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genome -I $idir/$bam -o $tdir/$out -XL na.intervals -gt_mode DISCOVERY -hets 0.001 -mbq 20 -out_mode EMIT_VARIANTS_ONLY -ploidy 2 -stand_call_conf 50 -pcrModel AGGRESSIVE --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000\n";
        system "cp $tdir/$out $odir/$out\n";


        $pm->finish;
}

$pm->wait_all_children;

```
GenotypeGVCFs was then used to for joint variant calling.
```
module load gatk

cd /uufs/chpc.utah.edu/common/home/u6000989/data/timema/combind_wgs_dovetailV3/variants_tcr

java -Xmx460g -jar ~/bin/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /uufs/chpc.utah.edu/common/home/u6000989/data/timema/tcrDovetail/version3/mod_map_timema_06Jun2016_RvNkF702.fasta -hets 0.001 -nt 32 --variant combinded_1.g.vcf --variant combinded_2.g.vcf --variant combinded_3.g.vcf --variant combinded_4.g.vcf --variant re_uniqe_merge_timemaHVC_8021008.g.vcf --variant re_uniqe_merge_timemaHVC_8021018.g.vcf -o tcr_wgs_variants_x.vcf

```

3. Genotype estimation

Next, I obtained maximum likelihood estimates of allele frequencies for all experimental samples using an expectation-maximization (EM) algorithm. `estpEM` is compiled from `main_estpEM.C` and `func_estpEM.C`.

```
estpEM -i mod_filtered1X_tcrExperimentVariants.gl -o pem_tcr_ecol_all_1X.txt -e 0.001 -m 30 -h 1
```

4. Estimation of total selection (selection differentials), ABC model fitting, and inference of LD were done with the C++ programs included in this repository.
