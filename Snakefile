# Snakefile for the 22126 NGS final project (group 14)
#
# Basically all the commands from our notes glued together so we don't have to
# run them one sample at a time and copy paste the paths every single time.
# Everything is hardcoded for the pupil server (/home/ctools, /home/databases,
# ...) so it will not run anywhere else unless the paths below are changed.
#
# how we ran it:
#   snakemake -j 1 -p --keep-going
#
# the fastq files in raw_data/ are already subsampled to the first 1M reads,
# see the README for that command.

SAMPLES = ["1GC", "3GC", "4GC", "6GC", "8GC", "9GC"]

PROJ = "/home/projects/22126_NGS/projects/group14/final_project"

# we mapped everything twice, once against the genome and once against the
# transcriptome, so most of the rules have a {ref} wildcard
REFS = {
    "genome": "/home/databases/references/human/GRCh38_full_analysis_set_plus_decoy_hla.fa",
    "transcriptome": PROJ + "/ref_transcriptome/transcriptome.fa",
}

DBSNP = "/home/databases/databases/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
MAPPABILITY = "/home/databases/databases/GRCh38/filter99.bed.gz"
SNPEFF_DATA = "/home/databases/databases/snpEff/"

# nothing is in $PATH on the server so we need the full path for most tools
FASTQC = "/home/ctools/FastQC/fastqc"
TRIMMOMATIC = "/home/ctools/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS = "/home/ctools/Trimmomatic-0.39/adapters/TruSeq2-PE.fa"
PICARD = "/home/ctools/picard_2.23.8/picard.jar"
GATK = "/home/ctools/gatk-4.6.2.0/gatk"
TABIX = "/home/ctools/htslib-1.20/tabix"
BGZIP = "/home/ctools/htslib-1.20/bgzip"
BCFTOOLS = "/home/ctools/bcftools-1.23/bcftools"
SNPEFF = "/home/ctools/snpEff/snpEff.jar"
# bwa, samtools and bedtools were already there so we just call them directly

# without this the two fastqc rules are ambiguous, because "1_pair" also
# matches the {r} wildcard
wildcard_constraints:
    sample = "[0-9]+GC",
    r = "[12]",
    ref = "genome|transcriptome"


rule all:
    input:
        expand("fastqc/{sample}_sub{r}_fastqc.html", sample=SAMPLES, r=[1, 2]),
        expand("fastqc/{sample}_sub{r}_pair_fastqc.html", sample=SAMPLES, r=[1, 2]),
        expand("annotation/{sample}_{ref}_annotation.vcf.gz", sample=SAMPLES, ref=REFS),
        expand("stats/{sample}_{ref}_filtered.txt", sample=SAMPLES, ref=REFS),
        expand("stats/{sample}_{ref}_pass_map99.txt", sample=SAMPLES, ref=REFS),


# --------------------------------------------------------- QC and trimming

rule fastqc_raw:
    input:
        "raw_data/{sample}_sub{r}.fastq.gz"
    output:
        "fastqc/{sample}_sub{r}_fastqc.html",
        "fastqc/{sample}_sub{r}_fastqc.zip"
    shell:
        "{FASTQC} {input} -o fastqc"


rule trimmomatic:
    input:
        r1 = "raw_data/{sample}_sub1.fastq.gz",
        r2 = "raw_data/{sample}_sub2.fastq.gz"
    output:
        p1 = "trimmed/{sample}_sub1_pair.fastq.gz",
        u1 = "trimmed/{sample}_sub1_unpair.fastq.gz",
        p2 = "trimmed/{sample}_sub2_pair.fastq.gz",
        u2 = "trimmed/{sample}_sub2_unpair.fastq.gz"
    shell:
        "java -jar {TRIMMOMATIC} PE -threads 1 -phred33 "
        "{input.r1} {input.r2} "
        "{output.p1} {output.u1} {output.p2} {output.u2} "
        "ILLUMINACLIP:{ADAPTERS}:2:30:10 "
        "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:50"


# same thing as above but on the trimmed reads, to check the adapters are gone
rule fastqc_trimmed:
    input:
        "trimmed/{sample}_sub{r}_pair.fastq.gz"
    output:
        "fastqc/{sample}_sub{r}_pair_fastqc.html",
        "fastqc/{sample}_sub{r}_pair_fastqc.zip"
    shell:
        "{FASTQC} {input} -o fastqc"


# --------------------------------------------------------- mapping

rule bwa_mem:
    input:
        p1 = "trimmed/{sample}_sub1_pair.fastq.gz",
        p2 = "trimmed/{sample}_sub2_pair.fastq.gz"
    output:
        "mapped/{sample}_{ref}.bam"
    params:
        fa = lambda w: REFS[w.ref]
    shell:
        "bwa mem {params.fa} {input.p1} {input.p2} "
        "| samtools view -uS - | samtools sort /dev/stdin > {output}"


# this sorts a second time which is a bit pointless since the pipe above
# already sorts, but we kept it because that is how we did it in the exercise
rule sort:
    input:
        "mapped/{sample}_{ref}.bam"
    output:
        "sorted/{sample}_{ref}_sorted.bam"
    shell:
        "samtools sort {input} > {output}"


# --------------------------------------------------------- postprocessing

rule mark_duplicates:
    input:
        "sorted/{sample}_{ref}_sorted.bam"
    output:
        bam = "dedup/{sample}_{ref}_duplicates.bam",
        metrics = "dedup/{sample}_{ref}_metrics.txt"
    shell:
        "java -jar {PICARD} MarkDuplicates "
        "-I {input} -O {output.bam} -M {output.metrics}"


rule add_read_groups:
    input:
        "dedup/{sample}_{ref}_duplicates.bam"
    output:
        "dedup/{sample}_{ref}_duplicates.RG.bam"
    shell:
        "java -jar {PICARD} AddOrReplaceReadGroups -I {input} -O {output} "
        "-RGID 1 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM {wildcards.sample}"


rule index_bam:
    input:
        "dedup/{sample}_{ref}_duplicates.RG.bam"
    output:
        "dedup/{sample}_{ref}_duplicates.RG.bam.bai"
    shell:
        "samtools index {input}"


# --------------------------------------------------------- variant calling

# in the notes we used the genome fasta here even for the transcriptome bam,
# that was a copy paste mistake, so this takes whichever reference the reads
# were actually mapped against
rule haplotype_caller:
    input:
        bam = "dedup/{sample}_{ref}_duplicates.RG.bam",
        bai = "dedup/{sample}_{ref}_duplicates.RG.bam.bai"
    output:
        "gvcf/{sample}_{ref}.gvcf.gz"
    params:
        fa = lambda w: REFS[w.ref]
    shell:
        '{GATK} --java-options "-Xmx10g" HaplotypeCaller '
        "-R {params.fa} -I {input.bam} -O {output} "
        "--dbsnp {DBSNP} -ERC GVCF"


rule index_gvcf:
    input:
        "gvcf/{sample}_{ref}.gvcf.gz"
    output:
        "gvcf/{sample}_{ref}.gvcf.gz.tbi"
    shell:
        "{TABIX} -f -p vcf {input}"


rule genotype_gvcf:
    input:
        gvcf = "gvcf/{sample}_{ref}.gvcf.gz",
        tbi = "gvcf/{sample}_{ref}.gvcf.gz.tbi"
    output:
        "vcf/{sample}_{ref}.vcf.gz"
    params:
        fa = lambda w: REFS[w.ref]
    shell:
        "{GATK} GenotypeGVCFs -R {params.fa} -V {input.gvcf} -O {output} "
        "--dbsnp {DBSNP}"


# --------------------------------------------------------- filtering

rule hard_filter:
    input:
        "vcf/{sample}_{ref}.vcf.gz"
    output:
        "hard_filtering/{sample}_{ref}_filtering.vcf.gz"
    shell:
        "{GATK} VariantFiltration -V {input} -O {output} "
        '-filter "DP < 10.0" --filter-name "DP" '
        '-filter "QUAL < 30.0" --filter-name "QUAL30" '
        '-filter "SOR > 3.0" --filter-name "SOR3" '
        '-filter "FS > 60.0" --filter-name "FS60" '
        '-filter "MQ < 40.0" --filter-name "MQ40"'


# how many sites were filtered out, and which filters failed.
# the "|| true" is there because grep exits with 1 when it finds nothing and
# then snakemake kills the whole run
rule filter_stats:
    input:
        "hard_filtering/{sample}_{ref}_filtering.vcf.gz"
    output:
        "stats/{sample}_{ref}_filtered.txt"
    shell:
        "( {BCFTOOLS} view -H {input} | grep -v PASS | wc -l ; "
        "{BCFTOOLS} view -H {input} | grep -v PASS | cut -f7 "
        "| sort | uniq -c | sort -n ) > {output} || true"


rule mappability_filter:
    input:
        "hard_filtering/{sample}_{ref}_filtering.vcf.gz"
    output:
        "hard_filtering/{sample}_{ref}_filtering_map99.vcf.gz"
    shell:
        "bedtools intersect -header -a {input} -b {MAPPABILITY} "
        "| {BGZIP} -c > {output}"


rule pass_after_mappability:
    input:
        "hard_filtering/{sample}_{ref}_filtering_map99.vcf.gz"
    output:
        "stats/{sample}_{ref}_pass_map99.txt"
    shell:
        "{BCFTOOLS} view -H {input} | grep PASS | wc -l > {output} || true"


# --------------------------------------------------------- annotation

rule snpeff:
    input:
        "hard_filtering/{sample}_{ref}_filtering_map99.vcf.gz"
    output:
        vcf = "annotation/{sample}_{ref}_annotation.vcf.gz",
        html = "annotation/{sample}_{ref}_annotation.html"
    shell:
        "java -jar {SNPEFF} eff -dataDir {SNPEFF_DATA} "
        "-htmlStats {output.html} GRCh38.99 {input} "
        "| {BGZIP} -c > {output.vcf}"


# we also counted the total number of annotations per sample but that we did
# by hand afterwards, see count_annotations.sh
