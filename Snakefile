#Snakefile for the 22126 NGS final project (group 14)
#All the commands from our notes put together so we don't have to run them
#one sample at a time and copy paste the paths every time.

#Everything is hardcoded for the pupil server (/home/ctools, /home/databases),
#so it only runs there.

#How we ran it:
#	snakemake -np
#	snakemake --jobs 1

#The files in raw_data/ are already cut down to the first 1M reads, see README.

#We mapped everything twice, once against the genome and once against the
#transcriptome, so most steps are in here twice.

ids = ["1GC", "3GC", "4GC", "6GC", "8GC", "9GC"]

rule all:
	input:
		expand("fastqc/raw/{sample}_sub1_fastqc.html", sample=ids),
		expand("fastqc/trimmed/{sample}_sub1_pair_fastqc.html", sample=ids),
		expand("annotation/{sample}_genome_annotation.vcf.gz", sample=ids),
		expand("annotation/{sample}_transcriptome_annotation.vcf.gz", sample=ids)


#### Quality control and trimming ####

rule quality_check_raw:
	input:
		read_1="raw_data/{sample}_sub1.fastq.gz",
		read_2="raw_data/{sample}_sub2.fastq.gz"
	output:
		"fastqc/raw/{sample}_sub1_fastqc.html",
		"fastqc/raw/{sample}_sub2_fastqc.html"
	params:
		outdir="fastqc/raw"
	message: "Quality check of the raw reads of {wildcards.sample} with FastQC"
	shell:
		"""
		/home/ctools/FastQC/fastqc {input} -o {params.outdir}
		"""

rule trimming:
	input:
		read_1="raw_data/{sample}_sub1.fastq.gz",
		read_2="raw_data/{sample}_sub2.fastq.gz"
	output:
		pair_1="trimmed/{sample}_sub1_pair.fastq.gz",
		unpair_1="trimmed/{sample}_sub1_unpair.fastq.gz",
		pair_2="trimmed/{sample}_sub2_pair.fastq.gz",
		unpair_2="trimmed/{sample}_sub2_unpair.fastq.gz"
	params:
		adapters="/home/ctools/Trimmomatic-0.39/adapters/TruSeq2-PE.fa"
	message: "Trimming the adapters and bad quality bases of {wildcards.sample}"
	shell:
		"""
		java -jar /home/ctools/Trimmomatic-0.39/trimmomatic-0.39.jar \
		    PE -threads 1 -phred33 \
		    {input.read_1} {input.read_2} \
		    {output.pair_1} {output.unpair_1} \
		    {output.pair_2} {output.unpair_2} \
		    ILLUMINACLIP:{params.adapters}:2:30:10 \
		    LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:50
		"""

#Same as above but on the trimmed reads, to see if the adapters are gone.
rule quality_check_trimmed:
	input:
		read_1="trimmed/{sample}_sub1_pair.fastq.gz",
		read_2="trimmed/{sample}_sub2_pair.fastq.gz"
	output:
		"fastqc/trimmed/{sample}_sub1_pair_fastqc.html",
		"fastqc/trimmed/{sample}_sub2_pair_fastqc.html"
	params:
		outdir="fastqc/trimmed"
	message: "Quality check of the trimmed reads of {wildcards.sample} with FastQC"
	shell:
		"""
		/home/ctools/FastQC/fastqc {input} -o {params.outdir}
		"""


#### Mapping ####

rule mapping_genome:
	input:
		read_1="trimmed/{sample}_sub1_pair.fastq.gz",
		read_2="trimmed/{sample}_sub2_pair.fastq.gz"
	output:
		"mapped/{sample}_genome.bam"
	params:
		ref="/home/databases/references/human/GRCh38_full_analysis_set_plus_decoy_hla.fa"
	message: "Mapping {wildcards.sample} against the genome with bwa mem"
	shell:
		"""
		bwa mem {params.ref} {input.read_1} {input.read_2} \
		    | samtools view -uS - | samtools sort /dev/stdin > {output}
		"""

rule mapping_transcriptome:
	input:
		read_1="trimmed/{sample}_sub1_pair.fastq.gz",
		read_2="trimmed/{sample}_sub2_pair.fastq.gz"
	output:
		"mapped/{sample}_transcriptome.bam"
	params:
		ref="/home/projects/22126_NGS/projects/group14/final_project/ref_transcriptome/transcriptome.fa"
	message: "Mapping {wildcards.sample} against the transcriptome with bwa mem"
	shell:
		"""
		bwa mem {params.ref} {input.read_1} {input.read_2} \
		    | samtools view -uS - | samtools sort /dev/stdin > {output}
		"""


#### Postprocessing ####

#This sorts a second time even though the mapping already sorts, but we kept
#it because that is how we did it in the exercise.

rule sort_genome:
	input:
		"mapped/{sample}_genome.bam"
	output:
		"sorted/{sample}_genome_sorted.bam"
	message: "Sorting the genome bam file of {wildcards.sample}"
	shell:
		"""
		samtools sort {input} > {output}
		"""

rule sort_transcriptome:
	input:
		"mapped/{sample}_transcriptome.bam"
	output:
		"sorted/{sample}_transcriptome_sorted.bam"
	message: "Sorting the transcriptome bam file of {wildcards.sample}"
	shell:
		"""
		samtools sort {input} > {output}
		"""

rule mark_duplicates_genome:
	input:
		"sorted/{sample}_genome_sorted.bam"
	output:
		bam="dedup/{sample}_genome_duplicates.bam",
		metrics="dedup/{sample}_genome_metrics.txt"
	message: "Marking the duplicates in the genome bam file of {wildcards.sample}"
	shell:
		"""
		java -jar /home/ctools/picard_2.23.8/picard.jar MarkDuplicates \
		    -I {input} -O {output.bam} -M {output.metrics}
		"""

rule mark_duplicates_transcriptome:
	input:
		"sorted/{sample}_transcriptome_sorted.bam"
	output:
		bam="dedup/{sample}_transcriptome_duplicates.bam",
		metrics="dedup/{sample}_transcriptome_metrics.txt"
	message: "Marking the duplicates in the transcriptome bam file of {wildcards.sample}"
	shell:
		"""
		java -jar /home/ctools/picard_2.23.8/picard.jar MarkDuplicates \
		    -I {input} -O {output.bam} -M {output.metrics}
		"""

rule read_groups_genome:
	input:
		"dedup/{sample}_genome_duplicates.bam"
	output:
		"dedup/{sample}_genome_duplicates.RG.bam"
	message: "Adding the read groups to the genome bam file of {wildcards.sample}"
	shell:
		"""
		java -jar /home/ctools/picard_2.23.8/picard.jar AddOrReplaceReadGroups \
		    -I {input} -O {output} \
		    -RGID 1 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM {wildcards.sample}
		"""

rule read_groups_transcriptome:
	input:
		"dedup/{sample}_transcriptome_duplicates.bam"
	output:
		"dedup/{sample}_transcriptome_duplicates.RG.bam"
	message: "Adding the read groups to the transcriptome bam file of {wildcards.sample}"
	shell:
		"""
		java -jar /home/ctools/picard_2.23.8/picard.jar AddOrReplaceReadGroups \
		    -I {input} -O {output} \
		    -RGID 1 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM {wildcards.sample}
		"""

rule index_bam_genome:
	input:
		"dedup/{sample}_genome_duplicates.RG.bam"
	output:
		"dedup/{sample}_genome_duplicates.RG.bam.bai"
	message: "Indexing the genome bam file of {wildcards.sample}"
	shell:
		"""
		samtools index {input}
		"""

rule index_bam_transcriptome:
	input:
		"dedup/{sample}_transcriptome_duplicates.RG.bam"
	output:
		"dedup/{sample}_transcriptome_duplicates.RG.bam.bai"
	message: "Indexing the transcriptome bam file of {wildcards.sample}"
	shell:
		"""
		samtools index {input}
		"""


#### Variant calling ####

#In our notes we had the genome fasta here for the transcriptome as well,
#but that was a copy paste mistake, so each one uses its own reference.

rule haplotype_caller_genome:
	input:
		bam="dedup/{sample}_genome_duplicates.RG.bam",
		bai="dedup/{sample}_genome_duplicates.RG.bam.bai"
	output:
		"gvcf/{sample}_genome.gvcf.gz"
	params:
		ref="/home/databases/references/human/GRCh38_full_analysis_set_plus_decoy_hla.fa",
		dbsnp="/home/databases/databases/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
	message: "Calling the variants of {wildcards.sample} against the genome"
	shell:
		"""
		/home/ctools/gatk-4.6.2.0/gatk --java-options "-Xmx10g" HaplotypeCaller \
		    -R {params.ref} \
		    -I {input.bam} \
		    -O {output} \
		    --dbsnp {params.dbsnp} \
		    -ERC GVCF
		"""

rule haplotype_caller_transcriptome:
	input:
		bam="dedup/{sample}_transcriptome_duplicates.RG.bam",
		bai="dedup/{sample}_transcriptome_duplicates.RG.bam.bai"
	output:
		"gvcf/{sample}_transcriptome.gvcf.gz"
	params:
		ref="/home/projects/22126_NGS/projects/group14/final_project/ref_transcriptome/transcriptome.fa",
		dbsnp="/home/databases/databases/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
	message: "Calling the variants of {wildcards.sample} against the transcriptome"
	shell:
		"""
		/home/ctools/gatk-4.6.2.0/gatk --java-options "-Xmx10g" HaplotypeCaller \
		    -R {params.ref} \
		    -I {input.bam} \
		    -O {output} \
		    --dbsnp {params.dbsnp} \
		    -ERC GVCF
		"""

rule index_gvcf_genome:
	input:
		"gvcf/{sample}_genome.gvcf.gz"
	output:
		"gvcf/{sample}_genome.gvcf.gz.tbi"
	message: "Indexing the genome gvcf file of {wildcards.sample}"
	shell:
		"""
		/home/ctools/htslib-1.20/tabix -f -p vcf {input}
		"""

rule index_gvcf_transcriptome:
	input:
		"gvcf/{sample}_transcriptome.gvcf.gz"
	output:
		"gvcf/{sample}_transcriptome.gvcf.gz.tbi"
	message: "Indexing the transcriptome gvcf file of {wildcards.sample}"
	shell:
		"""
		/home/ctools/htslib-1.20/tabix -f -p vcf {input}
		"""

rule genotype_genome:
	input:
		gvcf="gvcf/{sample}_genome.gvcf.gz",
		tbi="gvcf/{sample}_genome.gvcf.gz.tbi"
	output:
		"vcf/{sample}_genome.vcf.gz"
	params:
		ref="/home/databases/references/human/GRCh38_full_analysis_set_plus_decoy_hla.fa",
		dbsnp="/home/databases/databases/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
	message: "Making the genotype file of {wildcards.sample} for the genome"
	shell:
		"""
		/home/ctools/gatk-4.6.2.0/gatk GenotypeGVCFs \
		    -R {params.ref} \
		    -V {input.gvcf} \
		    -O {output} \
		    --dbsnp {params.dbsnp}
		"""

rule genotype_transcriptome:
	input:
		gvcf="gvcf/{sample}_transcriptome.gvcf.gz",
		tbi="gvcf/{sample}_transcriptome.gvcf.gz.tbi"
	output:
		"vcf/{sample}_transcriptome.vcf.gz"
	params:
		ref="/home/projects/22126_NGS/projects/group14/final_project/ref_transcriptome/transcriptome.fa",
		dbsnp="/home/databases/databases/GRCh38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
	message: "Making the genotype file of {wildcards.sample} for the transcriptome"
	shell:
		"""
		/home/ctools/gatk-4.6.2.0/gatk GenotypeGVCFs \
		    -R {params.ref} \
		    -V {input.gvcf} \
		    -O {output} \
		    --dbsnp {params.dbsnp}
		"""


#### Filtering ####

rule hard_filtering_genome:
	input:
		"vcf/{sample}_genome.vcf.gz"
	output:
		"hard_filtering/{sample}_genome_filtering.vcf.gz"
	message: "Hard filtering the genome variants of {wildcards.sample}"
	shell:
		"""
		/home/ctools/gatk-4.6.2.0/gatk VariantFiltration -V {input} -O {output} \
		    -filter "DP < 10.0" --filter-name "DP" \
		    -filter "QUAL < 30.0" --filter-name "QUAL30" \
		    -filter "SOR > 3.0" --filter-name "SOR3" \
		    -filter "FS > 60.0" --filter-name "FS60" \
		    -filter "MQ < 40.0" --filter-name "MQ40"
		"""

rule hard_filtering_transcriptome:
	input:
		"vcf/{sample}_transcriptome.vcf.gz"
	output:
		"hard_filtering/{sample}_transcriptome_filtering.vcf.gz"
	message: "Hard filtering the transcriptome variants of {wildcards.sample}"
	shell:
		"""
		/home/ctools/gatk-4.6.2.0/gatk VariantFiltration -V {input} -O {output} \
		    -filter "DP < 10.0" --filter-name "DP" \
		    -filter "QUAL < 30.0" --filter-name "QUAL30" \
		    -filter "SOR > 3.0" --filter-name "SOR3" \
		    -filter "FS > 60.0" --filter-name "FS60" \
		    -filter "MQ < 40.0" --filter-name "MQ40"
		"""

rule mappability_genome:
	input:
		"hard_filtering/{sample}_genome_filtering.vcf.gz"
	output:
		"hard_filtering/{sample}_genome_filtering_map99.vcf.gz"
	params:
		bed="/home/databases/databases/GRCh38/filter99.bed.gz"
	message: "Filtering the genome variants of {wildcards.sample} by mappability"
	shell:
		"""
		bedtools intersect -header -a {input} -b {params.bed} \
		    | /home/ctools/htslib-1.20/bgzip -c > {output}
		"""

rule mappability_transcriptome:
	input:
		"hard_filtering/{sample}_transcriptome_filtering.vcf.gz"
	output:
		"hard_filtering/{sample}_transcriptome_filtering_map99.vcf.gz"
	params:
		bed="/home/databases/databases/GRCh38/filter99.bed.gz"
	message: "Filtering the transcriptome variants of {wildcards.sample} by mappability"
	shell:
		"""
		bedtools intersect -header -a {input} -b {params.bed} \
		    | /home/ctools/htslib-1.20/bgzip -c > {output}
		"""


#### Annotation ####

rule annotation_genome:
	input:
		"hard_filtering/{sample}_genome_filtering_map99.vcf.gz"
	output:
		vcf="annotation/{sample}_genome_annotation.vcf.gz",
		html="annotation/{sample}_genome_annotation.html"
	params:
		datadir="/home/databases/databases/snpEff/"
	message: "Annotating the genome variants of {wildcards.sample} with snpEff"
	shell:
		"""
		java -jar /home/ctools/snpEff/snpEff.jar eff \
		    -dataDir {params.datadir} \
		    -htmlStats {output.html} \
		    GRCh38.99 {input} \
		    | /home/ctools/htslib-1.20/bgzip -c > {output.vcf}
		"""

rule annotation_transcriptome:
	input:
		"hard_filtering/{sample}_transcriptome_filtering_map99.vcf.gz"
	output:
		vcf="annotation/{sample}_transcriptome_annotation.vcf.gz",
		html="annotation/{sample}_transcriptome_annotation.html"
	params:
		datadir="/home/databases/databases/snpEff/"
	message: "Annotating the transcriptome variants of {wildcards.sample} with snpEff"
	shell:
		"""
		java -jar /home/ctools/snpEff/snpEff.jar eff \
		    -dataDir {params.datadir} \
		    -htmlStats {output.html} \
		    GRCh38.99 {input} \
		    | /home/ctools/htslib-1.20/bgzip -c > {output.vcf}
		"""
