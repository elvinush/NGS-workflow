# NGS final project - group 14

Course 22126 (Next Generation Sequencing Analysis) at DTU. This is the pipeline
we used for the final project, rewritten as a Snakefile so it runs over all six
samples instead of us typing every command per sample.

Samples: `1GC`, `3GC`, `4GC`, `6GC`, `8GC`, `9GC` (paired end, gastric cancer
RNA-seq from SRA).

The reads are mapped twice, once against human genome (GRCh38) and once
against a transcriptome reference, both go through the same variant calling
and filtering.

## What it does

```
raw fastq
  -> fastqc
  -> trimmomatic (adapter + quality trimming)
  -> fastqc again
  -> bwa mem (genome and transcriptome)
  -> samtools sort
  -> picard MarkDuplicates
  -> picard AddOrReplaceReadGroups
  -> samtools index
  -> gatk HaplotypeCaller (gvcf)
  -> tabix
  -> gatk GenotypeGVCFs
  -> gatk VariantFiltration (hard filtering)
  -> bedtools intersect with the mappability track (filter99)
  -> snpEff annotation
```

## Running it

Everything is hardcoded for the DTU pupil server, the tool paths are written
into the rules. To run somewhere else, these would have to be adjusted
accordingly (and the reference/dbsnp paths in the params).

```
snakemake -np
snakemake --jobs 1
```

The first one is a dry run to check what it would do, the second one actually
runs it. We only had one thread on the server so there is no point in giving it
more jobs.

Expected layout:

```
final_project/
  Snakefile
  raw_data/     <- subsampled fastq files here
  fastqc/
  trimmed/
  mapped/
  sorted/
  dedup/
  gvcf/
  vcf/
  hard_filtering/
  annotation/
```

## Things we did by hand and did not put in the Snakefile

Subsampling the downloaded SRA files to the first 1M reads, otherwise the whole
thing takes forever:

```
zcat SRR20074880_2.fastq.gz | head -n 4000000 | gzip > 1GC_sub1.fastq.gz
```

Copying files to the server from windows (wsl):

```
scp /mnt/c/Users/<you>/Desktop/8GC_sub2.fastq.gz <student-id>@pupil2.healthtech.dtu.dk:/home/projects/22126_NGS/projects/group14/final_project/raw_data
```

Making the project folder writable for everyone in the group:

```
chmod -R a+rw /home/projects/22126_NGS/projects/group14/final_project
find /home/projects/22126_NGS/projects/group14/final_project -type d -exec chmod 1777 {} \;
find /home/projects/22126_NGS/projects/group14/final_project -type f -exec chmod 666 {} \;
```

Counting how many sites did not pass the filters, and which filters they failed:

```
/home/ctools/bcftools-1.23/bcftools view -H 1GC_genome_filtering.vcf.gz | grep -v PASS | wc -l
/home/ctools/bcftools-1.23/bcftools view -H 1GC_genome_filtering.vcf.gz | grep -v PASS | cut -f7 | sort | uniq -c | sort -n
```

Counting how many variants are left after the mappability filter:

```
/home/ctools/bcftools-1.23/bcftools view -H 1GC_genome_filtering_map99.vcf.gz | grep PASS | wc -l
```

Counting annotations, see `count_annotations.sh`.

## Results we wrote down

Duplicates marked by picard:

| sample | transcriptome | genome |
| --- | --- | --- |
| 1GC | 5557 | 2450 |
| 3GC | 5462 | 2092 |
| 4GC | 5520 | - |
| 6GC | 5767 | - |
| 8GC | 18205 | 7123 |
| 9GC | 16928 | 6476 |

Sites that did not pass the hard filters:

| sample | filtered out |
| --- | --- |
| 1GC | 28875 |
| 3GC | 3261798 |
| 4GC | 16753 |
| 6GC | 13562 |
| 8GC | 2673671 |
| 9GC | 2093108 |

Variants left after also filtering on mappability (filter99):

| sample | PASS |
| --- | --- |
| 1GC | 312 |
| 3GC | 3588 |
| 4GC | 3698 |
| 6GC | 3204 |
| 8GC | 2804 |
| 9GC | 3509 |

Note that the numbers for 3GC, 8GC and 9GC are much bigger because those were
run with a QUAL30 filter and the others were not, we only noticed that later.

## Notes / known issues

- The transcriptome branch calls variants against the transcriptome fasta. In
  our original notes we had the genome fasta there by accident.
- `bwa mem` already pipes into `samtools sort`, so the extra `sort` rule is
  redundant. Left in because that is how the exercise did it.
- No conda envs, the tools are just whatever was installed on the server.
