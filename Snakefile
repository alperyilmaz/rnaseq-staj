import pandas as pd


sample_table    = pd.read_table('annotation', sep='\t', lineterminator='\n')

sample_table.rename( columns={'Unnamed: 0':'SampleName'}, inplace=True )
sample_table    = sample_table.drop_duplicates(subset='SampleName', keep='first', inplace=False)
sample_table    = sample_table.dropna()
sample_table.set_index('SampleName',inplace=True)
SAMPLES=sample_table.index.values

rule all:
    input:
        "heatmap2.png",
        expand("{sample}_genecounts",sample=SAMPLES),
        expand("{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)

rule prepare_genome:
    output: 
        "dmel_6_star_index/SAindex"
    conda:
        "envs/star.yaml"
    threads: 16
    shell:
        """
        mkdir -p dmel_6_star_index
        wget http://ftp.flybase.net/releases/FB2022_01/dmel_r6.44/fasta/dmel-all-chromosome-r6.44.fasta.gz
        gunzip dmel-all-chromosome-r6.44.fasta.gz
        wget http://ftp.flybase.net/releases/FB2022_01/dmel_r6.44/gtf/dmel-all-r6.44.gtf.gz
        gunzip dmel-all-r6.44.gtf.gz
        head -1719359 dmel-all-chromosome-r6.44.fasta >| dmel-trimmed.fasta
        awk 'length($1) < 3' dmel-all-r6.44.gtf > dmel-trimmed.gtf
        STAR --runThreadN 16 \
          --runMode genomeGenerate \
          --genomeDir dmel_6_star_index \
          --genomeFastaFiles dmel-trimmed.fasta \
          --sjdbGTFfile dmel-trimmed.gtf \
          --genomeSAindexNbases 12 \
          --sjdbOverhang 149
        """

rule download:
    output:
       "fastqdump/{sample}_1.fastq",
       "fastqdump/{sample}_2.fastq"
    conda:
       "envs/sratools.yaml"
    threads: 6
    shell:
       "fasterq-dump -e {threads} -O fastqdump/ --split-files {wildcards.sample}"

rule fastqc:
     input: 
         "fastqdump/{sample}_1.fastq"
     output:
         "fastqdump/{sample}_1_fastqc.html"
     conda:
         "envs/fastqc.yaml"
     shell:
         "fastqc {input}"

rule starmap:
     input:
         r1="fastqdump/{sample}_1.fastq",
         r2="fastqdump/{sample}_2.fastq",
         index="dmel_6_star_index/SAindex"
     output:
         alignment="{sample}_Aligned.sortedByCoord.out.bam",
         genecount="{sample}_ReadsPerGene.out.tab"
     threads: 32
     conda:
         "envs/star.yaml"
     shell:
         """
         STAR --genomeDir dmel_6_star_index \
           --runThreadN {threads} \
           --readFilesIn {input.r1} {input.r2} \
           --outFileNamePrefix {wildcards.sample}_ \
           --outSAMtype BAM SortedByCoordinate \
           --outSAMunmapped Within \
           --outSAMattributes Standard \
           --quantMode GeneCounts \
           #--readFilesCommand zcat
         """

rule extractcounts:
     input:
         "{sample}_ReadsPerGene.out.tab"
     output:
         "{sample}_genecounts"
     shell:
         """
         awk -vSRR={wildcards.sample} 'BEGIN{{print "gene",SRR}} NR>4 {{print $1,$2}}' {input} > {output}
         """

rule mergedcounts:
     input:
         counts=expand("{sample}_genecounts",sample=SAMPLES)
     output:
         "merged_counts"
     shell:
        """
        paste {input.counts} | tr " " "\t" | awk '{{ printf"%s\t",$1; for (i=2;i<=NF;i+=2){{printf"%s\t",$i}}; print"" }}' | sed  's/\s*$//' > {output}
        """

rule deseq:
    input:
        data="merged_counts",
        script="scripts/dmelrna.r",
        anno="annotation"
    output:
        "deseq_results"
    conda:
        "envs/deseq.yaml"
    shell:
        "Rscript {input.script}"

rule heatmap:
    input:
        data="deseq_results",
        script="scripts/draw_heatmap.R"
    output:
        "heatmap2.png"
    conda:
        "envs/deseq.yaml"
    shell:
        "Rscript {input.script}"
