## Snakemake_workflow

### Summary 

Snakemake is a workflow engine that provides a readable Phython-based workflow definition language and a powooerful execution environment that scales from single-core workstations to compute clusters without modifying the workflow.

### Snakemake Language

The workflow is denifined in a 'snakefile' through a domain-specific language that is close to standard Python syntax.

The workflow is implied by dependencies between the rules that arise from one rule needing an output file of another as an input file.

Rule definition specifies (i) a name, (ii) any number of input and output file and (iii)  either a shell command or Python code that creates the output from the input.

### Installing Snakemake

Snakemake can be installed using Bioconda with the shell command `conda install snakemake`

### Steps for installation

`conda install snakemake`
`conda activate snakemake` 

### Creating a snakefile
```rule all:
  input:
    "results/call/all.vcf",
    "results/plots/quals.svg"
    
rule map_reads:
  input:
    "data/genome.fa",
    "data/{sample}.fastq"
  output:
  "results/{sample}A.bam"
  conda:
    "envs/mapping.yaml"
  shell:
  "bwa mem {input} | samtools view -b - >{output}"
 
 rule sort_alignments:
  input:
    "results/mapped/{samples}.bam"
  output:
    "results/{sample}.sorted.bam"
  conda:
    "envs/mapping.yaml"
  shell:
    "samtools sort -o {output}{input}"
    
  rule call_variants:
    input:
      fa="data/genome.fa"
      bam=expand("results/mapped/{sample}.sorted.bam",sample=SAMPLES)
   output:
      "results/calls/all.vcf"
    conda:
      "envs/calling.yaml"
    shell:
      "bcftools mpileup -f {input.fa}{output}.bam | bcftools call -mv - >{output}"
      
    rule plot_quals:
      input:
        "results/calls/all.vcf"
      output:
        "results/plots/quals.svg"
      conda:
        "envs/stats.yaml"
      notebook:
        "notebooks/plot-quals.py.ipynb"
      
```
