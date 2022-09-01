## Snakemake_workflow

### Summary 

Snakemake is a workflow engine that provides a readable Phython-based workflow definition language and a powooerful execution environment that scales from single-core workstations to compute clusters without modifying the workflow.

### Snakemake Language

The workflow is denifined in a `snakefile` through a domain-specific language that is close to standard Python syntax.

The workflow is implied by dependencies between the rules that arise from one rule needing an output file of another as an input file.

Rule definition specifies (i) a name, (ii) any number of input and output file and (iii)  either a shell command or Python code that creates the output from the input.

### Installing Snakemake

Snakemake can be installed using Bioconda with the shell command `conda install snakemake`

### Steps for installation

`conda install snakemake`
`conda activate snakemake` 

To visualize the process or what a snakefile rule does, we run `snakemake` with the `-np`

#### DAG
Also known as directed acrylic graphs

This defines the route map into which snakemake reasons.

The input of to `snakemake` are a destination to cater for the next move by the snakemake.

#### Working with dependencies

In most cases when working on this workflow your analysis such `fastqc` may require tools which aren't installed in your working directory.

This now allows you to create `conda` environment that will provide the software dependencies.

A file e.g `fastqc.yaml` is created containing the channels and dependencies necessary for acccessing the required process.

Using `conda env export` we now activate the file/env.

The file may look like this in example:

```
name: fastqc
channels:
  -bioconda
  -conda-forge
  -defaults
dependencies:
  -fastqc=0.11.7=5
  -libgcc-ng=7.2.0=hdf63c60_3
  -perl=5.26.2=h470a237_1
  -openjdk=8.0.152=h46b5887_1
  ```


### Creating a snakefile

A simple example of a snakefile should look like this
 
##### Example 1
 
```
rule complement:                            //this rule declares the start of a new execution called complement
  input:                                    //gives run for user to define the input file
     "data/text.txt"                        //this is the input file
  output:                                   //allows you to define your output file
     "complement_data/text.txt"             //the outputfile path saved
  shell:                                    //the shell command is defined here
      "cat{input}|tr atcg tagc > {output}"      //this command with take the input file, then pipe it to 'tr' for 'translate' by replacing atcg with tagc which will                                                  result to a complement saved into output 
```
##### Example 2

Using wildcards

```
samples = ["A", "B", "C"]

rule all:
  input:
    "results/call/all.vcf",
    "results/plots/quals.svg"
    
rule map_reads:   //we replace the concrete sample with a name 'A' with a wildcard {sample} then we create an output from this input and name it w.r.t 'sample'
  input:
    "data/genome.fa",
    "data/{sample}.fastq"
  output:
  "results/{sample}A.bam"
  conda:
    "envs/mapping.yaml"
  shell:
  "bwa mem {input} | samtools view -b - >{output}"
 
 rule sort_alignments:     //this rule sorts the obtained '.bam' file by genomic cordinate
  input:
    "results/mapped/{samples}.bam"
  output:
    "results/{sample}.sorted.bam"
  conda:
    "envs/mapping.yaml"
  shell:
    "samtools sort -o {output}{input}"
    
  rule call_variants:    //we aggregate over all samples perform a joint calling of genomic variants
    input:
      fa="data/genome.fa"
      bam=expand("results/mapped/{sample}.sorted.bam",sample=SAMPLES) //joint calling is performed by `expand` command
 
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
        "notebooks/plot-quals.py.ipynb" //instead of shell command , we use jupyter notebook intergration
      
```

The above conda: envs/mapping.yaml contains the following tools/dependencies; and also the channels
```

channels:
  -bioconda
  -conda-forge
dependencies:
  -bwa =0.7.17
  -samtools =1.9
  //upon execution, snakemake will create that environment.
  ```
  
  
