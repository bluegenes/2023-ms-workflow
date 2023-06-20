import os
import pandas as pd

# read config file
configfile: "inputs/config.yml"

# get/set output directory
outdir= config.get('output_dir', 'output.ms')
logdir = os.path.join(outdir, 'logs') 

# get and prep sample info
sample_csv = config['samples']
sample_info = pd.read_csv(sample_csv, sep = ',')
SAMPLES = sample_info['name'].tolist()
sample_info.set_index('name', inplace=True)

# which samples to use for assembly?
# by default, use all. 
assembly_default = {"ms": SAMPLES}
assembly_info = config.get('assembly_info', assembly_default)

# do we also want to assemble all samples together?
if config.get("also_assemble_all", False):
    assembly_info.update(assembly_default)

# get/set some useful params
TRINITY_MAX_MEMORY = config.get('trinity_max_memory', 150000)

#### define workflow endpoints (targets) ####
rule all:
    message:
        "RNA-Seq workflow: QC Reads, Assemble, and Quantify"
    input:
        outdir + "/fastqc/multiqc_report.html",
        expand(outdir + "trinity/{assembly_basename}_assembly.fasta", assembly_basename=assembly_info.keys()),
        expand(outdir + "/{assembly_basename}/quant/{sample}_quant/quant.sf", assembly_basename=assembly_info.keys(), sample=SAMPLES),

rule download_reads:
    message:
        "Download all reads"
    input:
        expand("inputs/{sample}_{end}.fastq.gz", sample=SAMPLES, end = [1,2]),

rule qc_reads:
    message:
        "QC reads"
    input:
        outdir + "/fastqc/multiqc_report.html",

rule assemble:
    message:
        "Assemble reads"
    input:
        expand(outdir + "trinity/{assembly_basename}_assembly.fasta", assembly_basename=assembly_info.keys()),

rule quantify:
    message:
        "Quantify reads"
    input:
        expand(outdir + "/{assembly_basename}/quant/{sample}_quant/quant.sf", assembly_basename=assembly_info.keys(), sample=SAMPLES),


# #### actual workflow rules ####

rule download_readfile:
    message:
        "download each read file to 'inputs' folder"
    output:
        r1=protected("inputs/{sample}_{end}.fastq.gz"),
    params:
        outdir = "inputs",
        url = lambda w: sample_info.at[w.sample, f'link{w.end}'],
    shell:
        """
        curl -O -J -L {params.url} --output-dir {params.outdir} # keeps name from download. errors here mean I got the download link wrong
        """
        # wget -O {output.r1} {params.url}


# Preprocess Reads #
rule trim_for_adapters_and_quality:
    message:
        "Trimming adapters and very low quality bases with fastp"
    input:
        r1 = lambda w: ancient(sample_info.at[w.sample, 'read1']),
        r2 = lambda w: ancient(sample_info.at[w.sample, 'read2']),
    output:
        r1 = outdir + '/trim/{sample}_1.trim.fq.gz',
        r2 = outdir + '/trim/{sample}_2.trim.fq.gz',
        json=outdir + "/trim/{sample}.trim.json",
        html=outdir + "/trim/{sample}.trim.html",
    conda: "conf/env/trim.yml"
    threads: 4
    resources:
        mem_mb=10000,
        time=600,
        partition='med2',
    log: logdir + "/trim/{sample}.trim.log"
    benchmark: logdir + "/trim/{sample}.trim.benchmark.txt"
    shell: """
        fastp --in1 {input.r1} --in2 {input.r2} \
             --detect_adapter_for_pe  --qualified_quality_phred 4 \
             --length_required 25 --correction --thread {threads} \
             --json {output.json} --html {output.html} \
             --low_complexity_filter --out1 {output.r1} --out2 {output.r2} 2> {log} 
    """

rule fastqc_raw:
    message:
        "Running fastqc on raw data"
    input:
        r1 = lambda w: ancient(sample_info.at[w.sample, 'read1']),
        r2 = lambda w: ancient(sample_info.at[w.sample, 'read2']),
    output: outdir + "/fastqc/{sample}.fastqc.html",
    params:
        outdir="rnaseq/raw_data/fastqc"
    conda: "conf/env/trim.yml"
    resources:
        mem_mb=10000,
        time=600,
        partition='med2',
    shell:
        """
        fastqc {input} --outdir {params.outdir}
        """

rule fastqc_trimmed:
    input:
        r1 = outdir + '/trim/{sample}_1.trim.fq.gz',
        r2 = outdir + '/trim/{sample}_2.trim.fq.gz',
    output:
        outdir + "/fastqc/{sample}.trim.fastqc.html",
    params:
        outdir= outdir + "fastqc"
    conda: "conf/env/trim.yml"
    resources:
        mem_mb=10000,
        time=600,
        partition='med2',
    shell:
        """
        fastqc {input} --outdir {params.outdir}
        """

rule multiqc:
    input: 
        raw=expand(outdir + "/fastqc/{sample}.fastqc.html", sample=SAMPLES),
        trimmed=expand(outdir + "/fastqc/{sample}.trim.fastqc.html", sample=SAMPLES)
    output: outdir + "/fastqc/multiqc_report.html"
    params:
        fastqc_dir=f"{outdir}/fastqc",
    conda: "conf/env/trim.yml"
    resources:
        mem_mb=10000,
        time=600,
        partition='med2',
    shell:
        """
        multiqc -f {params.fastqc_dir} --filename {output}
        """


localrules: write_trinity_samplesfile # fast/low mem, so run on current machine rather than submitting to cluster

rule write_trinity_samplesfile:
    message:
        "Writing samples file for trinity"
    input: lambda w: expand(outdir + "/trim/{sample}_{end}.trim.fq.gz", sample=assembly_info[w.assembly_basename], end=[1,2])
    output: outdir + "/{assembly_basename}.trinity_info.tsv"
    run:
        # default assembly info: all files together
        with open(output, 'w') as f:
            for inF in input:
                if '_1' in inF:
                    r1 = inF
                    r2 = inF.replace("_1", "_2")
                    sample = os.path.basename(inF).rsplit("_", 1)[0] # requires last _ in filename to come after sample name
                    f.write(f"{sample}\t{r1}\t{r2}\n")
                    

rule trinity_assembly:
    input: 
        sample_info = outdir + "/{assembly_basename}.trinity_info.tsv" 
    output:
        assembled_transcripts = outdir + "trinity/{assembly_basename}_assembly.fasta",
    params:
        max_memory=TRINITY_MAX_MEMORY,
    threads: 20
    resources:
        mem_mb=max(150000, TRINITY_MAX_MEMORY),
        time=6000,
        partition='bmh',
    conda: "conf/env/trinity.yml"
    shell:
        "Trinity --seqType fq --samples_file {input.sample_info} --CPU {threads} --output {output.assembled_transcripts}"


### index assembly and quantify reads with salmon ###
rule salmon_index:
    message:
        "Indexing transcripts with salmon"
    input: outdir + "trinity/{assembly_basename}_assembly.fasta" 
    output: directory(outdir + "/{assembly_basename}/quant/{assembly_basename}_index")
    conda: "conf/env/salmon.yml"
    resources:
        mem_mb=10000,
        time=600,
        partition='med2',
    log: logdir + "/{assembly_basename}/salmon_index.log"
    benchmark: logdir + "/{assembly_basename}/salmon_index.benchmark"
    shell:
        """
        salmon index --index {output} --transcripts {input} -k 31 2> {log}
        """

rule salmon_quantify:
    input:
        r1 = outdir + '/trim/{sample}_1.trim.fq.gz',
        r2 = outdir + '/trim/{sample}_2.trim.fq.gz',
        index_dir=outdir + "/{assembly_basename}/quant/{assembly_basename}_index",
    output: outdir + "/{assembly_basename}/quant/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "rnaseq/quant/" + wildcards.sample + "_quant"
    conda: "conf/env/salmon.yml"
    resources:
        mem_mb=10000,
        time=600,
        partition='med2',
    log: logdir + "/{assembly_basename}/{sample}.salmon.log"
    benchmark: logdir + "/{assembly_basename}/{sample}.salmon.benchmark"
    shell:
        """
        salmon quant -i {input.index_dir} --libType A -1 {input.r1} \
                     -2 {input.r2} -o {params.outdir} --seqBias \
                     --gcBias --validateMappings 2> {log}
        """




# Differential Expression: DESeq2
# rule knit:
#     message:
#         "Differential Expression Analysis using DESeq2"
#     input:
#         "rnaseq-workflow.html",
#        # "rnaseq-workflow.pdf",
# rule knit_actual:
#     input:
#         # create a new filename for every entry in SAMPLES,
#         # replacing {name} with each entry.
#         expand("rnaseq/quant/{name}_quant/quant.sf", name=SAMPLES),
#         "rnaseq-workflow.Rmd",
#     output:
#         "rnaseq-workflow.{format}",
#     shell:
#         "./knit-Rmd.R {wildcards.format}_document"

