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
        expand(outdir + "/{assembly_basename}.trinity.Trinity.fasta", assembly_basename=assembly_info.keys()),
        expand(outdir + "/raw_assemblies/{assembly_basename}.trinity.Trinity.fasta", assembly_basename=assembly_info.keys()),
        expand(outdir + "/{assembly_basename}.rnaspades/transcripts.fasta", assembly_basename=assembly_info.keys()),
        expand(outdir + "/raw_assemblies/{assembly_basename}.rnaspades/transcripts.fasta", assembly_basename=assembly_info.keys()),
        expand(outdir + "/{assembly_basename}.{assembler}/quant/{sample}_quant/quant.sf", assembly_basename=assembly_info.keys(), sample=SAMPLES, assembler = ['trinity']),
        #expand(outdir + "/{assembly_basename}/dammit/{assembly_basename}.dammit.gff3", assembly_basename=assembly_info.keys()),

wildcard_constraints:
    assembly_basename = '\w+'


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

rule annotate:
    message:
        "Annotate assembly"
    input:
        expand(outdir + "/{assembly_basename}/dammit/{assembly_basename}.dammit.gff3", assembly_basename=assembly_info.keys()),

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
        lambda w: ancient(sample_info.at[w.sample, f'read{w.end}']),
    output: outdir + "/fastqc/{sample}_{end}_fastqc.html",
    params:
        outdir=f"{outdir}/fastqc"
    conda: "conf/env/trim.yml"
    resources:
        mem_mb=3000,
        time=60,
        partition='low2',
    log: logdir + "/fastqc/{sample}_{end}.fastqc.log"
    benchmark: logdir + "/fastqc/{sample}_{end}.fastqc.benchmark.txt"
    shell:
        """
        fastqc {input} --outdir {params.outdir} 2> {log}
        """

rule fastqc_trimmed:
    input:
        outdir + '/trim/{sample}_{end}.trim.fq.gz',
    output:
        outdir + "/fastqc/{sample}_{end}.trim_fastqc.html",
    params:
        outdir=f"{outdir}/fastqc"
    conda: "conf/env/trim.yml"
    resources:
        mem_mb=3000,
        time=60,
        partition='low2',
    log: logdir + "/fastqc/{sample}_{end}.trim.fastqc.log"
    benchmark: logdir + "/fastqc/{sample}_{end}.trim.fastqc.benchmark.txt"
    shell:
        """
        fastqc {input} --outdir {params.outdir} 2> {log}
        """

rule multiqc:
    input: 
        raw=expand(outdir + "/fastqc/{sample}_{end}_fastqc.html", sample=SAMPLES, end = [1,2]),
        trimmed=expand(outdir + "/fastqc/{sample}_{end}.trim_fastqc.html", sample=SAMPLES, end = [1,2])
    output: outdir + "/fastqc/multiqc_report.html"
    params:
        fastqc_dir=f"{outdir}/fastqc",
    conda: "conf/env/trim.yml"
    resources:
        mem_mb=3000,
        time=60,
        partition='low2',
    log: logdir + "/multiqc/multiqc.log"
    benchmark: logdir + "/multiqc/multiqc.benchmark.txt"
    shell:
        """
        multiqc -f {params.fastqc_dir} --filename {output} 2> {log}
        """


localrules: write_trinity_samplesfile # fast/low mem, so run on current machine rather than submitting to cluster

rule write_trinity_samplesfile:
    message:
        "Writing samples file for trinity"
    input: lambda w: expand(outdir + "/trim/{sample}_{end}.trim.fq.gz", sample=assembly_info[w.assembly_basename], end=[1,2])
    output: outdir + "/{assembly_basename}.trinity_info.tsv"
    run:
        # default assembly info: all files together
        with open(str(output), 'w') as f:
            for inF in input:
                if '_1' in inF:
                    r1 = str(inF)
                    r2 = r1.replace("_1", "_2")
                    sample = os.path.basename(r1).rsplit("_", 1)[0] # requires last _ in filename to come after sample name
                    treatment = sample_info.at[sample, 'treatment']
                    section = sample_info.at[sample, 'section']
                    f.write(f"{treatment}\t{sample}\t{r1}\t{r2}\n")
                    
        #sample_info = outdir + "/{assembly_basename}.trinity_info.tsv",

rule trinity_assembly:
    input: 
        left = lambda w: expand(outdir + "/trim/{sample}_1.trim.fq.gz", sample=assembly_info[w.assembly_basename]),
        right = lambda w: expand(outdir + "/trim/{sample}_2.trim.fq.gz", sample=assembly_info[w.assembly_basename]),
    output:
        assembled_transcripts = outdir + "/{assembly_basename}.trinity.Trinity.fasta",
    resources:
        mem_mb= lambda wildcards, attempt: TRINITY_MAX_MEMORY * attempt,
        time=6000,
        partition='bmh',
    params:
        out_dir = lambda w: f"{outdir}/{w.assembly_basename}.trinity",
        mem_gb = str(int(TRINITY_MAX_MEMORY / 1000)) + 'G',
        left_str = lambda w: ','.join(expand(outdir + "/trim/{sample}_1.trim.fq.gz", sample=assembly_info[w.assembly_basename])),
        right_str = lambda w: ','.join(expand(outdir + "/trim/{sample}_2.trim.fq.gz", sample=assembly_info[w.assembly_basename])),
    log: logdir + "/trinity/{assembly_basename}.trinity.log"
    benchmark: logdir + "/trinity/{assembly_basename}.trinity.benchmark.txt"
    conda: "conf/env/trinity.yml"
    threads: 16
    shell:
        """
        Trinity --seqType fq --left {params.left_str} --right {params.right_str} \
                --CPU {threads} --max_memory {params.mem_gb} \
                --output {params.out_dir} > {log}
        """
        #Trinity --seqType fq --single reads.fq --max_memory 10G
        #Trinity --seqType fq --samples_file {input.sample_info} --CPU {threads} \
        #        --max_memory {resources.mem_mb}M --output {params.out_dir} 2> {log}


rule rnaspades_assembly:
    input:
        r1 = lambda w: expand(outdir + "/trim/{sample}_1.trim.fq.gz", sample=assembly_info[w.assembly_basename]),
        r2 = lambda w: expand(outdir + "/trim/{sample}_2.trim.fq.gz", sample=assembly_info[w.assembly_basename]),
    output:
        transcripts = outdir + "/{assembly_basename}.rnaspades/transcripts.fasta",    
    threads: 32
    params:
        spades_output_dir = lambda w: f"{outdir}/{w.assembly_basename}.rnaspades",
        rnaspades_tmp_dir = lambda w: f"{outdir}/{w.assembly_basename}.rnaspades/tmp",
        r1 = lambda w: ["-1 " + x for x in expand(outdir + "/trim/{sample}_1.trim.fq.gz ", sample=assembly_info[w.assembly_basename])],
        r2 = lambda w: ["-2 " + x for x in expand(outdir + "/trim/{sample}_2.trim.fq.gz ", sample=assembly_info[w.assembly_basename])],

    resources:
        mem_mb = 100000,
        nodes = 1,
        time = 600,
        partition = "bmh"
    conda: "conf/env/spades.yml"
    log: logdir + "/rnaspades/{assembly_basename}.log"
    benchmark: logdir + "/rnaspades/{assembly_basename}.benchmark.txt"
    shell:
        """
        rnaspades.py {params.r1} {params.r2} -t {threads} \
                     --tmp-dir {params.rnaspades_tmp_dir} \
                     -o {params.spades_output_dir} 2> {log}
        """

rule raw_trinity_assembly:
    input:
        left = lambda w: expand("inputs/{sample}_1.fastq.gz", sample=assembly_info[w.assembly_basename]),
        right = lambda w: expand("inputs/{sample}_2.fastq.gz", sample=assembly_info[w.assembly_basename]),
        #left = lambda w: ancient(expand(sample_info.at["{sample}", 'read1'], sample=assembly_info[w.assembly_basename])),
        #right = lambda w: ancient(expand(sample_info.at["{sample}", 'read2'], sample=assembly_info[w.assembly_basename])),
    output:
        assembled_transcripts = outdir + "/raw_assemblies/{assembly_basename}.trinity.Trinity.fasta",
    resources:
        mem_mb= lambda wildcards, attempt: TRINITY_MAX_MEMORY * attempt,
        time=6000,
        partition='bmh',
    params:
        out_dir = lambda w: f"{outdir}/raw_assemblies/{w.assembly_basename}.trinity",
        mem_gb = str(int(TRINITY_MAX_MEMORY / 1000)) + 'G',
        left_str = lambda w: ','.join(expand("inputs/{sample}_1.fastq.gz", sample=assembly_info[w.assembly_basename])),
        right_str = lambda w: ','.join(expand("inputs/{sample}_2.fastq.gz", sample=assembly_info[w.assembly_basename])),
    log: logdir + "/trinity/raw.{assembly_basename}.trinity.log"
    benchmark: logdir + "/trinity/raw.{assembly_basename}.trinity.benchmark.txt"
    conda: "conf/env/trinity.yml"
    threads: 16
    shell:
        """
        Trinity --seqType fq --left {params.left_str} --right {params.right_str} \
                --CPU {threads} --max_memory {params.mem_gb} --trimmomatic \
                --output {params.out_dir} > {log}
        """

rule raw_rnaspades_assembly:
    input:
        r1 = lambda w: expand("inputs/{sample}_1.fastq.gz", sample=assembly_info[w.assembly_basename]),
        r2 = lambda w: expand("inputs/{sample}_2.fastq.gz", sample=assembly_info[w.assembly_basename]),
        #r1 = lambda w: ancient(expand(sample_info.at["{sample}", 'read1'], sample=assembly_info[w.assembly_basename])),
        #r2 = lambda w: ancient(expand(sample_info.at["{sample}", 'read2'], sample=assembly_info[w.assembly_basename])),
    output:
        transcripts = outdir + "/raw_assemblies/{assembly_basename}.rnaspades/transcripts.fasta",
    threads: 32
    params:
        spades_output_dir = lambda w: f"{outdir}/raw_assemblies/{w.assembly_basename}.rnaspades",
        rnaspades_tmp_dir = lambda w: f"{outdir}/raw_assemblies/{w.assembly_basename}.rnaspades/tmp",
        r1 = lambda w: ["-1 " + x for x in expand("inputs/{sample}_1.fastq.gz", sample=assembly_info[w.assembly_basename])],
        r2 = lambda w: ["-2 " + x for x in expand("inputs/{sample}_2.fastq.gz", sample=assembly_info[w.assembly_basename])],
    resources:
        mem_mb = 100000,
        nodes = 1,
        time = 600,
        partition = "bmh"
    conda: "conf/env/spades.yml"
    log: logdir + "/rnaspades/raw.{assembly_basename}.log"
    benchmark: logdir + "/rnaspades/raw.{assembly_basename}.benchmark.txt"
    shell:
        """
        rnaspades.py {params.r1} {params.r2} -t {threads} \
                     --tmp-dir {params.rnaspades_tmp_dir} \
                     -o {params.spades_output_dir} 2> {log}
        """

#rule dammit_annotate:
#    input:
#        assembly = outdir + "trinity/{assembly_basename}.trinity.fasta",
#    output:
#        dammit_gff3 = outdir + "/{assembly_basename}/dammit/{assembly_basename}.dammit.gff3",
#        dammit_fasta = outdir + "/{assembly_basename}/dammit/{assembly_basename}.dammit.fasta",
#        dammit_report = outdir + "{assembly_basename}/dammit/{assembly_basename}.dammit.report.html",
#    params:
#        max_memory=TRINITY_MAX_MEMORY,
#        outdir= lambda w: f"{outdir}/{w.assembly_basename}/dammit",
#        database_dir='databases',
#    threads: 20
#    resources:
#        mem_mb=max(150000, TRINITY_MAX_MEMORY),
#        time=6000,
#        partition='bmh',
#    log: logdir + "/dammit/{assembly_basename}.dammit.log"
#    benchmark: logdir + "/dammit/{assembly_basename}.dammit.benchmark.txt"
#    conda: "conf/env/dammit.yml"
#    shell:
#        """
#        dammit annotate {input.assembly} --busco-group eukaryota --n_threads {threads} \
#               --database-dir {params.database_dir} --out-dir {params.outdir} \
#               --full -o {wildcards.assembly_basename} --no-rename \
#               --profile slurm 2> {log}
#        """


### index assembly and quantify reads with salmon ###
rule salmon_index:
    message:
        "Indexing transcripts with salmon"
    input: outdir + "/{assembly_name}.Trinity.fasta" 
    output: directory(outdir + "/{assembly_name}/quant/{assembly_name}_index")
    conda: "conf/env/salmon.yml"
    resources:
        mem_mb=10000,
        time=600,
        partition='med2',
    log: logdir + "/salmon/{assembly_name}.salmon_index.log"
    benchmark: logdir + "/salmon/{assembly_name}.salmon_index.benchmark"
    shell:
        """
        salmon index --index {output} --transcripts {input} -k 31 2> {log}
        """

rule salmon_quantify:
    input:
        r1 = outdir + '/trim/{sample}_1.trim.fq.gz',
        r2 = outdir + '/trim/{sample}_2.trim.fq.gz',
        index_dir=outdir + "/{assembly_name}/quant/{assembly_name}_index",
    output:
        outdir + "/{assembly_name}/quant/{sample}_quant/quant.sf"
    params:
        outdir= lambda w: f"{outdir}/{w.assembly_name}/quant/" + w.sample + "_quant"
    conda: "conf/env/salmon.yml"
    resources:
        mem_mb=10000,
        time=600,
        partition='med2',
    log: logdir + "/salmon/{assembly_name}/{sample}.quant.log"
    benchmark: logdir + "/salmon/{assembly_name}/{sample}.quant.benchmark"
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

