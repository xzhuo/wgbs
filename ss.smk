import re
import socket
from datetime import datetime
import fnmatch
import gzip
from math import floor

FASTQ = "fastq"
TMP_DIR = "tmp"
QC_DIR = "0_fastqc"
TRIM_DIR = "1_trim"
PHIX_DIR = "1_phix"
BISMARK_LAMBDA = "2_bismark_lambda"
BISMARK = "2_bismark"
PRESEQ = "3_preseq"
INSERT_CPG_BIAS = "4_insert_cpg_bias"
COVERAGE = "5_coverage"
TRACKS = "6_tracks"

threads_fastqc = 16
threads_trim = 16
threads_phix = 4
threads_lambda = 4
threads_bismark = 12
threads_coverage = 4
threads_preseq = 4

mem_bismark = "70G"

hostname = socket.gethostname()
if len(re.findall("compute1", hostname)) > 0:
    server = "ris"
elif len(re.findall("^n\\d+$", hostname)) == 0:
    server = "wanglab"
else:
    server = "htcf"

print(server)
configfile: "config.yaml"
genome = config["genome"]
if "params" in config.keys():
    include: config["params"]

if server == "ris":
    PHIX_REF = "/storage1/fs1/hprc/Active/xzhuo/genomes/phiX174/bwa_index/phiX174.fa"
    LAMBDA_DIR = "/storage1/fs1/hprc/Active/xzhuo/genomes/lambda"
    REF_DIR = "/storage1/fs1/hprc/Active/xzhuo/genomes/" + genome
    pipe_path = "/storage1/fs1/hprc/Active/xzhuo/github/wgbs"
elif server == "htcf":
    PHIX_REF = "/scratch/twlab/fanc/genomes/phiX174/bwa_index/phiX174.fa"
    LAMBDA_DIR = "/scratch/twlab/fanc/genomes/lambda"
    REF_DIR = "/scratch/twlab/fanc/genomes/" + genome + "/bismark"
    shell.prefix("module load cutadapt/1.16-python3 trimgalore/0.6.0-python-3.6.5-java-11 samtools/1.9 bwa/0.7.15 bowtie2/2.3.4.1 preseq/2.0.2 bedtools/2.27.1 htslib/1.3.1 r/3.6.3-python-3.6.5-java-11;")
    pipe_path = "/home/fanc/software/wgbs"
else:
    PHIX_REF = "/scratch/genomes/phiX174/bwa_index/phiX174.fa"
    LAMBDA_DIR = "/scratch/genomes/lambda"
    REF_DIR = "/bar/cfan/genomes/" + genome + "/bismark"
    shell.prefix("module load FastQC/0.11.5 multiqc/1.7 trim_galore/0.6.6 samtools/1.9 bwa/0.7.15 bismark/0.23.1 preseq/3.1.2 bedtools/2.27.1 htslib/1.3.1 R/3.6.1;")
    pipe_path = "/bar/cfan/software/wgbs"

trash = os.system("rm -rf scripts")
trash = os.system("ln -s " + pipe_path + "/scripts ./")

SUFFIX = '.fastq.gz'
suffix_length = len(SUFFIX) + 3
# SAMPLES = set(map(lambda x: x[:-suffix_length], filter(lambda y: y.endswith(SUFFIX), os.listdir("."))))
SAMPLES = set(map(lambda x: x[:-suffix_length], fnmatch.filter(os.listdir(FASTQ), '*.fastq.gz')))
print(SAMPLES)
rule mode_full:
    input:
        expand(QC_DIR + "/{sample}_R{n}_fastqc.zip", sample=SAMPLES, n=[1, 2]),
        expand(TRIM_DIR + "/{sample}_R{n}{ext}", sample=SAMPLES, n=[1, 2], ext=[".fq.gz", "_trimmed_fastqc.html", "_trimming_report.txt"]),
        expand(TRIM_DIR + "/{sample}_R{Rn}_val_{Rn}_fastqc.zip", sample=SAMPLES, Rn=[1, 2]),
        expand(PHIX_DIR + "/{sample}.bwa_phiX.{ext}", sample=SAMPLES, ext=["bam", "txt"]),
        expand(BISMARK_LAMBDA + "/{sample}_bismark_bt2.CXme.txt", sample=SAMPLES),
        expand(BISMARK + "/{sample}_bismark_bt2{ext}", sample=SAMPLES, ext=[".CXme.txt", "_pe.bam"]),
        expand(BISMARK + "/{sample}_bismark_bt2_PE_report.html", sample=SAMPLES),
        expand(PRESEQ + "/{sample}.preseq_lc_extrap.txt", sample=SAMPLES, allow_missing=True),
        expand(INSERT_CPG_BIAS + "/{sample}.insert_length.txt", sample=SAMPLES),
        expand(COVERAGE + "/{sample}.genome_cov.txt", sample=SAMPLES),
        expand(TRACKS + "/{sample}.{ext}", sample=SAMPLES, ext=["cov.bg.gz", "CG.methylC.gz"])
    # output:
    #     "{sample}_mode_full.txt"
    # shell:
    #     "echo 'full run' > {output}"

rule mode_shallow:
    input:
        expand(QC_DIR + "/{sample}_R{n}_fastqc.zip", n=[1, 2], allow_missing=True),
        expand(TRIM_DIR + "/{sample}_R{n}{ext}", n=[1, 2], ext=[".fq.gz", "_trimmed_fastqc.html", "_trimming_report.txt"], allow_missing=True),
        expand(TRIM_DIR + "/{sample}_R{Rn}_val_{Rn}_fastqc.zip", Rn=[1, 2], allow_missing=True),
        expand(PHIX_DIR + "/{sample}.bwa_phiX.{ext}", ext=["bam", "txt"], allow_missing=True),
        BISMARK_LAMBDA + "/{sample}_bismark_bt2.CXme.txt",
        BISMARK + "/{sample}_bismark_bt2_PE_report.html",
        INSERT_CPG_BIAS + "/{sample}.insert_length.txt",
    output:
        "{sample}_mode_shallow.txt"
    shell:
        "echo 'shallow run' > {output}"

rule del_all:
    output:
        "del_all.out"
    params:
        TMP_DIR,
        QC_DIR,
        TRIM_DIR,
        PHIX_DIR,
        BISMARK_LAMBDA,
        BISMARK,
        PRESEQ,
        INSERT_CPG_BIAS,
        COVERAGE,
        TRACKS,   
    shell:
        "rm -rf {params} && "
        "echo 'deleted all directories' > {output}"

rule fastqc:
    input:
        expand(FASTQ + "/{sample}_R{n}" + SUFFIX, n=[1, 2], allow_missing=True)
    threads:
        threads_fastqc
    params:
        dir = directory(QC_DIR)
    output:
        expand(QC_DIR + "/{sample}_R{n}_fastqc.zip", n=[1, 2], allow_missing=True)
    log:
        QC_DIR + "/{sample}.fastqc.log"
    shell:
        "fastqc -o {params.dir} --noextract --nogroup -t {threads} {input[0]} {input[1]} &> {log}"

rule trim:
    input:
        expand(FASTQ + "/{sample}_R{n}" + SUFFIX, n=[1, 2], allow_missing=True)
    output:
        trim_fq = expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True),
        html = expand(TRIM_DIR + "/{sample}_R{Rn}_trimmed_fastqc.html", Rn=[1, 2], allow_missing=True),
        zips = expand(TRIM_DIR + "/{sample}_R{Rn}_val_{Rn}_fastqc.zip", Rn=[1, 2], allow_missing=True),
        report = expand(TRIM_DIR + "/{sample}_R{Rn}_trimming_report.txt", Rn=[1, 2], allow_missing=True)
    log:
        TRIM_DIR + "/{sample}.trim_galore.log"
    threads:
        threads_trim
    params:
        dir = directory(TRIM_DIR),
        trim_base1 = 10,
        trim_base2 = 15,
        tmp = expand(TRIM_DIR + "/{sample}_R{rn}_val_{vn}.fq.gz", zip, rn=[1, 2], vn=[1, 2], allow_missing=True),
        # tmp_r2 = TRIM_DIR + "/{sample}_R2_val_2.fq.gz",
        html = TRIM_DIR + "/{sample}_R[12]_val_[12]_fastqc.html",
        report = expand(TRIM_DIR + "/{sample}_R{Rn}" + SUFFIX + "_trimming_report.txt", Rn=[1, 2], allow_missing=True)

    run:
        shell("""trim_galore -q 20 --phred33 --fastqc --fastqc_args "-o {params.dir} --noextract --nogroup" \
            --illumina --stringency 1 -e 0.1 --length 20 \
            --clip_R1 {params.trim_base1} --clip_R2 {params.trim_base2} \
            -o {params.dir} \
            -j {threads} \
            --paired --retain_unpaired -r1 21 -r2 21 {input[0]} {input[1]} &>{log}""")
        os.rename(params.tmp[0], output.trim_fq[0])
        os.rename(params.tmp[1], output.trim_fq[1])
        # shell("""mv {params.tmp_r1} {output.trim1_fq} \
        # && mv {params.tmp_r2} {output.trim2_fq}""")
        shell("""rename "s/_val_[12]/_trimmed/g" {params.html} \
            && rename "s/{SUFFIX}//g" {params.report}""")

rule phix:
    input:
        ref = PHIX_REF,
        reads = expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True)
        # R2_fq = TRIM_DIR + "/{sample}_R2.fq.gz"
    params:
        dir = directory(PHIX_DIR),
    output:
        bam = PHIX_DIR + "/{sample}.bwa_phiX.bam",
        txt = PHIX_DIR + "/{sample}.bwa_phiX.txt"
    log:
        PHIX_DIR + "/{sample}.bwa_phiX.log"
    threads:
        threads_phix
    run:
        total = int(shell("""zcat {input.reads[0]} | grep -c "^@" """, read=True).rstrip())
        shell("bwa mem -t {threads} {input.ref} {input.reads[0]} {input.reads[1]} 2> {log} | samtools view -b -o {output.bam} -F 4 -@ {threads} -")
        phi = int(shell("samtools view -c -@ {threads} {output.bam}", read=True).rstrip())
        # phi_rate = shell("""awk -v c={phi} -v t={total} "BEGIN{{ printf(\\"%.7g\\", c/(t*2))}}" """, read=True).rstrip()
        # shell("""echo -e "{wildcards.sample}\t{total}\t{phi}\t{phi_rate}" > {output.txt}""")
        phi_rate = phi / (total * 2)
        with open(output.txt, "w") as f:
            # print("sample\ttotal\tphi\tphi_rate", file = f)
            print("{:s}\t{:d}\t{:d}\t{:.7g}".format(wildcards.sample, total, phi, phi_rate), file=f)

rule bismark_lambda:
    input:
        expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True),
        # R1_fq = TRIM_DIR + "/{sample}_R1.fq.gz",
        # R2_fq = TRIM_DIR + "/{sample}_R2.fq.gz"
    threads:
        threads_lambda
    params:
        ref_dir = directory(LAMBDA_DIR),
        min_insert = 0,
        max_insert = 2000,
        out_dir = directory(BISMARK_LAMBDA),
        report = "{sample}_bismark_bt2_PE_report.html",
        temp = BISMARK_LAMBDA + "/{sample}_R1_bismark_bt2_*",
        cg_pe = BISMARK_LAMBDA + "/CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        ch_pe = BISMARK_LAMBDA + "/Non_CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        merged = BISMARK_LAMBDA + "/{sample}_bismark_bt2.extracted.txt.gz",
        bedGraph = "{sample}.bedGraph.gz",
        cov = BISMARK_LAMBDA + "/{sample}.bismark.cov.gz",
        cx_report = "{sample}.CX_report.txt.gz",
        parallel = 2,
        bt_p = floor(threads_lambda / 2),
        tmp_dir = directory(TMP_DIR + "/bismark_lambda")

    log:
        bismark_pe = BISMARK_LAMBDA + "/{sample}.bismark.log",
        dedup_pe = BISMARK_LAMBDA + "/{sample}.deduplicate_bismark.log",
        methx_pe = BISMARK_LAMBDA + "/{sample}_pe.bismark_methylation_extractor.log",
        html = BISMARK_LAMBDA + "/{sample}_pe.bismark2report.log",
        bismark2bg = BISMARK_LAMBDA + "/{sample}.bismark2bedGraph.log",
        cov2c = BISMARK_LAMBDA + "/{sample}.coverage2cytosine.log"

    output:
        bam_pe = BISMARK_LAMBDA + "/{sample}_bismark_bt2_pe.bam",
        bam_dedup_pe = BISMARK_LAMBDA + "/{sample}_bismark_bt2_pe.deduplicated.bam",
        report = BISMARK_LAMBDA + "/{sample}_bismark_bt2_PE_report.html",
        alignment_report = BISMARK_LAMBDA + "/{sample}_bismark_bt2_PE_report.txt",
        dedup_report = BISMARK_LAMBDA + "/{sample}_bismark_bt2_pe.deduplication_report.txt",
        splitting_report = BISMARK_LAMBDA + "/{sample}_bismark_bt2_pe.deduplicated_splitting_report.txt",
        mbias_report = BISMARK_LAMBDA + "/{sample}_bismark_bt2_pe.deduplicated.M-bias.txt",
        nucleotide_report = BISMARK_LAMBDA + "/{sample}_bismark_bt2_pe.nucleotide_stats.txt",
        cx_me = BISMARK_LAMBDA + "/{sample}_bismark_bt2.CXme.txt",
        cx_report = BISMARK_LAMBDA + "/{sample}.CX_report.txt.gz"

    run:
        print("Started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        # Mapping with bismark/bowtie2
        # Note --bowtie2 and -p $nthreads are both SLOWER than single threaded bowtie1
        print("1. Mapping to reference with bismark/bowtie2... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("mkdir -p {params.tmp_dir}")
        shell("bismark -q -I {params.min_insert} -X {params.max_insert} --parallel {params.parallel} -p {params.bt_p} \
            --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 \
            -o {params.out_dir} --temp_dir {params.tmp_dir} --gzip --nucleotide_coverage \
            {params.ref_dir} -1 {input[0]} -2 {input[1]} &>{log.bismark_pe}")
        shell("""rename "s/_R1_bismark_bt2/_bismark_bt2/g" {params.temp}""")

        # Deduplicate reads
        print("-- 2. Deduplicating aligned reads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("deduplicate_bismark -p --output_dir {params.out_dir} --bam {output.bam_pe} &>{log.dedup_pe}")

        # Run methylation extractor for the sample
        print("-- 3. Analyse methylation in " + output.bam_dedup_pe + " using " + str(threads) + " threads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report \
            -o {params.out_dir} --gzip --parallel {threads} \
            {output.bam_dedup_pe} &>{log.methx_pe}")

        # Generate HTML Processing Report
        print("-- 4. Generate bismark HTML processing report file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark2report -o {params.report} --dir {params.out_dir} \
           --alignment_report {output.alignment_report} \
           --dedup_report {output.dedup_report} \
           --splitting_report {output.splitting_report} \
           --mbias_report {output.mbias_report} \
           --nucleotide_report {output.nucleotide_report} &>{log.html}")

        # Generate bedGraph file
        print("-- 5. Generate bedGraph file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("mv {params.ch_pe} {params.merged}")
        shell("cat {params.cg_pe} >>{params.merged}")

        shell("bismark2bedGraph --dir {params.out_dir} --cutoff 1 --CX_context --buffer_size=75G --scaffolds \
            -o {params.bedGraph} {params.merged} &>{log.bismark2bg}")
        shell("rm {params.out_dir}/{params.bedGraph}")  # $merged

        # Calculate average methylation levels per each CN context
        print("-- 6. Generate cytosine methylation file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("coverage2cytosine -o {wildcards.sample} --dir {params.out_dir} --genome_folder {params.ref_dir} --CX_context --gzip \
              {params.cov} &>{log.cov2c}")

        shell("rm {params.cov}")
        shell("""zcat {params.out_dir}/{params.cx_report} | \
            awk "BEGIN{{ca=0;cc=0;cg=0;ct=0;mca=0;mcc=0;mcg=0;mct=0}} \
            \$7~/^CA/ {{ca+=\$5; mca+=\$4}} \
            \$7~/^CC/ {{cc+=\$5; mcc+=\$4}} \
            \$7~/^CG/ {{cg+=\$5; mcg+=\$4}} \
            \$7~/^CT/ {{ct+=\$5; mct+=\$4}} \
            END{{printf(\\"CA\\t%d\\t%d\\t%.3f\\n\\", ca, mca, mca/(ca+mca)); \
            printf(\\"CC\\t%d\\t%d\\t%.3f\\n\\", cc, mcc, mcc/(cc+mcc)); \
            printf(\\"CG\\t%d\\t%d\\t%.3f\\n\\", cg, mcg, mcg/(cg+mcg)); \
            printf(\\"CT\\t%d\\t%d\\t%.3f\\n\\", ct, mct, mct/(ct+mct));}}" >{output.cx_me}""")

        # Print the files generated
        print("-- The results...")
        shell("ls -l {params.out_dir}/*{wildcards.sample}*")
        print("-- Finished on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

rule bismark_align:
    input:
        expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True)
    output:
        bam_pe = BISMARK + "/{sample}_bismark_bt2_pe.bam",
        alignment_report = BISMARK + "/{sample}_bismark_bt2_PE_report.txt",
        nucleotide_report = BISMARK + "/{sample}_bismark_bt2_pe.nucleotide_stats.txt"
    params:
        ref_dir = directory(REF_DIR),
        min_insert = 0,
        max_insert = 700,
        out_dir = directory(BISMARK),
        temp = BISMARK + "/{sample}_R1_bismark_bt2_*",
        parallel = 2,
        bt_p = floor(threads_bismark / 2),
        tmp_dir = directory(TMP_DIR + "/bismark")
    threads:
        threads_bismark
    log:
        bismark_pe = BISMARK + "/{sample}.bismark.log"
    run:
        print("Bismark: Started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        print("1. Mapping to reference with bismark/bowtie2... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("mkdir -p {params.tmp_dir}")
        shell("bismark -q -I {params.min_insert} -X {params.max_insert} --parallel {params.parallel} -p {params.bt_p} \
            --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 \
            -o {params.out_dir} --temp_dir {params.tmp_dir} --gzip --nucleotide_coverage \
            {params.ref_dir} -1 {input[0]} -2 {input[1]} &>{log.bismark_pe}")
        shell("""rename "s/_R1_bismark_bt2/_bismark_bt2/g" {params.temp}""")

rule bismark_dedup:
    input:
        bam_pe = BISMARK + "/{sample}_bismark_bt2_pe.bam"
    output:
        bam_dedup_pe = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.bam",
        dedup_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplication_report.txt",
    params:
        out_dir = directory(BISMARK),
    log:
        dedup_pe = BISMARK + "/{sample}.deduplicate_bismark.log",
    threads:
        threads_bismark
    run:
        print("-- 2. Deduplicating aligned reads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("deduplicate_bismark -p --output_dir {params.out_dir} --bam {input.bam_pe} &>{log.dedup_pe}")

rule bismark_extract:
    input:
        bam_dedup_pe = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.bam",
    output:
        cg_pe = BISMARK + "/CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        ch_pe = BISMARK + "/Non_CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        splitting_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated_splitting_report.txt",
        mbias_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.M-bias.txt",
    params:
        out_dir =  directory(BISMARK),
    log:
        methx_pe = BISMARK + "/{sample}_pe.bismark_methylation_extractor.log",
    threads:
        threads_bismark
    run:
        print("-- 3. Analyse methylation in " + input.bam_dedup_pe + " using " + str(threads) + " threads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report \
            -o {params.out_dir} --gzip --parallel {threads} \
            {input.bam_dedup_pe} &>{log.methx_pe}")

rule bismark_report:
    input:
        alignment_report = BISMARK + "/{sample}_bismark_bt2_PE_report.txt",
        dedup_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplication_report.txt",
        splitting_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated_splitting_report.txt",
        mbias_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.M-bias.txt",
        nucleotide_report = BISMARK + "/{sample}_bismark_bt2_pe.nucleotide_stats.txt"
    output:
        html = BISMARK + "/{sample}_bismark_bt2_PE_report.html"
    params:
        out_dir = directory(BISMARK),
        html = "{sample}_bismark_bt2_PE_report.html",
        # note bismark will automatically add directory name... that's why html is repeated in output and params
    log:
        html = BISMARK + "/{sample}_pe.bismark2report.log",
    threads:
        1
    run:
        print("-- 4. Generate bismark HTML processing report file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark2report -o {params.html} --dir {params.out_dir} \
           --alignment_report {input.alignment_report} \
           --dedup_report {input.dedup_report} \
           --splitting_report {input.splitting_report} \
           --mbias_report {input.mbias_report} \
           --nucleotide_report {input.nucleotide_report} &>{log.html}")

rule bismark_bdg:
    input:
        cg_pe = BISMARK + "/CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        ch_pe = BISMARK + "/Non_CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz"
    output:
        bedGraph = BISMARK + "/{sample}.bedGraph.gz",
        cov = BISMARK + "/{sample}.bismark.cov.gz",
    params:
        bedGraph = "{sample}.bedGraph.gz", # again, repeated this one from output because bismark require filename without dirname
        merged = BISMARK + "/{sample}_bismark_bt2.extracted.txt.gz",
        out_dir = directory(BISMARK),
        mem = mem_bismark,
    log:
        bismark2bg = BISMARK + "/{sample}.bismark2bedGraph.log",
    threads:
        threads_bismark
    run:
        print("-- 5. Generate bedGraph file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("cat {input.ch_pe} {input.cg_pe} >> {params.merged}")
        shell("bismark2bedGraph --dir {params.out_dir} --cutoff 1 --CX_context --buffer_size={params.mem} \
            -o {params.bedGraph} {params.merged} &>{log.bismark2bg}")
rule bismark_c2c:
    input:
        cov = BISMARK + "/{sample}.bismark.cov.gz",
    output:
        cx_report = BISMARK + "/{sample}.CX_report.txt.gz",
    params:
        out_dir = directory(BISMARK),
        ref_dir = directory(REF_DIR),
    log:
        cov2c = BISMARK + "/{sample}.coverage2cytosine.log"
    threads:
        threads_bismark
    run:
        print("-- 6. Generate cytosine methylation file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("coverage2cytosine -o {wildcards.sample} --dir {params.out_dir} --genome_folder {params.ref_dir} --CX_context --gzip \
              {input.cov} &>{log.cov2c}")
        
rule bismark_c2c_hj:
    input:
        cx_report = BISMARK + "/{sample}.CX_report.txt.gz"
    output:
        cx_me = BISMARK + "/{sample}_bismark_bt2.CXme.txt",
    params:
    log:
        BISMARK + "/{sample}_bismark_bt2.CXme.log",
    threads:
        1
    run:
        print("-- 6. Hyung Joo's script to reformat c2c results... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("""zcat {input.cx_report} | \
            awk "BEGIN{{ca=0;cc=0;cg=0;ct=0;mca=0;mcc=0;mcg=0;mct=0}} \
                \$7~/^CA/ {{ca+=\$5; mca+=\$4}} \
                \$7~/^CC/ {{cc+=\$5; mcc+=\$4}} \
                \$7~/^CG/ {{cg+=\$5; mcg+=\$4}} \
                \$7~/^CT/ {{ct+=\$5; mct+=\$4}} \
                END{{printf(\\"CA\\t%d\\t%d\\t%.3f\\n\\", ca, mca, mca/(ca+mca)); \
                printf(\\"CC\\t%d\\t%d\\t%.3f\\n\\", cc, mcc, mcc/(cc+mcc)); \
                printf(\\"CG\\t%d\\t%d\\t%.3f\\n\\", cg, mcg, mcg/(cg+mcg)); \
                printf(\\"CT\\t%d\\t%d\\t%.3f\\n\\", ct, mct, mct/(ct+mct));}}" >{output.cx_me} 2>{log}""")
        
        print("-- Finished on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

rule preseq:
    input:
        BISMARK + "/{sample}_bismark_bt2_pe.bam"
    params:
        dir = directory(PRESEQ),
        bam_sorted = PRESEQ + "/{sample}.sorted.bam",
        temptDir = TMP_DIR
    output:
        PRESEQ + "/{sample}.preseq_lc_extrap.txt"
    log:
        PRESEQ + "/{sample}.preseq_lc_extrap.log"
    threads:
        threads_preseq
    run:
        shell("samtools sort -m 2G -o {params.bam_sorted} -T {params.temptDir}/preseq_{wildcards.sample} -@ {threads} {input}")
        try:
            shell("preseq lc_extrap -o {output} -B -P {params.bam_sorted} 2>{log}")
        except subprocess.CalledProcessError:
            shell("preseq lc_extrap -o {output} -B -P -D {params.bam_sorted} 2>{log}")

rule track_coverage:
    input:
        bam_in = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.bam",
    params:
        dir = directory(TRACKS),
        temptDir = TMP_DIR
    threads:
        threads_coverage
    output:
        bam_sorted = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.sorted.bam",
        cov_out = TRACKS + "/{sample}.cov.bg.gz"
    run:
        shell("samtools sort -m 2G -o {output.bam_sorted} -T {params.temptDir}/track_{wildcards.sample} -@ {threads} {input.bam_in}")
        shell("samtools index -@ {threads} {output.bam_sorted}")
        shell("bedtools genomecov -bg -ibam {output.bam_sorted} | \
            bgzip > {output.cov_out} && tabix -p bed {output.cov_out}")

rule track_mergedCG:
    input:
        BISMARK + "/{sample}.CX_report.txt.gz"
    params:
        dir = directory(TRACKS)
    output:
        TRACKS + "/{sample}.CG.methylC.gz"
    shell:
        """zcat {input} | \
            awk -F"\\t" "BEGIN{{OFS=FS}} \$6==\\"CG\\" && \$4+\$5>0 {{ if (\$3==\\"+\\") {{print \$1,\$2-1,\$2+1,\$4,\$5}} if (\$3==\\"-\\") {{print \$1,\$2-2,\$2,\$4,\$5}} }}" | \
            sort -k1,1 -k2,2n | bedtools groupby -g 1,2,3 -c 4,5 -o sum,sum | \
            awk -F"\\t" "BEGIN{{OFS=FS}} {{mcg=sprintf(\\"%.3f\\", \$4/(\$4+\$5)); print \$1,\$2,\$3,\\"CG\\",mcg,\\"+\\",\$4+\$5 }}" | \
            bgzip > {output} && tabix -p bed {output}"""

rule insert_cpg_bias:
    input:
        bam_dedup = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.bam",
        bg_chr1_1kb = REF_DIR + "/CpGs.chr1.1kb_win.bg.gz",
        rscript_insert = "scripts/density_insert_length.R",
        rscript_cpgbias = "scripts/CpGbias_1kb.R"
    params:
        dir = directory(INSERT_CPG_BIAS),
        bam_tmp = INSERT_CPG_BIAS + "/{sample}.tmp.bam",
        insert_tmp = INSERT_CPG_BIAS + "/{sample}.tmp.insert.txt.gz",
        cov_chr1 = INSERT_CPG_BIAS + "/{sample}.tmp.CpG.cov_chr1_1kb_win.txt.gz",
        sam_threads = 5
    threads:
        1
    output:
        INSERT_CPG_BIAS + "/{sample}.insert_length.txt"
    log:
        insert_len = INSERT_CPG_BIAS + "/{sample}.density_insert_length.R.log",
        cpgbias = INSERT_CPG_BIAS + "/{sample}.CpGbias_1kb.R.log"
    run:
        shell("""bedtools bamtobed -bedpe -i {input.bam_dedup} | awk -v OFS="\t" "{{print \$1,\$2,\$6,\$6-\$2}}" | gzip -nc > {params.insert_tmp}""")

        # select only chr1 
        shell("""bedtools bamtobed -bedpe -i {input.bam_dedup} | awk "\$1==\\"chr1\\"" | \
            bedtools coverage -counts -a {input.bg_chr1_1kb} -b stdin | awk "\$4>0 && \$5>0" | gzip -nc > {params.cov_chr1}""")

        # run rscripts
        shell("Rscript {input.rscript_insert} {params.insert_tmp} {output} &> {log.insert_len}")
        shell("Rscript {input.rscript_cpgbias} {params.cov_chr1} {wildcards.sample} &> {log.cpgbias}")

        # remove temporary files
        shell("rm {params.insert_tmp} {params.cov_chr1}")

rule qc_cal_genome_cov:
    input:
        bismark_cx = BISMARK + "/{sample}_bismark_bt2.CXme.txt",
        dedup_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplication_report.txt"
    params:
        directory(COVERAGE)
    output:
        COVERAGE + "/{sample}.genome_cov.txt"
    run:
        dedup_pattern = re.compile("Total count of deduplicated leftover sequences:\s+(\d+)")
        dedup_reads = find_number(dedup_pattern, input.dedup_report)
        genome_size = config["genome_size"]
        read_length = config["read_length"]
        max_coverage = dedup_reads * (read_length * 2 - 25) / genome_size / 2
        cnt_c = config["num_C"] + config["num_G"]
        cnt_cg = config["num_CG"] * 2
        s_cnt_c = 0
        s_cnt_cg = 0
        with open(input.bismark_cx, "r") as f:
            for line in f:
                cols = line.split("\t")
                c = int(cols[1]) + int(cols[2])
                s_cnt_c += c
                if cols[0] == "CG":
                    s_cnt_cg += c
        c_cov = s_cnt_c / cnt_c
        cg_cov = s_cnt_cg / cnt_cg

        with open(str(output), "w") as f:
            print("{:s}\t{:.0f}\t{:.3f}\t{:.3f}\t{:.3f}".format(wildcards.sample, dedup_reads, max_coverage, c_cov, cg_cov), file=f)


def find_number(pattern, file):
    with open(file, "r") as f:
        for line in f:
            match = re.match(pattern, line)
            if match:
                number = int(match.group(1))
    return number


