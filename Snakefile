import re
from datetime import datetime

PHIX_REF = "/echofs1/data/xiaoyu/genomes/phiX174/phiX174.fa"
LAMBDA_DIR = "/echofs1/data/xiaoyu/genomes/lambda"
HUMAN_DIR = "/echofs1/data/xiaoyu/genomes/hg38"
TMP_DIR = "/tmp"
QC_DIR = "0_fastqc"
TRIM_DIR = "1_trim"
PHIX_DIR = "1_phix"
BISMARK_LAMBDA = "2_bismark_lambda"
BISMARK = "2_bismark"
PRESEQ = "3_preseq"

shell.prefix("module load FastQC/0.11.5 multiqc/1.7 trim_galore/0.6.6 samtools/1.9 bwa/0.7.15 bismark/0.18.1 preseq/3.1.2 bedtools/2.27.1 R/3.6.1;")

# pattern = re.compile(r'\.fastq.gz$')
SAMPLES = list(map(lambda x: x[:-12], filter(lambda y: y.endswith('.fastq.gz'), os.listdir("."))))
print(SAMPLES)
rule all:
    input:
        expand(QC_DIR + "/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand(QC_DIR + "/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand(QC_DIR + "/{sample}_R1.html", sample=SAMPLES),
        expand(QC_DIR + "/{sample}_R2.html", sample=SAMPLES),
        expand(TRIM_DIR + "/{sample}_R1.fq.gz", sample=SAMPLES),
        expand(TRIM_DIR + "/{sample}_R2.fq.gz", sample=SAMPLES),
        expand(TRIM_DIR + "/{sample}_R1_trimmed_fastqc.zip", sample=SAMPLES),
        expand(TRIM_DIR + "/{sample}_R2_trimmed_fastqc.zip", sample=SAMPLES),
        expand(TRIM_DIR + "/{sample}_R1_trimming_report.txt", sample=SAMPLES),
        expand(TRIM_DIR + "/{sample}_R2_trimming_report.txt", sample=SAMPLES),
        expand(PHIX_DIR + "/{sample}.bwa_phiX.{ext}", sample=SAMPLES, ext=["bam", "txt"]),
        expand(BISMARK_LAMBDA + "/{sample}_bismark_bt2.CXme.txt", sample=SAMPLES),
        expand(BISMARK + "/{sample}_bismark_bt2.CXme.txt", sample=SAMPLES),
        expand(BISMARK + "/{sample}_bismark_bt2_pe.bam", sample=SAMPLES)



rule fastqc:
    input:
        R1 = "{sample}_R1.fastq.gz",
        R2 = "{sample}_R2.fastq.gz"
    threads:
        16
    params:
        dir = directory(QC_DIR)
    output:
        R1 = QC_DIR + "/{sample}_R1_fastqc.zip",
        R2 = QC_DIR + "/{sample}_R2_fastqc.zip"
    log:
        QC_DIR + "/{sample}.fastqc.log"
    shell:
        "fastqc -o {params.dir} --noextract --nogroup -t {threads} {input.R1} {input.R2} &> {log}"

rule multiqc:
    input:
        R1 = QC_DIR + "/{sample}_R1_fastqc.zip",
        R2 = QC_DIR + "/{sample}_R2_fastqc.zip"
    params:
        R1 = QC_DIR + "/{sample}_R1",
        R2 = QC_DIR + "/{sample}_R2"
    output:
        R1 = QC_DIR + "/{sample}_R1.html", 
        R2 = QC_DIR + "/{sample}_R2.html"
    log:
        R1 = QC_DIR + "/{sample}_R1.multiqc_fastqc.log",
        R2 = QC_DIR + "/{sample}_R2.multiqc_fastqc.log"
    run:
        shell("multiqc -m fastqc -n {params.R1} -v {input.R1} &> {log.R1}")
        shell("multiqc -m fastqc -n {params.R2} -v {input.R2} &> {log.R2}")

rule trim:
    input:
        R1 = "{sample}_R1.fastq.gz",
        R2 = "{sample}_R2.fastq.gz"
    output:
        trim1_fq = TRIM_DIR + "/{sample}_R1.fq.gz",
        trim2_fq = TRIM_DIR + "/{sample}_R2.fq.gz",
        # outs_R1 = multiext(TRIM_DIR + "/{sample}_R1", ".fq.gz_trimming_report.txt", "_val_1_fastqc.html", "_val_1_fastqc.zip"),
        # outs_R2 = multiext(TRIM_DIR + "/{sample}_R2", ".fq.gz_trimming_report.txt", "_val_2_fastqc.html", "_val_2_fastqc.zip"),
        out = expand(TRIM_DIR + "/{sample}_R{Rn}_trimmed_fastqc.{Ext}", Rn=[1, 2], Ext=["zip", "html"], allow_missing=True),
        report = expand(TRIM_DIR + "/{sample}_R{Rn}_trimming_report.txt", Rn=[1, 2], allow_missing=True)
    log:
        TRIM_DIR + "/{sample}.trim_galore.log"
    threads:
        16
    params:
        dir = directory(TRIM_DIR),
        trim_base1 = 10,
        trim_base2 = 15,
        tmp_r1 = TRIM_DIR + "/{sample}_R1_val_1.fq.gz",
        tmp_r2 = TRIM_DIR + "/{sample}_R2_val_2.fq.gz",
        R1 = expand(TRIM_DIR + "/{sample}_R1_val_1_fastqc.{Ext}", Ext=["zip", "html"], allow_missing=True),
        R2 = expand(TRIM_DIR + "/{sample}_R2_val_2_fastqc.{Ext}", Ext=["zip", "html"], allow_missing=True),
        report = expand(TRIM_DIR + "/{sample}_R{Rn}.fq.gz_trimming_report.txt", Rn=[1, 2], allow_missing=True)

    run:
        shell("""trim_galore -q 20 --phred33 --fastqc --fastqc_args "-o {params.dir} --noextract --nogroup" \
            --illumina --stringency 1 -e 0.1 --length 20 \
            --clip_R1 {params.trim_base1} --clip_R2 {params.trim_base1} \
            -o {params.dir} \
            -j {threads} \
            --paired --retain_unpaired -r1 21 -r2 21 {input.R1} {input.R2} &>{log}""")
        shell("mv {params.tmp_r1} {output.trim1_fq}")
        shell("mv {params.tmp_r2} {output.trim2_fq}")
        shell("""rename "s/_val_[12]/_trimmed/g" {params.R1} {params.R2}""")
        shell("""rename "s/.fq.gz//g" {params.report}""")

# rule trim_rename:
#     input:
#         R1 = expand(TRIM_DIR + "/{sample}_R1_val_1_fastqc.{Ext}", Ext=["zip", "html"], allow_missing=True),
#         R2 = expand(TRIM_DIR + "/{sample}_R2_val_2_fastqc.{Ext}", Ext=["zip", "html"], allow_missing=True),
#         report = expand(TRIM_DIR + "/{sample}_R{Rn}.fq.gz_trimming_report.txt", Rn=[1, 2], allow_missing=True)
#     output:
#         out = expand(TRIM_DIR + "/{sample}_R{Rn}_trimmed_fastqc.{Ext}", Rn=[1, 2], Ext=["zip", "html"], allow_missing=True),
#         report = expand(TRIM_DIR + "/{sample}_R{Rn}_trimming_report.txt", Rn=[1, 2], allow_missing=True)
#     run:
#         shell("""rename "s/_val_[12]/_trimmed/g" {input.R1} {input.R2}""")
#         shell("""rename "s/.fq.gz//g" {input.report}""")

rule phix:
    input:
        ref = PHIX_REF,
        R1_fq = TRIM_DIR + "/{sample}_R1.fq.gz",
        R2_fq = TRIM_DIR + "/{sample}_R2.fq.gz"
    params:
        dir = directory(PHIX_DIR),
    output:
        bam = PHIX_DIR + "/{sample}.bwa_phiX.bam",
        output = PHIX_DIR + "/{sample}.bwa_phiX.txt"
    log:
        PHIX_DIR + "/{sample}.bwa_phiX.log"
    threads:
        4
    run:
        total = shell("""zcat {input.R1_fq} | grep -c "^@" """, read=True).rstrip()
        shell("bwa mem -t {threads} {input.ref} {input.R1_fq} {input.R2_fq} 2> {log} | samtools view -b -o {output.bam} -F 4 -@ {threads} -")
        phi = shell("samtools view -c -@ {threads} {output.bam}", read=True).rstrip()
        # phi_rate = shell("""awk -v c={phi} -v t={total} "BEGIN{{ printf(\\"%.7g\\", c/(t*2))}}" """, read=True).rstrip()
        # shell("""echo -e "{wildcards.sample}\t{total}\t{phi}\t{phi_rate}" > {output.output}""")
        phi_rate = phi / (total * 2)
        with open("{output.output}", "a") as f:
            print("{:s}\t{:d}\t{:d}\t{:.7g}".format(wildcards.smaple, total, phi, phi_rate), file=f)
            # print("{wildcard.sample}\t{total}\t{phi}\t{phi_rate}", file=f)

rule bismark_lambda:
    input:
        # fqs = expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True),
        R1_fq = TRIM_DIR + "/{sample}_R1.fq.gz",
        R2_fq = TRIM_DIR + "/{sample}_R2.fq.gz"
    threads:
        6
    params:
        ref_dir = directory(LAMBDA_DIR),
        min_insert = 0,
        max_insert = 2000,
        out_dir = directory(BISMARK_LAMBDA),
        report = "{sample}_bismark_bt2_PE_report.html",
        cg_pe = BISMARK_LAMBDA + "/CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        ch_pe = BISMARK_LAMBDA + "/Non_CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        merged = BISMARK_LAMBDA + "/{sample}_bismark_bt2.extracted.txt.gz",
        bedGraph = "{sample}.bedGraph.gz",
        cov = "{sample}.bismark.cov.gz",
        cx_report = "{sample}.CX_report.txt.gz"

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
        cx_me = BISMARK_LAMBDA + "/{sample}_bismark_bt2.CXme.txt"

    run:
        print("Started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        # Mapping with bismark/bowtie2
        # Note --bowtie2 and -p $nthreads are both SLOWER than single threaded bowtie1
        print("1. Mapping to reference with bismark/bowtie2... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark -q -I {min_insert} -X {max_insert} --parallel 2 -p {threads} \
            --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 \
            -o {params.out_dir} --temp_dir {TMP_DIR} --gzip --nucleotide_coverage \
            {ref_dir} -1 {R1_fq} -2 {R2_fq} &>{log.bismark_pe}")
        for f in (filter(lambda x: x.startswith("{wildcards.sample}_R1_bismark_bt2_"), os.listdir("{params.out_dir}"))):
            shell("""rename "s/_R1_bismark_bt2/_bismark_bt2/g" {params.out_dir}/{f}""")

        # Deduplicate reads
        print("-- 2. Deduplicating aligned reads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("deduplicate_bismark -p --bam {output.bam_pe}  &>${log.dedup_pe}")

        # Run methylation extractor for the sample
        print("-- 3. Analyse methylation in {output.bam_dedup_pe} using $CPU threads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report \
            -o {out_dir} --gzip --parallel {threads} \
            {output.bam_dedup_pe} &>{log.methx_pe}")

        # Generate HTML Processing Report
        print("-- 4. Generate bismark HTML processing report file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark2report -o {params.report} --dir {out_dir} \
           --alignment_report {output.alignment_report} \
           --dedup_report {output.dedup_report} \
           --splitting_report {output.splitting_report} \
           --mbias_report {output.mbias_report} \
           --nucleotide_report {output.nucleotide_report} &>{log.html}")

        # Generate bedGraph file
        print("-- 5. Generate bedGraph file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("mv {params.ch_pe} {params.merged}")
        shell("cat {params.cg_pe} >>{params.merged}")

        shell("bismark2bedGraph --dir {out_dir} --cutoff 1 --CX_context --buffer_size=75G --scaffolds \
            -o {params.bedGraph} {params.merged} &>{log.bismark2bg}")
        shell("rm {out_dir}/{params.bedGraph}")  # $merged

        # Calculate average methylation levels per each CN context
        print("-- 6. Generate cytosine methylation file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("coverage2cytosine -o {params.cx_report} --dir {out_dir} --genome_folder {ref_dir} --CX_context --gzip \
              {params.cov} &>{log.cov2c}")

        shell("rm {out_dir}/{cov}")
        shell("""zcat {out_dir}/{params.cx_report} | \
            awk "BEGIN{{ca=0;cc=0;cg=0;ct=0;mca=0;mcc=0;mcg=0;mct=0}} \
                 $7~/^CA/ {{ca+=$5; mca+=$4}} \
                 $7~/^CC/ {{cc+=$5; mcc+=$4}} \
                 $7~/^CG/ {{cg+=$5; mcg+=$4}} \
                 $7~/^CT/ {{ct+=$5; mct+=$4}} \
                 END{{printf(\\"CA\t%d\t%d\t%.3f\n\\", ca, mca, mca/(ca+mca)); \
                     printf(\\"CC\t%d\t%d\t%.3f\n\\", cc, mcc, mcc/(cc+mcc)); \
                     printf(\\"CG\t%d\t%d\t%.3f\n\\", cg, mcg, mcg/(cg+mcg)); \
                     printf(\\"CT\t%d\t%d\t%.3f\n\\", ct, mct, mct/(ct+mct));}" >{output.cx_me}""")

        # Print the files generated
        print("-- The results...")
        shell("ls -l {out_dir}/*{wildcards.sample}*")
        print("-- Finished on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

rule bismark:
    input:
        # fqs = expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True),
        R1_fq = TRIM_DIR + "/{sample}_R1.fq.gz",
        R2_fq = TRIM_DIR + "/{sample}_R2.fq.gz"
    threads:
        6
    params:
        ref_dir = directory(LAMBDA_DIR),
        min_insert = 0,
        max_insert = 2000,
        out_dir = directory(BISMARK),
        report = "{sample}_bismark_bt2_PE_report.html",
        cg_pe = BISMARK + "/CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        ch_pe = BISMARK + "/Non_CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        merged = BISMARK + "/{sample}_bismark_bt2.extracted.txt.gz",
        bedGraph = "{sample}.bedGraph.gz",
        cov = "{sample}.bismark.cov.gz",
        cx_report = "{sample}.CX_report.txt.gz"

    log:
        bismark_pe = BISMARK + "/{sample}.bismark.log",
        dedup_pe = BISMARK + "/{sample}.deduplicate_bismark.log",
        methx_pe = BISMARK + "/{sample}_pe.bismark_methylation_extractor.log",
        html = BISMARK + "/{sample}_pe.bismark2report.log",
        bismark2bg = BISMARK + "/{sample}.bismark2bedGraph.log",
        cov2c = BISMARK + "/{sample}.coverage2cytosine.log"

    output:
        read1_unmapped_fq = BISMARK + "/{sample}_R1.fq.gz_unmapped_reads_1.fq.gz",
        read2_unmapped_fq = BISMARK + "/{sample}_R2.fq.gz_unmapped_reads_2.fq.gz",
        read_se1_fq = BISMARK + "/{sample}_R1.fq.gz",
        read_se2_fq = BISMARK + "/{sample}_R2.fq.gz",
        bam_pe = BISMARK + "/{sample}_bismark_bt2_pe.bam",
        bam_dedup_pe = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.bam",
        report = BISMARK + "/{sample}_bismark_bt2_PE_report.html",
        alignment_report = BISMARK + "/{sample}_bismark_bt2_PE_report.txt",
        dedup_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplication_report.txt",
        splitting_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated_splitting_report.txt",
        mbias_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.M-bias.txt",
        nucleotide_report = BISMARK + "/{sample}_bismark_bt2_pe.nucleotide_stats.txt",
        cx_me = BISMARK + "/{sample}_bismark_bt2.CXme.txt"

    run:
        print("Started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        # Mapping with bismark/bowtie2
        # Note --bowtie2 and -p $nthreads are both SLOWER than single threaded bowtie1
        print("1. Mapping to reference with bismark/bowtie2... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark -q -I {min_insert} -X {max_insert} --parallel 2 -p {threads} \
            --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 \
            -o {params.out_dir} --temp_dir {TMP_DIR} --gzip --nucleotide_coverage \
            {ref_dir} -1 {R1_fq} -2 {R2_fq} &>{log.bismark_pe}")
        for f in (filter(lambda x: x.startswith("{wildcards.sample}_R1_bismark_bt2_"), os.listdir("{params.out_dir}"))):
            shell("""rename "s/_R1_bismark_bt2/_bismark_bt2/g" {params.out_dir}/{f}""")

        # Deduplicate reads
        print("-- 2. Deduplicating aligned reads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("deduplicate_bismark -p --bam {output.bam_pe}  &>${log.dedup_pe}")

        # Run methylation extractor for the sample
        print("-- 3. Analyse methylation in {output.bam_dedup_pe} using $CPU threads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report \
            -o {out_dir} --gzip --parallel {threads} \
            {output.bam_dedup_pe} &>{log.methx_pe}")

        # Generate HTML Processing Report
        print("-- 4. Generate bismark HTML processing report file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark2report -o {params.report} --dir {out_dir} \
           --alignment_report {output.alignment_report} \
           --dedup_report {output.dedup_report} \
           --splitting_report {output.splitting_report} \
           --mbias_report {output.mbias_report} \
           --nucleotide_report {output.nucleotide_report} &>{log.html}")

        # Generate bedGraph file
        print("-- 5. Generate bedGraph file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("mv {params.ch_pe} {params.merged}")
        shell("cat {params.cg_pe} >>{params.merged}")

        shell("bismark2bedGraph --dir {out_dir} --cutoff 1 --CX_context --buffer_size=75G --scaffolds \
            -o {params.bedGraph} {params.merged} &>{log.bismark2bg}")
        shell("rm {out_dir}/{params.bedGraph}")  # $merged

        # Calculate average methylation levels per each CN context
        print("-- 6. Generate cytosine methylation file... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("coverage2cytosine -o {params.cx_report} --dir {out_dir} --genome_folder {ref_dir} --CX_context --gzip \
              {params.cov} &>{log.cov2c}")

        shell("rm {out_dir}/{cov}")
        shell("""zcat {out_dir}/{params.cx_report} | \
            awk "BEGIN{{ca=0;cc=0;cg=0;ct=0;mca=0;mcc=0;mcg=0;mct=0}} \
                 $7~/^CA/ {{ca+=$5; mca+=$4}} \
                 $7~/^CC/ {{cc+=$5; mcc+=$4}} \
                 $7~/^CG/ {{cg+=$5; mcg+=$4}} \
                 $7~/^CT/ {{ct+=$5; mct+=$4}} \
                 END{{printf(\\"CA\t%d\t%d\t%.3f\n\\", ca, mca, mca/(ca+mca)); \
                     printf(\\"CC\t%d\t%d\t%.3f\n\\", cc, mcc, mcc/(cc+mcc)); \
                     printf(\\"CG\t%d\t%d\t%.3f\n\\", cg, mcg, mcg/(cg+mcg)); \
                     printf(\\"CT\t%d\t%d\t%.3f\n\\", ct, mct, mct/(ct+mct));}" >{output.cx_me}""")

        # Print the files generated
        print("-- The results...")
        shell("ls -l {out_dir}/*{wildcards.sample}*")
        print("-- Finished on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

rule preseq:
    input:
        BISMARK + "/{sample}_bismark_bt2_pe.bam"
    params:
        dir = directory(PRESEQ),
        bam_sorted = PRESEQ + "/{sample}.sorted.bam"
    output:
        PRESEQ + "/{sample}.preseq_lc_extrap.txt"
    log:
        PRESEQ + "/{sample}.preseq_lc_extrap.log"
    threads:
        workflow.cores * 0.5 - 1
    run:
        shell("samtools sort -m 2G -o {params.bam_sorted} -T /tmp/{wildcards.sample} -@ threads {input}")
        shell("preseq lc_extrap -o {output} -B -P -D {params.bam_sorted} 2>{log}")
        shell("rm {params.bam_sorted}")

rule insert_cpg_bias:



# rule jsfile:
#     input:
#         TEMPLATE_JS
#     output:
#         REF + "/" + REF + ".js"
#     params:
#         ref = REF, server = SERVER
#     shell:
#         """sed "s:refgenome:{params.ref}:g;s|server|{params.server}|g" {input} > {output}"""
