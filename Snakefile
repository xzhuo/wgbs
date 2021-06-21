import re
from datetime import datetime

# PHIX_REF = "/echofs1/data/xiaoyu/genomes/phiX174/phiX174.fa"
# LAMBDA_DIR = "/echofs1/data/xiaoyu/genomes/lambda"
# HUMAN_DIR = "/echofs1/data/xiaoyu/genomes/hg38"
# HG_QC = "/echofs1/data/xiaoyu/genomes/hg38/wgbs_qc"
PHIX_REF = "/scratch/genomes/phiX174/bwa_index/phiX174.fa"
LAMBDA_DIR = "/scratch/genomes/lambda"
HUMAN_DIR = "/scratch/genomes/hg38/bismark"
HG_QC = "/scratch/genomes/hg38/wgbs_qc"
TMP_DIR = "/tmp"
QC_DIR = "0_fastqc"
TRIM_DIR = "1_trim"
PHIX_DIR = "1_phix"
BISMARK_LAMBDA = "2_bismark_lambda"
BISMARK = "2_bismark"
PRESEQ = "3_preseq"
INSERT_CPG_BIAS = "4_insert_cpg_bias"
COVERAGE = "5_coverage"
TRACKS = "6_tracks"

shell.prefix("module load FastQC/0.11.5 multiqc/1.7 trim_galore/0.6.6 samtools/1.9 bwa/0.7.15 bismark/0.18.1 preseq/3.1.2 bedtools/2.27.1 htslib/1.3.1 R/3.6.1;")

# pattern = re.compile(r'\.fastq.gz$')
SUFFIX = '.fq.gz'
suffix_length = len(SUFFIX) + 3
SAMPLES = set(map(lambda x: x[:-suffix_length], filter(lambda y: y.endswith(SUFFIX), os.listdir("."))))
# print(SAMPLES)
rule all:
    input:
        expand(QC_DIR + "/{sample}_R{n}{ext}", sample=SAMPLES, n=[1, 2], ext=["_fastqc.zip", ".html"]),
        # expand(QC_DIR + "/{sample}_R2_fastqc.zip", sample=SAMPLES),
        # expand(QC_DIR + "/{sample}_R{n}.html", sample=SAMPLES, n=[1, 2]),
        # expand(QC_DIR + "/{sample}_R2.html", sample=SAMPLES),
        expand(TRIM_DIR + "/{sample}_R{n}{ext}", sample=SAMPLES, n=[1, 2], ext=[".fq.gz", "_trimmed_fastqc.zip", "_trimmed_fastqc.html", "_trimming_report.txt"]),
        # expand(TRIM_DIR + "/{sample}_R2.fq.gz", sample=SAMPLES),
        # expand(TRIM_DIR + "/{sample}_R{n}_trimmed_fastqc.zip", sample=SAMPLES, n=[1, 2]),
        # # expand(TRIM_DIR + "/{sample}_R2_trimmed_fastqc.zip", sample=SAMPLES),
        # expand(TRIM_DIR + "/{sample}_R{n}_trimming_report.txt", sample=SAMPLES, n=[1, 2]),
        # expand(TRIM_DIR + "/{sample}_R2_trimming_report.txt", sample=SAMPLES),
        expand(PHIX_DIR + "/{sample}.bwa_phiX.{ext}", sample=SAMPLES, ext=["bam", "txt"]),
        expand(BISMARK_LAMBDA + "/{sample}_bismark_bt2.CXme.txt", sample=SAMPLES),
        expand(BISMARK + "/{sample}_bismark_bt2{ext}", sample=SAMPLES, ext=[".CXme.txt", "_pe.bam"]),
        expand(PRESEQ + "/{sample}.preseq_lc_extrap.txt", sample=SAMPLES),
        expand(INSERT_CPG_BIAS + "/{sample}.insert_length.txt", sample=SAMPLES),
        expand(COVERAGE + "/{sample}.genome_cov.txt", sample=SAMPLES),
        expand(TRACKS + "/{sample}.{ext}", sample=SAMPLES, ext=["cov.bg.gz", "CG.methylC.gz"])


rule fastqc:
    input:
        expand("{sample}_R{n}" + SUFFIX, n=[1, 2], allow_missing=True)
    threads:
        16
    params:
        dir = directory(QC_DIR)
    output:
        expand(QC_DIR + "/{sample}_R{n}_fastqc.zip", n=[1, 2], allow_missing=True)
    log:
        QC_DIR + "/{sample}.fastqc.log"
    shell:
        "fastqc -o {params.dir} --noextract --nogroup -t {threads} {input[0]} {input[1]} &> {log}"

rule multiqc:
    input:
        QC_DIR + "/{sample}_R{n}_fastqc.zip"
    params:
        reads = QC_DIR + "/{sample}_R{n}"
    output:
        QC_DIR + "/{sample}_R{n}.html"
    log:
        QC_DIR + "/{sample}_R{n}.multiqc_fastqc.log"
    run:
        shell("multiqc -m fastqc -n {params} -v {input} &> {log}")

rule trim:
    input:
        expand("{sample}_R{n}" + SUFFIX, n=[1, 2], allow_missing=True)
        # R2 = "{sample}_R2" + SUFFIX
    output:
        trim_fq = expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True),
        # trim2_fq = TRIM_DIR + "/{sample}_R2.fq.gz",
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
        tmp = expand(TRIM_DIR + "/{sample}_R{rn}_val_{vn}.fq.gz", zip, rn=[1, 2], vn=[1, 2], allow_missing=True),
        # tmp_r2 = TRIM_DIR + "/{sample}_R2_val_2.fq.gz",
        out = TRIM_DIR + "/{sample}_R[12]_val_[12]_fastqc.*",
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
        shell("""rename "s/_val_[12]/_trimmed/g" {params.out} \
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
        4
    run:
        total = int(shell("""zcat {input.reads[0]} | grep -c "^@" """, read=True).rstrip())
        shell("bwa mem -t {threads} {input.ref} {input.reads[0]} {input.reads[0]} 2> {log} | samtools view -b -o {output.bam} -F 4 -@ {threads} -")
        phi = int(shell("samtools view -c -@ {threads} {output.bam}", read=True).rstrip())
        # phi_rate = shell("""awk -v c={phi} -v t={total} "BEGIN{{ printf(\\"%.7g\\", c/(t*2))}}" """, read=True).rstrip()
        # shell("""echo -e "{wildcards.sample}\t{total}\t{phi}\t{phi_rate}" > {output.txt}""")
        phi_rate = phi / (total * 2)
        with open(output.txt, "w") as f:
            print("{:s}\t{:d}\t{:d}\t{:.7g}".format(wildcards.sample, total, phi, phi_rate), file=f)

rule bismark_lambda:
    input:
        expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True),
        # R1_fq = TRIM_DIR + "/{sample}_R1.fq.gz",
        # R2_fq = TRIM_DIR + "/{sample}_R2.fq.gz"
    threads:
        6
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
        cx_me = BISMARK_LAMBDA + "/{sample}_bismark_bt2.CXme.txt",
        cx_report = BISMARK_LAMBDA + "/{sample}.CX_report.txt.gz"

    run:
        print("Started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        # Mapping with bismark/bowtie2
        # Note --bowtie2 and -p $nthreads are both SLOWER than single threaded bowtie1
        print("1. Mapping to reference with bismark/bowtie2... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark -q -I {params.min_insert} -X {params.max_insert} --parallel 2 -p {threads} \
            --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 \
            -o {params.out_dir} --temp_dir {TMP_DIR} --gzip --nucleotide_coverage \
            {params.ref_dir} -1 {input[0]} -2 {input[1]} &>{log.bismark_pe}")
        shell("""rename "s/_R1_bismark_bt2/_bismark_bt2/g" {params.temp}""")

        # Deduplicate reads
        print("-- 2. Deduplicating aligned reads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("deduplicate_bismark -p -o {output.bam_dedup_pe} --bam {output.bam_pe} &>{log.dedup_pe}")

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
        shell("coverage2cytosine -o {params.cx_report} --dir {params.out_dir} --genome_folder {params.ref_dir} --CX_context --gzip \
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

rule bismark:
    input:
        expand(TRIM_DIR + "/{sample}_R{n}.fq.gz", n=[1, 2], allow_missing=True),
        # R1_fq = TRIM_DIR + "/{sample}_R1.fq.gz",
        # R2_fq = TRIM_DIR + "/{sample}_R2.fq.gz"
    threads:
        6
    params:
        ref_dir = directory(HUMAN_DIR),
        min_insert = 0,
        max_insert = 2000,
        out_dir = directory(BISMARK),
        report = "{sample}_bismark_bt2_PE_report.html",
        temp = BISMARK + "/{sample}_R1_bismark_bt2_*",
        cg_pe = BISMARK + "/CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        ch_pe = BISMARK + "/Non_CpG_context_{sample}_bismark_bt2_pe.deduplicated.txt.gz",
        merged = BISMARK + "/{sample}_bismark_bt2.extracted.txt.gz",
        bedGraph = "{sample}.bedGraph.gz",
        cov = BISMARK + "/{sample}.bismark.cov.gz",
        cx_report = "{sample}.CX_report.txt.gz"

    log:
        bismark_pe = BISMARK + "/{sample}.bismark.log",
        dedup_pe = BISMARK + "/{sample}.deduplicate_bismark.log",
        methx_pe = BISMARK + "/{sample}_pe.bismark_methylation_extractor.log",
        html = BISMARK + "/{sample}_pe.bismark2report.log",
        bismark2bg = BISMARK + "/{sample}.bismark2bedGraph.log",
        cov2c = BISMARK + "/{sample}.coverage2cytosine.log"

    output:
        # read1_unmapped_fq = BISMARK + "/{sample}_R1.fq.gz_unmapped_reads_1.fq.gz",
        # read2_unmapped_fq = BISMARK + "/{sample}_R2.fq.gz_unmapped_reads_2.fq.gz",
        # read_se1_fq = BISMARK + "/{sample}_R1.fq.gz",
        # read_se2_fq = BISMARK + "/{sample}_R2.fq.gz",
        bam_pe = BISMARK + "/{sample}_bismark_bt2_pe.bam",
        bam_dedup_pe = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.bam",
        report = BISMARK + "/{sample}_bismark_bt2_PE_report.html",
        alignment_report = BISMARK + "/{sample}_bismark_bt2_PE_report.txt",
        dedup_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplication_report.txt",
        splitting_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated_splitting_report.txt",
        mbias_report = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.M-bias.txt",
        nucleotide_report = BISMARK + "/{sample}_bismark_bt2_pe.nucleotide_stats.txt",
        cx_me = BISMARK + "/{sample}_bismark_bt2.CXme.txt",
        cx_report = BISMARK + "/{sample}.CX_report.txt.gz"

    run:
        print("Started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        # Mapping with bismark/bowtie2
        # Note --bowtie2 and -p $nthreads are both SLOWER than single threaded bowtie1
        print("1. Mapping to reference with bismark/bowtie2... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("bismark -q -I {params.min_insert} -X {params.max_insert} --parallel 2 -p {threads} \
            --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 \
            -o {params.out_dir} --temp_dir {TMP_DIR} --gzip --nucleotide_coverage \
            {params.ref_dir} -1 {input[0]} -2 {input[1]} &>{log.bismark_pe}")
        shell("""rename "s/_R1_bismark_bt2/_bismark_bt2/g" {params.temp}""")

        # Deduplicate reads
        print("-- 2. Deduplicating aligned reads... started on " + datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        shell("deduplicate_bismark -p -o {output.bam_dedup_pe} --bam {output.bam_pe} &>{log.dedup_pe}")

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
        shell("coverage2cytosine -o {params.cx_report} --dir {params.out_dir} --genome_folder {params.ref_dir} --CX_context --gzip \
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
    input:
        bam_dedup = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.bam",
        bg_chr1_1kb = HG_QC + "/CpGs.hg38_chr1_1kb_win.bg.gz",
        rscript_insert = HG_QC + "/density_insert_length.R",
        rscript_cpgbias = HG_QC + "/CpGbias_1kb.R"
    params:
        dir = directory(INSERT_CPG_BIAS),
        bam_tmp = INSERT_CPG_BIAS + "/{sample}.tmp.bam",
        insert_tmp = INSERT_CPG_BIAS + "/{sample}.tmp.insert.txt.gz",
        cov_chr1 = INSERT_CPG_BIAS + "/{sample}.tmp.CpG.cov_chr1_1kb_win.txt.gz"
    threads:
        5
    output:
        INSERT_CPG_BIAS + "/{sample}.insert_length.txt"
    log:
        insert_len = INSERT_CPG_BIAS + "/{sample}.density_insert_length.R.log",
        cpgbias = INSERT_CPG_BIAS + "/{sample}.CpGbias_1kb.R.log"
    run:
        # select the first 100K alignments
        shell("samtools view -h -@ {threads} {input.bam_dedup} | \
            head -100000197 | samtools view -b -o {params.bam_tmp} -@ {threads}")
        # make temporary insert length txt file
        shell("""bamToBed -bedpe -i {params.bam_tmp} | awk -v OFS="\t" "{{print \$1,\$2,\$6,\$6-\$2}}" | gzip -nc > {params.insert_tmp}""")

        # select only chr1 
        shell("""bamToBed -bedpe -i {params.bam_tmp} | awk "\$1==\\"chr1\\"" | \
            coverageBed -counts -a {input.bg_chr1_1kb} -b stdin | awk "\$4>0 && \$5>0" | gzip -nc > {params.cov_chr1}""")

        # run rscripts
        shell("Rscript {input.rscript_insert} {params.insert_tmp} {output} &> {log.insert_len}")
        shell("Rscript {input.rscript_cpgbias} {params.cov_chr1} {wildcards.sample} &> {log.cpgbias}")

        # remove temporary files
        shell("rm {params.bam_tmp} {params.insert_tmp} {params.cov_chr1}")

rule qc_cal_genome_cov:
    input:
        BISMARK + "/{sample}_bismark_bt2.CXme.txt"
    params:
        directory(COVERAGE)
    output:
        COVERAGE + "/{sample}.genome_cov.txt"
    run:
        cnt_c = 598683433 + 600854940 - 171823 * 2
        cnt_cg = 29303965 * 2
        s_cnt_c = 0
        s_cnt_cg = 0
        with open(str(input), "r") as f:
            for line in f:
                cols = line.split("\t")
                c = int(cols[1]) + int(cols[2])
                s_cnt_c += c
                if cols[0] == "CG":
                    s_cnt_cg += c
        c_cov = s_cnt_c / cnt_c
        cg_cov = s_cnt_cg / cnt_cg
        with open(str(output), "w") as f:
            print("{:s}\t{:.3f}\t{:.3f}".format(wildcards.sample, c_cov, cg_cov), file=f)

rule track_coverage:
    input:
        bam_in = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.bam",
    params:
        dir = directory(TRACKS)
    threads:
        8
    output:
        bam_sorted = BISMARK + "/{sample}_bismark_bt2_pe.deduplicated.sorted.bam",
        cov_out = TRACKS + "/{sample}.cov.bg.gz"
    run:
        shell("samtools sort -m 2G -o {output.bam_sorted} -T /tmp/{wildcards.sample} -@ {threads} {input.bam_in}")
        shell("genomeCoverageBed -bg -ibam {output.bam_sorted} | \
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
            sort -k1,1 -k2,2n | groupBy -g 1,2,3 -c 4,5 -o sum,sum | \
            awk -F"\\t" "BEGIN{{OFS=FS}} {{mcg=sprintf(\\"%.3f\\", \$4/(\$4+\$5)); print \$1,\$2,\$3,\\"CG\\",mcg,\\"+\\",\$4+\$5 }}" | \
            bgzip > {output} && tabix -p bed {output}"""


# rule run_distCGme_cov10:

# dir_qc=/scratch/genomes/hg38/wgbs_qc
# rscript=${dir_qc}/distCGme_cov10.R


# # INPUT
# list=${workdir}/samples.txt
# sample=$( cat $list | sed "${ID}q;d" )


# dir_in=${workdir}/dmr/1_dss_input
# dss=${dir_in}/${sample}.dss.txt.gz


# # OUTPUT
# dir_out=${workdir}/7_distCGme
# out=${dir_out}/${sample}_CGme_density_cov10.txt
# log=${dir_out}/distCGme_cov10.R.${sample}.log

# # COMMANDS
# mkdir -p $dir_out
# $rscript $dss $out &>$log

