import re
from os.path import exists
import socket
from datetime import datetime
import fnmatch
import gzip
from math import floor
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

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

SUFFIX = '.fq.gz'

suffix_length = len(SUFFIX) + 3
SAMPLES = list(set(map(lambda x: x[:-suffix_length], filter(lambda y: y.endswith(SUFFIX), os.listdir(FASTQ)))))
eprint(SAMPLES)

hostname = socket.gethostname()
if len(re.findall("compute1", hostname)) > 0:
    server = "ris"
elif len(re.findall("^n\\d+$", hostname)) == 0:
    server = "wanglab"
else:
    server = "htcf"

if "params" in config.keys():
    include: config["params"]

if server == "htcf":
    shell.prefix("module load cutadapt/1.16-python3 trimgalore/0.6.0-python-3.6.5-java-11 samtools/1.9 bwa/0.7.15 bowtie2/2.3.4.1 preseq/2.0.2 bedtools/2.27.1 htslib/1.3.1 r/3.6.3-python-3.6.5-java-11;")
    pipe_path = "/home/fanc/software/wgbs"
elif server == "wanglab":
    shell.prefix("module load FastQC/0.11.5 multiqc/1.7 trim_galore/0.6.6 samtools/1.9 bwa/0.7.15 bismark/0.23.1 preseq/3.1.2 bedtools/2.27.1 htslib/1.3.1 R/3.6.1;")
    pipe_path = "/bar/cfan/software/wgbs"
else:
    pipe_path = "/home/fanc/software/wgbs"

trash = os.system("rm -rf scripts")
trash = os.system("ln -s " + pipe_path + "/scripts ./")

rule all:
    input:
        PRESEQ + "/aggregate_preseq.txt",
        INSERT_CPG_BIAS + "/insert_length_all.txt",
        TRACKS + "/distribution.sum.txt",
        QC_DIR + "/fq_aggregated/summary.txt",
        QC_DIR + "/fq_aggregated/statistics.txt",
        TRIM_DIR + "/fq_aggregated/summary.txt",
        TRIM_DIR + "/fq_aggregated/statistics.txt",
        TRIM_DIR + "/multiqc_report.html",
        QC_DIR + "/multiqc_report.html",
        "all.genome_cov.txt",
        "alignment.txt",
        "conversion.txt",
        "reads_number.pdf",
        "reads_perc.pdf",
        "conversion_rate.pdf",
        "library_complexity.pdf",
        "insert_length_all.pdf",
        "CpGbias_all.txt",
        "distribution.pdf",

rule shallow:
    input:
        QC_DIR + "/fq_aggregated/summary.txt",
        QC_DIR + "/fq_aggregated/statistics.txt",
        TRIM_DIR + "/fq_aggregated/summary.txt",
        TRIM_DIR + "/fq_aggregated/statistics.txt",
        TRIM_DIR + "/multiqc_report.html",
        QC_DIR + "/multiqc_report.html",
        "alignment.txt",
        "conversion.txt",
        "reads_number.pdf",
        "reads_perc.pdf",
        "conversion_rate.pdf",
        "insert_length_all.pdf",
        "CpGbias_all.txt",

rule delete_aggr:
    output:
        "del_aggr.out"
    params:
        PRESEQ + "/aggregate_preseq.txt",
        INSERT_CPG_BIAS + "/insert_length_all.txt",
        TRACKS + "/distribution.sum.txt",
        QC_DIR + "/fq_aggregated/summary.txt",
        QC_DIR + "/fq_aggregated/statistics.txt",
        TRIM_DIR + "/fq_aggregated/summary.txt",
        TRIM_DIR + "/fq_aggregated/statistics.txt",
        "all.genome_cov.txt",
        "alignment.txt",
        "conversion.txt",
        "reads_number.pdf",
        "reads_perc.pdf",
        "conversion_rate.pdf",
        "library_complexity.pdf",
        "insert_length_all.pdf",
        "CpGbias_all.txt",
        "distribution.pdf"
    shell:
        "rm -rf {params} && "
        "echo 'aggr deleted' > {output}"       

rule multiqc:
    input:
        fastqc = QC_DIR,
        trim = TRIM_DIR
    output:
        QC_DIR + "/multiqc_report.html",
        TRIM_DIR + "/multiqc_report.html"
    shell:
        "cd {input.fastqc} && multiqc -f . && "
        "cd ../{input.trim} && multiqc -f . &&  cd ../"
rule extract_fastqc:
    input:
        expand(QC_DIR + "/{sample}_R{Rn}_fastqc.zip", Rn=[1, 2], sample=SAMPLES)
    params:
        QC_DIR
    output:
        summary = QC_DIR + "/fq_aggregated/summary.txt",
        stats = QC_DIR + "/fq_aggregated/statistics.txt"
    shell:
        "scripts/extract_fastqc.sh {params}"

rule extract_trim_fastqc:
    input:
        expand(TRIM_DIR + "/{sample}_R{Rn}_val_{Rn}_fastqc.zip", Rn=[1, 2], sample=SAMPLES)
    params:
        TRIM_DIR
    output:
        summary = TRIM_DIR + "/fq_aggregated/summary.txt",
        stats = TRIM_DIR + "/fq_aggregated/statistics.txt"
    shell:
        "scripts/extract_fastqc.sh {params}"

rule aggregation:
    input:
        trim_reports = expand(TRIM_DIR + "/{sample}_R2_trimming_report.txt", sample=SAMPLES),
        bismark_reports = expand(BISMARK + "/{sample}_bismark_bt2_PE_report.txt", sample=SAMPLES),
        dedup_reports = expand(BISMARK + "/{sample}_bismark_bt2_pe.deduplication_report.txt", sample=SAMPLES),
        phix_txts = expand(PHIX_DIR + "/{sample}.bwa_phiX.txt", sample=SAMPLES),
        lambda_cxs = expand(BISMARK_LAMBDA + "/{sample}_bismark_bt2.CXme.txt", sample=SAMPLES),
        # bismark_cxs = expand(BISMARK + "/{sample}_bismark_bt2.CXme.txt", sample=SAMPLES),

    output:
        alignment_sum = "alignment.txt",
        conversion_sum = "conversion.txt",
    run:
        total_pattern = re.compile("Total number of sequences analysed for the sequence pair length validation:\s+(\d+)")
        trimmed_pattern = re.compile("Number of sequence pairs removed because at least one read was shorter than the length cutoff.+:\s+(\d+)")
        after_trim_pattern = re.compile("Sequence pairs analysed in total:\s+(\d+)$")
        uniq_pattern = re.compile("Number of paired-end alignments with a unique best hit:\s+(\d+)$")
        multi_pattern = re.compile("Sequence pairs did not map uniquely:\s+(\d+)$")
        bam_pattern = re.compile("Total number of alignments analysed in .+:\s+(\d+)$")
        dup_pattern = re.compile("Total number duplicated alignments removed:\s+(\d+)")
        dedup_pattern = re.compile("Total count of deduplicated leftover sequences:\s+(\d+)")
        second_pattern = re.compile("^\S+\s+(\d+)\s+\d+\s+")
        third_pattern = re.compile("^\S+\s+\d+\s+(\d+)\s+")

        ca_pattern = re.compile("^CA\s+(\d+)\s+\d+")
        mca_pattern = re.compile("^CA\s+\d+\s+(\d+)")
        cc_pattern = re.compile("^CC\s+(\d+)\s+\d+")
        mcc_pattern = re.compile("^CC\s+\d+\s+(\d+)")
        cg_pattern = re.compile("^CG\s+(\d+)\s+\d+")
        mcg_pattern = re.compile("^CG\s+\d+\s+(\d+)")
        ct_pattern = re.compile("^CT\s+(\d+)\s+\d+")
        mct_pattern = re.compile("^CT\s+\d+\s+(\d+)")

        with open(output.alignment_sum, "w") as a:
            print("sample\ttotal\ttrimmed\tafter_trim\tmapped\tmapped_perc\tmulti\tuniq\tunique_perc\tbam_reads\tdup_reads\tdedup\tdedup_perc", file=a)
            for sample, trim_report, bismark_report, dedup_report in zip(SAMPLES, input.trim_reports, input.bismark_reports, input.dedup_reports):
                total_reads = find_number(total_pattern, trim_report)
                trimmed_reads = find_number(trimmed_pattern, trim_report)
                after_trim_reads = find_number(after_trim_pattern, bismark_report)
                uniq_reads = find_number(uniq_pattern, bismark_report)
                multi_reads = find_number(multi_pattern, bismark_report)
                mapped_reads = uniq_reads + multi_reads
                bam_reads = find_number(bam_pattern, dedup_report)
                dup_reads = find_number(dup_pattern, dedup_report)
                dedup_reads = find_number(dedup_pattern, dedup_report)
                mapped_perc = mapped_reads / total_reads
                uniq_perc = uniq_reads / mapped_reads
                dedup_perc = dedup_reads / bam_reads
                print("{:s}\t{:d}\t{:d}\t{:d}\t{:d}\t{:.2%}\t{:d}\t{:d}\t{:.2%}\t{:d}\t{:d}\t{:d}\t{:.2%}".format(
                    sample, total_reads, trimmed_reads, after_trim_reads, mapped_reads, mapped_perc, multi_reads, uniq_reads, uniq_perc, bam_reads, dup_reads, dedup_reads, dedup_perc), file=a)
        bismark_cxs = expand(BISMARK + "/{sample}_bismark_bt2.CXme.txt", sample=SAMPLES)
        if exists(bismark_cxs[0]):
            with open(output.conversion_sum, "w") as c:
                print("sample\tafter_trim_phix\tphix_reads\tphix_perc\tlambda_total\tlambda_conv\tlambda_rate\tCA_met\tCC_met\tCT_met\tCH_conv_rate", file=c)
                for sample, phix, lambda_cx, bismark_cx in zip(SAMPLES, input.phix_txts, input.lambda_cxs, bismark_cxs):
                    after_trim_phix = find_number(second_pattern, phix)
                    phix_reads = find_number(third_pattern, phix)
                    phix_perc = phix_reads / after_trim_phix / 2
                    lambda_ca_reads = find_number(ca_pattern, lambda_cx)
                    lambda_mca_reads = find_number(mca_pattern, lambda_cx)
                    lambda_cc_reads = find_number(cc_pattern, lambda_cx)
                    lambda_mcc_reads = find_number(mcc_pattern, lambda_cx)
                    lambda_cg_reads = find_number(cg_pattern, lambda_cx)
                    lambda_mcg_reads = find_number(mcg_pattern, lambda_cx)
                    lambda_ct_reads = find_number(ct_pattern, lambda_cx)
                    lambda_mct_reads = find_number(mct_pattern, lambda_cx)
                    lambda_total = lambda_ca_reads + lambda_mca_reads + lambda_cc_reads + lambda_mcc_reads + lambda_cg_reads + lambda_mcg_reads + lambda_ct_reads + lambda_mct_reads
                    lambda_conv = lambda_ca_reads + lambda_cc_reads + lambda_cg_reads + lambda_ct_reads
                    lambda_rate = lambda_conv / lambda_total

                    ca_reads = find_number(ca_pattern, bismark_cx)
                    mca_reads = find_number(mca_pattern, bismark_cx)
                    cc_reads = find_number(cc_pattern, bismark_cx)
                    mcc_reads = find_number(mcc_pattern, bismark_cx)
                    cg_reads = find_number(cg_pattern, bismark_cx)
                    mcg_reads = find_number(mcg_pattern, bismark_cx)
                    ct_reads = find_number(ct_pattern, bismark_cx)
                    mct_reads = find_number(mct_pattern, bismark_cx)

                    ca_met = mca_reads / (mca_reads + ca_reads)
                    cc_met = mcc_reads / (mcc_reads + cc_reads)
                    cg_met = mcg_reads / (mcg_reads + cg_reads)
                    ct_met = mct_reads / (mct_reads + ct_reads)
                    ch_conv_rate = (ca_reads + cc_reads + ct_reads) / (ca_reads + mca_reads + cc_reads + mcc_reads + ct_reads + mct_reads)
                    print("{:s}\t{:d}\t{:d}\t{:.2e}\t{:d}\t{:d}\t{:.2%}\t{:.2%}\t{:.2%}\t{:.2%}\t{:.2%}".format(
                        sample, after_trim_phix, phix_reads, phix_perc, lambda_total, lambda_conv, lambda_rate, ca_met, cc_met, ct_met, ch_conv_rate), file=c)
        else:
            with open(output.conversion_sum, "w") as c:
                print("sample\tafter_trim_phix\tphix_reads\tphix_perc\tlambda_total\tlambda_conv\tlambda_rate", file=c)
                for sample, phix, lambda_cx in zip(SAMPLES, input.phix_txts, input.lambda_cxs):
                    after_trim_phix = find_number(second_pattern, phix)
                    phix_reads = find_number(third_pattern, phix)
                    phix_perc = phix_reads / after_trim_phix / 2
                    lambda_ca_reads = find_number(ca_pattern, lambda_cx)
                    lambda_mca_reads = find_number(mca_pattern, lambda_cx)
                    lambda_cc_reads = find_number(cc_pattern, lambda_cx)
                    lambda_mcc_reads = find_number(mcc_pattern, lambda_cx)
                    lambda_cg_reads = find_number(cg_pattern, lambda_cx)
                    lambda_mcg_reads = find_number(mcg_pattern, lambda_cx)
                    lambda_ct_reads = find_number(ct_pattern, lambda_cx)
                    lambda_mct_reads = find_number(mct_pattern, lambda_cx)
                    lambda_total = lambda_ca_reads + lambda_mca_reads + lambda_cc_reads + lambda_mcc_reads + lambda_cg_reads + lambda_mcg_reads + lambda_ct_reads + lambda_mct_reads
                    lambda_conv = lambda_ca_reads + lambda_cc_reads + lambda_cg_reads + lambda_ct_reads
                    lambda_rate = lambda_conv / lambda_total

                    print("{:s}\t{:d}\t{:d}\t{:.2e}\t{:d}\t{:d}\t{:.2%}".format(
                        sample, after_trim_phix, phix_reads, phix_perc, lambda_total, lambda_conv, lambda_rate), file=c)         


rule aggregation_plots:
    input:
        "alignment.txt",
        "conversion.txt"
    output:
        "reads_number.pdf",
        "reads_perc.pdf",
        "conversion_rate.pdf"
    script:
        "scripts/summary.plots.R"

rule aggregate_preseq:
    input:
        expand(PRESEQ + "/{sample}.preseq_lc_extrap.txt", sample=SAMPLES)
    output:
        PRESEQ + "/aggregate_preseq.txt"
    run:
        with open(str(output), "w") as o:
            print("SAMPLE\tTOTAL_READS\tEXPECTED_DISTINCT\tLOWER_0.95CI\tUPPER_0.95CI", file=o)
            for sample, preseq in zip(SAMPLES, input):
                with open(preseq, "r") as f:
                    for i, line in enumerate(f):
                        if i == 0:
                            continue
                        print("{:s}\t{:s}".format(sample, line), file=o, end='')
        shell("""mv {output} {output}.tempt && awk -F '\t' 'BEGIN{{OFS = FS}} {{print $1, $2, $3}} ' {output}.tempt > {output} """)

rule library_complexity_plot:
    input:
        PRESEQ + "/aggregate_preseq.txt"
    output:
        "library_complexity.pdf"
    log:
        "log.library_complexity.txt"
    script:
        "scripts/preseq.R"

rule aggregate_inser_len:
    input:
        expand(INSERT_CPG_BIAS + "/{sample}.insert_length.txt", sample=SAMPLES)
    output:
        INSERT_CPG_BIAS + "/insert_length_all.txt"
    run:
        with open(str(output), "w") as o:
            print("sample\tposition\tdensity", file=o)
            for sample, infile in zip(SAMPLES, input):
                with open(infile, "r") as f:
                    for i, line in enumerate(f):
                        if i == 0:
                            continue
                        pos = (i - 1) * 4
                        print("{:s}\t{:d}\t{:s}".format(sample, pos, line), file=o, end='')

rule insertion_length_plot:
    input:
        INSERT_CPG_BIAS + "/insert_length_all.txt"
    output:
        "insert_length_all.pdf"
    script:
        "scripts/insert_length.R"

rule aggregate_gc_bias:
    input:
        expand(INSERT_CPG_BIAS + "/{sample}.CpGbias_1kb.R.log", sample=SAMPLES)
    output:
        "CpGbias_all.txt"
    run:
        with open(str(output), "w") as o:
            print("sample\tSpearman_correlation\tPearson_correlation", file=o)
            for sample, infile in zip(SAMPLES, input):
                with open(infile, "r") as f:
                    context = f.read()
                    print(context)
                    cor = re.search(r'^\s+cor\s*\n(.+)\s*\n', context, re.MULTILINE).group(1)
                    rho = re.search(r'^\s+rho\s*\n(.+)\s*\n', context, re.MULTILINE).group(1)
                    print(cor, rho)
                    print("{:s}\t{:s}\t{:s}".format(sample, rho, cor), file=o)

rule cover_all:
    input:
        expand(COVERAGE + "/{sample}.genome_cov.txt", sample=SAMPLES)
    output:
        "all.genome_cov.txt"
    shell:
        "echo 'sample\tdedup_reads\tmax_coverage\tc_cov\tcg_cov' > {output} && "
        "cat {input} >> {output}"

rule distribution_table:
    input:
        expand(TRACKS + "/{sample}.CG.methylC.gz", sample=SAMPLES)
    output:
        TRACKS + "/distribution.sum.txt"
    run:
        with open(str(output), "w") as o:
            print("sample\tperc\tcount", file=o)
            for sample, infile in zip(SAMPLES, input):
                with gzip.open(infile, 'rt') as f:
                    summary_dict = dict()
                    for i, line in enumerate(f):
                        # print("{:s}\t{:s}".format(sample, line), file=o, end='')
                        split_line = line.rstrip().split()
                        perc = split_line[4]
                        if perc in summary_dict:
                            summary_dict[perc] += 1
                        else:
                            summary_dict[perc] = 1
                    for key in summary_dict:
                        print("{:s}\t{:s}\t{:.3f}".format(sample, key, summary_dict[key] / (i + 1)), file=o)

rule distribution_plot:
    input:
        TRACKS + "/distribution.sum.txt"
    output:
        "distribution.pdf"
    script:
        "scripts/distribution.R"


def find_number(pattern, file):
    with open(file, "r") as f:
        for line in f:
            match = re.match(pattern, line)
            if match:
                number = int(match.group(1))
    return number


