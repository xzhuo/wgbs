FROM ubuntu
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y wget r-base r-cran-ggplot2 \
    && apt-get clean

RUN Rscript -e 'install.packages(c("patchwork", "tidyr"))'

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/src/miniconda3 \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
ENV PATH=${PATH}:/usr/src/miniconda3/bin

RUN conda install -n base -c conda-forge mamba
RUN mamba create -q -y -c conda-forge -c bioconda -n snakemake snakemake=6.12.3 
# RUN bash -c "source activate snakemake && conda install -y -c bioconda fastqc samtools preseq bedtools htslib && conda clean --all -y"
## source activate only works in bash, not the default sh. 
WORKDIR /opt/apps

RUN bash -c "wget -c https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2 && bzip2 -d bwa-0.7.15.tar.bz2 && \
    tar -xf bwa-0.7.15.tar && cd bwa-0.7.15 && make && rm -rf ../bwa*tar"
RUN wget -c https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip \
    && unzip bowtie2-2.3.4.1-linux-x86_64.zip && rm -rf bowtie2-2.3.4.1-linux-x86_64.zip
RUN wget -c https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.6.tar.gz \
    && tar -xzf 0.6.6.tar.gz && rm -rf 0.6.6.tar.gz
RUN wget -c https://github.com/s-andrews/FastQC/archive/refs/tags/v0.11.9.tar.gz \
    && tar -xzf v0.11.9.tar.gz && rm -rf v0.11.9.tar.gz \
    && cd FastQC-0.11.9

RUN wget -c https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip && rm -rf fastqc_v0.11.9.zip \
    && chmod 755 FastQC/fastqc

RUN apt-get install -y default-jre && apt-get clean
RUN wget -c https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 \
    && bzip2 -d samtools-1.15.tar.bz2 && tar -xf samtools-1.15.tar && rm -rf samtools-1.15.tar \
    && cd samtools* && ./configure && make && make install && cd ../

RUN wget -c https://github.com/smithlabcode/preseq/releases/download/v3.2.0/preseq-3.2.0.tar.gz \
    && tar -xf preseq-3.2.0.tar.gz && cd preseq_v3.2.0 && make all # buildkit

# RUN wget -c http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 \
#     && bzip2 -d preseq_linux_v2.0.tar.bz2 && tar -xf preseq_linux_v2.0.tar && rm -rf preseq_linux_v2.0.tar
RUN wget -c https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
    && mkdir bin && mv bedtools.static.binary bin/bedtools && chmod 755 bin/bedtools
RUN wget -c https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2 \
    && bzip2 -d htslib-1.15.tar.bz2 && rm -rf htslib-1.15.tar.bz2 && tar -xf htslib-1.15.tar \
    && cd htslib* && ./configure && make && make install && cd ../

RUN wget -c https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz \
    && tar -xzf 0.6.6.tar.gz

RUN pip install multiqc
RUN pip install cutadapt
RUN rm -rf *.tar *.gz FastQC-0.11.9
ENV PATH="${PATH}:/opt/apps/FastQC:/opt/apps/TrimGalore-0.6.6:/opt/apps/bin:/opt/apps/bowtie2-2.3.4.1-linux-x86_64:/opt/apps/bwa-0.7.15:/opt/apps/htslib-1.15:/opt/apps/preseq_v3.2.0:/opt/apps/samtools-1.15:/opt/apps/TrimGalore-0.6.6"
RUN apt-get install -y sudo rename && apt-get clean

RUN wget -c https://github.com/FelixKrueger/Bismark/archive/refs/tags/0.23.1.tar.gz \
    && tar -xzf 0.23.1.tar.gz && rm -rf 0.23.1.tar.gz
ENV PATH="${PATH}:/opt/apps/Bismark-0.23.1"

RUN apt-get install -y libgsl-dev && apt-get clean

RUN ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.23.1.0 /usr/lib/x86_64-linux-gnu/libgsl.so.0