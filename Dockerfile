FROM ubuntu:16.04
MAINTAINER robbie.petit@gmail.com

# Bioconda (SPAdes, BBmap)
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 g++ gcc \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && apt-get -qq -y remove bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

RUN conda config --add channels r && \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -y spades==3.11.1 && \
    conda install -y bbmap==37.62

RUN pip install numpy

ENV PATH /opt/conda/bin:$PATH

# Nextflow
RUN cd /tmp && \
    curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin && \
    chmod 755 /usr/local/bin/nextflow

# FASTQ tools
COPY src/fastq-interleave.cpp /tmp/fastq-interleave.cpp
COPY src/fastq-stats.cpp /tmp/fastq-stats.cpp
RUN cd /tmp && \
    g++ -Wall -O3 -o /usr/local/bin/fastq-interleave fastq-interleave.cpp && \
    g++ -Wall -O3 -o /usr/local/bin/fastq-stats fastq-stats.cpp

COPY src/illumina-cleanup.nf /usr/local/bin/illumina-cleanup.nf
COPY src/illumina-cleanup.py /usr/local/bin/illumina-cleanup.py
COPY src/merge-json.py /usr/local/bin/merge-json.py

RUN mkdir -p /opt/references/ /data
COPY data/phiX-NC_001422.fasta /opt/references/phiX-NC_001422.fasta
COPY data/adapters.fasta /opt/references/adapters.fasta

WORKDIR /data

CMD ["illumina-cleanup.nf", "--help"]
