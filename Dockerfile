FROM rpetit3/nextconda-base

MAINTAINER robbie.petit@gmail.com

RUN apt-get -qq update && apt-get -qq -y install g++ gcc

# SPAdes and BBmap
RUN conda install -y spades==3.11.1 \
    && conda install -y bbmap==37.62 \
    && pip install numpy

# FASTQ tools
COPY src/fastq-interleave.cpp /tmp/fastq-interleave.cpp
COPY src/fastq-stats.cpp /tmp/fastq-stats.cpp
RUN cd /tmp && \
    g++ -Wall -O3 -o /usr/local/bin/fastq-interleave fastq-interleave.cpp && \
    g++ -Wall -O3 -o /usr/local/bin/fastq-stats fastq-stats.cpp

COPY src/illumina-cleanup.nf /usr/local/bin/illumina-cleanup.nf
COPY src/illumina-cleanup-noecc.nf /usr/local/bin/illumina-cleanup-noecc.nf
COPY src/illumina-cleanup.py /usr/local/bin/illumina-cleanup.py
COPY src/merge-json.py /usr/local/bin/merge-json.py

RUN mkdir -p /opt/references/ /data
COPY data/phiX-NC_001422.fasta /opt/references/phiX-NC_001422.fasta
COPY data/adapters.fasta /opt/references/adapters.fasta

WORKDIR /data

CMD ["illumina-cleanup.nf", "--help"]
