FROM kbase/sdkpython:3.8.0
MAINTAINER KBase Developer

ENV DEBIAN_FRONTEND=noninteractive

# Tools & build deps (unzip & ca-certs are important for the download step)
RUN apt-get -y update && apt-get -y install --no-install-recommends \
      wget ca-certificates unzip g++ zlib1g-dev make \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir coverage

WORKDIR /kb/module

# --- Bowtie2 2.5.4 (prebuilt linux x86_64) ---
ENV BT2_VER=2.5.4
RUN mkdir -p /opt /kb/module/bowtie2-bin \
 && wget -qO /tmp/bowtie2-${BT2_VER}-linux-x86_64.zip \
      "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BT2_VER}/bowtie2-${BT2_VER}-linux-x86_64.zip/download" \
 && unzip -q /tmp/bowtie2-${BT2_VER}-linux-x86_64.zip -d /opt \
 && ln -s /opt/bowtie2-${BT2_VER}-linux-x86_64 /opt/bowtie2 \
 # global PATH convenience
 && ln -s /opt/bowtie2/bowtie2 /usr/local/bin/bowtie2 \
 && ln -s /opt/bowtie2/bowtie2-build /usr/local/bin/bowtie2-build \
 && ln -s /opt/bowtie2/bowtie2-inspect /usr/local/bin/bowtie2-inspect \
 # keep your historical path used by the app
 && ln -s /opt/bowtie2/bowtie2 /kb/module/bowtie2-bin/bowtie2 \
 && ln -s /opt/bowtie2/bowtie2-build /kb/module/bowtie2-bin/bowtie2-build \
 && ln -s /opt/bowtie2/bowtie2-inspect /kb/module/bowtie2-bin/bowtie2-inspect \
 && rm -f /tmp/bowtie2-${BT2_VER}-linux-x86_64.zip \
 && bowtie2 --version && bowtie2-build --version
# --- end Bowtie2 ---

COPY ./ /kb/module
RUN mkdir -p /kb/module/work && chmod -R a+rw /kb/module

WORKDIR /kb/module
RUN make all

ENTRYPOINT ["./scripts/entrypoint.sh"]
CMD []
