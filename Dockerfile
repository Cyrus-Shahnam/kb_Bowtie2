FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------

# RUN apt-get update

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install coverage && apt-get -y update && apt-get -y install wget g++ zlib1g-dev

WORKDIR /kb/module

RUN VERSION='2.3.2' \
    && mkdir bowtie2-bin \
    && wget --no-verbose "http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${VERSION}/bowtie2-${VERSION}-source.zip" \
    && unzip -q bowtie2-${VERSION}-source.zip \
    && cd bowtie2-${VERSION} \
    && make NO_TBB=1 \
    && cp bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s \
          bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s ../bowtie2-bin \
    && cd .. \
    && rm -rf bowtie2-${VERSION}*


COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
