FROM quay.io/pacbio/pbmm2@sha256:b58eef0645dd9adca850c9b4811f289ddcd57f630596d9f524701bc787467f30

LABEL MAINTAINER="Bioinformatics <bioinformatics@dnastack.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME="${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG="${IMAGE_TAG}"

RUN apt-get -qq update && \
	apt-get -qq install \
        wget

ARG JASMINE_VERSION
RUN wget https://github.com/PacificBiosciences/jasmine/releases/download/v${JASMINE_VERSION}/jasmine.tar.xz && \
    tar xf jasmine.tar.xz -C /usr/local && \
    chown root:root /usr/local/bin/jasmine && \
    chown root:root /usr/local/lib/libonnxruntime.so.* && \
    chmod +x /usr/local/bin/jasmine
