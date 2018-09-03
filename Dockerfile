FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y mcl && \
    apt-get install -y mafft && \
    apt-get install -y fasttree && \
    ln -s /usr/bin/fasttree /usr/bin/FastTree

RUN which mafft && \
    mafft --version 2>&1 | grep MAFFT && \
    which FastTree && \
    FastTree -expert 2>&1 | grep Detailed

RUN mkdir -p /kb/deps
WORKDIR /kb/deps
RUN wget --quiet https://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz && \
    tar -zxvf diamond-linux64.tar.gz && \
    mv diamond /kb/deployment/bin/diamond

RUN which diamond && \
    diamond version

RUN wget --quiet https://github.com/davidemms/OrthoFinder/releases/download/v2.2.6/OrthoFinder-2.2.6_source.tar.gz && \
    tar -zxf OrthoFinder-2.2.6_source.tar.gz && \
    mv OrthoFinder-2.2.6_source/orthofinder /kb/deployment/bin/orthofinder

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
