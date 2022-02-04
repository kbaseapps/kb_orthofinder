FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && apt-get upgrade -y && apt-get install -y wget

RUN conda install -c conda-forge phantomjs

RUN pip install --upgrade pip && \
    pip install scipy && \
    pip install bokeh && \
    pip install selenium

RUN wget --quiet https://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-linux64.tar.gz && \
    tar -zxf diamond-linux64.tar.gz diamond && \
    mv diamond /usr/bin/diamond && \
    diamond version

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y grace && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y mcl

RUN apt-get install -y mafft && \
    mafft --version
    
RUN apt-get install -y fasttree && \
    ln -s /usr/bin/fasttree /usr/bin/FastTree && \
    FastTree -expert

RUN mkdir -p /kb/deps
WORKDIR /kb/deps

RUN wget --quiet https://github.com/davidemms/OrthoFinder/releases/download/2.5.2/OrthoFinder_source.tar.gz && \
    tar -zxf OrthoFinder_source.tar.gz && \
    mv OrthoFinder_source /kb/deployment/bin/orthofinder

# Using a fork with fixes for Python 3
RUN git clone https://github.com/samseaver/pygrace.git && \
    cd pygrace && \
    mv PyGrace /kb/deployment/lib/PyGrace

# Loading PlantSEED data
RUN git clone -b kbase_release https://github.com/ModelSEED/PlantSEED /kb/module/PlantSEED

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
