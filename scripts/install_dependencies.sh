#!/bin/bash

# Initialization
apt-get update
apt-get upgrade -y
apt-get install -y wget

# Main OrthoFinder Software
wget --quiet https://github.com/davidemms/OrthoFinder/releases/download/2.5.2/OrthoFinder_source.tar.gz
tar -zxf OrthoFinder_source.tar.gz
mv OrthoFinder_source /usr/bin/orthofinder

# OrthoFinder dependencies
DEBIAN_FRONTEND=noninteractive apt-get install -y mcl
apt-get install -y mafft
apt-get install -y fasttree
wget --quiet https://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-linux64.tar.gz
tar -zxf diamond-linux64.tar.gz diamond
mv diamond /usr/bin/diamond

# Dependcies for generating plot for main report
DEBIAN_FRONTEND=noninteractive apt-get install -y grace
# Using a fork with fixes for Python 3
git clone https://github.com/samseaver/pygrace.git
cd pygrace
mv PyGrace /opt/conda3/lib/python3.8/site-packages/PyGrace

# Loading PlantSEED data
git clone -b kbase_release https://github.com/ModelSEED/PlantSEED /kb/module/PlantSEED

# Dependencies for Bokeh JS Plots
conda install -c conda-forge phantomjs
pip install --upgrade pip
pip install scipy
pip install bokeh
pip install selenium
