## Using GPN for functional predicion

This paper uses a machine learning model to predict genomic properties.

https://doi.org/10.1073/pnas.2311219120

One application is detecting genomic features e.g. exons, introns.

##### setup

clone this repo https://github.com/songlab-cal/gpn

download fa and gff3 from Phytozome

set this to use scratch for intermediate files (managed by HuggingFace)

export HF_HOME=/path/to/new/cache_diri

conda environment: gpn_conda.yaml
