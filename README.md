## Using GPN for functional predicion

This paper uses a machine learning model to predict genomic properties.

https://doi.org/10.1073/pnas.2311219120

One application is detecting genomic features e.g. exons, introns.

#### setup

clone this repo https://github.com/songlab-cal/gpn

set this to use scratch for intermediate files (managed by HuggingFace)

`export HF_HOME=$PSCRATCH`
conda environment: `gpn_conda.yaml`


#### genomes

##### arabidopsis from ENSEMBL
https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gff3/arabidopsis_thaliana/

repeats - this site also has poplar
http://ucsc.gao-lab.org/cgi-bin/hgTables?hgsid=167291_E9nY5UIAQRUOAR01xJAsum4vDukw

introns are transcripts minus exons (this will include the rare introns that occur in UTRs)
`zgrep exon Arabidopsis_thaliana.TAIR10.60.chromosome.4.gff3.gz > Arabidopsis_thaliana.TAIR10.60.chromosome.4.exon.gff3`
`zgrep mRNA Arabidopsis_thaliana.TAIR10.60.chromosome.4.gff3.gz > Arabidopsis_thaliana.TAIR10.60.chromosome.4.transcript.gff3`
`bedtools subtract -a Arabidopsis_thaliana.TAIR10.60.chromosome.4.transcript.gff3 -b Arabidopsis_thaliana.TAIR10.60.chromosome.4.exon.gff3 > Arabidopsis_thaliana.TAIR10.60.chromosome.4.intron.gff3`

give each bin a single region label if possible. Ignore bins that stradle regions or have conflicting regions on different strands
`python gff_to_window.py TAIR10.chr4.json > TAIR10.chr4.anno.csv`

#### steps

get average embeddings for each 100nt window
`python run_gpn.py TAIR10.chr4.json`

make umap plot colored by labels

#### result
Putting the labels on the umap doesn't label cluseters. 

<image>

#### TODO
Try with an arabidopsis chromosome
- get the chrom size ahead of time so the files are the same dimension
- limit the number of unlabeled bins to ~ the number of labeled bins
- discard ambiguous bins (see paper)

