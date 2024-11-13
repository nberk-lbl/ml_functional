# Trying to replicate the umap plot with a poplar genome

embeddings = read.csv("data/averaged_embeddings.tsv", sep ="\t", head=F)
embeddings.umap = umap(embeddings)

gff_regions = read.csv("datachr19_region_windows.tsv", sep="\t", head=F)
#[1] "gene..."              "gene....1"            "mRNA..."              "mRNA....1"            "exon..."             
# [6] "exon....1"            "CDS..."               "CDS....1"             "five_prime_UTR..."    "five_prime_UTR....1" 
#[11] "three_prime_UTR..."   "three_prime_UTR....1"

colnames(gff_regions) = c(
    'gene.plus', 'gene.minus', 
    'mRNA.plus', 'mRNA.minus', 
    'exon.plus', 'exon.minus', 
    'CDS.plus', 'CDS.minus', 
    'five_prime_UTR.plus', 'five_prime_UTR.minus', 
    'three_prime_UTR.plus', 'three_prime_UTR.minus', 
    )

r = 3
c1 <- ifelse(gff_regions[,3] < 0.75, 0, 2)
c2 <- ifelse(gff_regions[,4] < 0.75, 0, 3)
c3 <- c1 + c2 + 1

plot(m$layout[,1], m$layout[,2], cex=.4, pch=19, col=c, xlim=c(-10,10), ylim=c(-10,10))

# This plot doesn't show feture separation. There is some structure but it doesn't line up with annotations
##  the GFF has fewer feature types
##  the Song lab paper filtered out repeats
##  also filtered out any "ambiguous" bins
##  the training genomes were related to the test genome.

# Try arabidopsis
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/

embeddings = read.csv("data/averaged_embeddings.tsv", sep ="\t", head=F)
embeddings.umap = umap(embeddings)



