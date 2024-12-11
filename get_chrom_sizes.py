from Bio import SeqIO
import gzip
import sys 

with gzip.open(sys.argv[1],'rt') as fa:
    for record in SeqIO.parse(fa, "fasta"):
        print(f'"{record.id}" :  {len(record.seq)}')