import sys
import json

embeddings_file = sys.argv[1]
annot_file = sys.argv[2]

out_embed_name = embeddings_file.split("\.")[0] + "_filtered.txt"
out_annot_name = annot_file.split("\.")[0] + "_filtered.txt"


unannotated_windows = []

with open(embeddings_file) as embf, open(annot_file) as annotf:
    embedding_line = embf.readline()
    annot_line = annotf.readline()

    while embedding_line:
     
