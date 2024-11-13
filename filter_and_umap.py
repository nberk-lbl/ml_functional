
# run from inside gpn/gpn
import gpn.model
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import torch
from transformers import AutoModel, AutoModelForMaskedLM, AutoTokenizer


import json

from Bio import SeqIO
import gzip
import sys
import os
os.environ['HF_HOME'] = '/pscratch/sd/n/nberk/gpn/gpn/tmp/'
run_cfg = json.loads(open(sys.argv[1], 'r').read())

def rc(seq):
    t = {
        'A' : 'T',
        'T' : 'A',
        'G' : 'C',
        'C' : 'G',
    }

    rc = ""
    for i in seq[::-1]:
        n = i
        if i in t:
            n = t[i]
        rc += n
    return(rc)

print('loading_fasta')
seq = ''
seq_len = 0
with gzip.open(run_cfg['fa_path'], "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id == run_cfg['test_seq']:
            seq = str(record.seq)
            seq_len = len(seq)


print("marking repeats")
repeat_bins = set()
with open(run_cfg["rep_bed"]) as repeat_bed:
    for b in repeat_bed:
        fields = b.rstrip().split()
        if fields[0] == run_cfg['test_seq']:
            start, end = int(fields[1]), int(fields[2])
            for i in range(start, end):
                r = int(i/run_cfg['window_size'])
                repeat_bins.add(r)

print("building region map")
nt_region_map = {}
with gzip.open(run_cfg["gff_path"], "rb") as gff:
    for f in gff:
        g = f.decode('utf-8')
        if not(g[0] == "#"): # skip header
            fields = g.rstrip().split()
            if (fields[0] == run_cfg['test_seq']):
                start, end = int(fields[3]), int(fields[4])
                region = fields[2]
                if (region not in ('chromosome')):
                    if not region in nt_region_map:
                        nt_region_map[region] = [0] * 18585056
                    for i in range(start, end):
                        nt_region_map[region][i] = 1

print("getting annotation bins")
region_bins = {}
ambiguous_bins = set() # bins on any annotation boundary

for k in nt_region_map:
    region_bins[k] = set()
    for i in range(0, int(seq_len/run_cfg['window_size'])):
        start = i * run_cfg['window_size']
        end = start + run_cfg['window_size']
        window = nt_region_map[k][start:end]
        if (sum(window) == len(window)):
            region_bins[k].add(i)
        elif (sum(window) > 1 and sum(window) < len(window)):
            ambiguous_bins.add(i)
regions = region_bins.keys()

print("running gpn")
model_path = "songlab/gpn-brassicales"

print("setting up tokenizer")
tokenizer = AutoTokenizer.from_pretrained(model_path)
tokenizer.get_vocab()

print("loading model")
model = AutoModel.from_pretrained(model_path)
model.eval()

final_averaged_embeddings = pd.DataFrame()

#for i in range(0, int(seq_len/run_cfg['window_size'])):
prefix =  f"{run_cfg['out_pfx']}_{run_cfg['test_seq']}"
with open(f"{prefix}_regions.csv", "w") as regions_file:
    print(",".join(regions))
    for i in range(int(seq_len / run_cfg['window_size'])):
        start = i * run_cfg['window_size']
        end = start + run_cfg['window_size']
        window = seq[start:end]
        if not i in repeat_bins and i not in ambiguous_bins and window.count("N") < 5:

            # apply the bin filters and track the region types
            reg = []
            for r in regions:
                x = "0"
                if i in region_bins[r]:
                    x = "1"
                reg.append(x)
            print(",".join(reg), file=regions_file)

            # fwd strand
            print(f"tokenizing {i} : {window}")
            input_ids = tokenizer(str(window), return_tensors="pt", return_attention_mask=False, return_token_type_ids=False)["input_ids"]
            #print(input_ids.shape)

            with torch.no_grad():
                embedding = model(input_ids=input_ids).last_hidden_state
            #print(embedding.shape)

            embedding_df = pd.DataFrame(StandardScaler().fit_transform(embedding[0].numpy()))
            #print(embedding_df.shape)

            averaged_embeddings = embedding_df.groupby(embedding_df.index // run_cfg['window_size']).mean()

            # rev strand
            rev_comp = rc(window)
            print(f"tokenizing {i} : {rev_comp}")
            rc_input_ids = tokenizer(str(rev_comp), return_tensors="pt", return_attention_mask=False, return_token_type_ids=False)["input_ids"]
            #print(rc_input_ids.shape)

            with torch.no_grad():
                rc_embedding = model(input_ids=rc_input_ids).last_hidden_state
            #print(rc_embedding.shape)

            rc_embedding_df = pd.DataFrame(StandardScaler().fit_transform(embedding[0].numpy()))
            #print(rc_embedding_df.shape)

            rc_averaged_embeddings = rc_embedding_df.groupby(rc_embedding_df.index // run_cfg['window_size']).mean()

            # avg
            strand_avg_embeddings = (averaged_embeddings + rc_averaged_embeddings) / 2
            final_averaged_embeddings = pd.concat([final_averaged_embeddings, strand_avg_embeddings], ignore_index=True)

final_averaged_embeddings.to_csv(f'{run_cfg["output_dir"]}/{prefix}_avg_embeddings.tsv', sep='\t', index=False)
