
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
fasta = run_cfg['']

model_path = "songlab/gpn-brassicales"
tokenizer = AutoTokenizer.from_pretrained(model_path)
tokenizer.get_vocab()

print("loading model")
model = AutoModel.from_pretrained(model_path)
model.eval()


"""
window_size = run_cfg['window_size']
seq_len = len(seq)
num_windows = int(seq_len / window_size)

# Tokenize the entire sequence once for both strands
input_ids = tokenizer(seq, return_tensors="pt", return_attention_mask=False, return_token_type_ids=False)["input_ids"]
rev_comp = rc(seq)  # Assuming `rc` is the function to generate reverse complement
rc_input_ids = tokenizer(rev_comp, return_tensors="pt", return_attention_mask=False, return_token_type_ids=False)["input_ids"]

# Get embeddings for forward and reverse strand
with torch.no_grad():
    embedding = model(input_ids=input_ids).last_hidden_state.squeeze(0)  # Shape: (seq_len, 512)
    rc_embedding = model(input_ids=rc_input_ids).last_hidden_state.squeeze(0)  # Shape: (seq_len, 512)

# Standardize embeddings for both strands
embedding = StandardScaler().fit_transform(embedding.numpy())
rc_embedding = StandardScaler().fit_transform(rc_embedding.numpy())

prefix = run_cfg['out_pfx']
final_averaged_embeddings = []
with open(f"{prefix}_regions.csv", "w") as regions_file:
    print(",".join(regions))
    for i in range(num_windows):
        start = i * window_size
        end = start + window_size
        window = seq[start:end]

        # Apply filters
        if i in repeat_bins or i in ambiguous_bins or window.count("N") >= 5:
            continue

        # Track region types
        reg = [("1" if i in region_bins[r] else "0") for r in regions]
        print(",".join(reg), file=regions_file)

        # Forward strand
        fwd_window_embeddings = pd.DataFrame(embedding[start:end]).groupby(lambda idx: idx // window_size).mean()

        # Reverse strand
        rc_window_embeddings = pd.DataFrame(rc_embedding[start:end]).groupby(lambda idx: idx // window_size).mean()

        # Average forward and reverse strand embeddings
        strand_avg_embeddings = (fwd_window_embeddings + rc_window_embeddings) / 2

        # Append to final list
        final_averaged_embeddings.append(strand_avg_embeddings)

    # Concatenate all windows into a single DataFrame
    final_averaged_embeddings = pd.concat(final_averaged_embeddings, ignore_index=True)

final_averaged_embeddings.to_csv(f'{run_cfg["output_dir"]}/{prefix}_avg_embeddings.tsv', sep='\t', index=False)
"""
