
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

print("loading model")
model_path = "songlab/gpn-brassicales"
model = AutoModelForMaskedLM.from_pretrained(model_path)

print("loading fasta")
sequences = {}
with gzip.open(f'{run_cfg["fa_path"]}','rt') as genome_fa:
    for record in SeqIO.parse(genome_fa, "fasta"):
        sequences[record.id] = record.seq

test_name = run_cfg['test_seq']
print(f"testing with {test_name}")

plt.figure(figsize=(10, 6))


print("setting up tokenizer")
tokenizer = AutoTokenizer.from_pretrained(model_path)
tokenizer.get_vocab()

print("loading model")
model = AutoModel.from_pretrained(model_path)
model.eval()

window_size = run_cfg['window_size']
chunk_size = window_size * 1000
pfx = run_cfg['out_pfx']

for rec in sequences.keys():
    seq = sequences[rec]
    final_averaged_embeddings = pd.DataFrame()
    for i in range(int(len(seq) / chunk_size) + 1):
        a = i * chunk_size
        b = min((i+1) * chunk_size, len(seq))
        
        print(f"computing {rec}: {a} - {b}")
        chunk = seq[a:b]
        input_ids = tokenizer(str(chunk), return_tensors="pt", return_attention_mask=False, return_token_type_ids=False)["input_ids"]
        print(input_ids.shape)

        with torch.no_grad():
            embedding = model(input_ids=input_ids).last_hidden_state
        print(embedding.shape)

        embedding_df = pd.DataFrame(StandardScaler().fit_transform(embedding[0].numpy()))
        averaged_embeddings = embedding_df.groupby(embedding_df.index // window_size).mean()
        final_averaged_embeddings = pd.concat([final_averaged_embeddings, averaged_embeddings], ignore_index=True)

    final_averaged_embeddings.to_csv(f'{run_cfg["output_dir"]}/{pfx}_{rec}_avg_embeddings.tsv', sep='\t', index=False)
    