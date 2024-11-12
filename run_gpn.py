
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

"""
Try visualize this gene

Chr19   phytozomev13    gene    688354  695128  .       +       .       ID=Potri.019G000800.v4.1;Name=Potri.019G000800;ancestorIdentifier=Potri.019G000800.v3.1
Chr19   phytozomev13    mRNA    688354  695128  .       +       .       ID=Potri.019G000800.1.v4.1;Name=Potri.019G000800.1;pacid=42774537;longest=1;ancestorIdentifier=Potri.019G000800.1.v3.1;Parent=Potri.019G000800.v4.1
Chr19   phytozomev13    exon    688354  689109  .       +       .       ID=Potri.019G000800.1.v4.1.exon.1;Parent=Potri.019G000800.1.v4.1;pacid=42774537
Chr19   phytozomev13    five_prime_UTR  688354  688446  .       +       .       ID=Potri.019G000800.1.v4.1.five_prime_UTR.1;Parent=Potri.019G000800.1.v4.1;pacid=42774537
Chr19   phytozomev13    CDS     688447  689109  .       +       0       ID=Potri.019G000800.1.v4.1.CDS.1;Parent=Potri.019G000800.1.v4.1;pacid=42774537
Chr19   phytozomev13    exon    690108  690393  .       +       .       ID=Potri.019G000800.1.v4.1.exon.2;Parent=Potri.019G000800.1.v4.1;pacid=42774537
Chr19   phytozomev13    CDS     690108  690393  .       +       0       ID=Potri.019G000800.1.v4.1.CDS.2;Parent=Potri.019G000800.1.v4.1;pacid=42774537
Chr19   phytozomev13    exon    694420  695128  .       +       .       ID=Potri.019G000800.1.v4.1.exon.3;Parent=Potri.019G000800.1.v4.1;pacid=42774537

another one to try - putative Rubisco methyltransferase
Chr14:12353883..12358380 forward
"""

print("setting up tokenizer")
tokenizer = AutoTokenizer.from_pretrained(model_path)
tokenizer.get_vocab()

print("loading model")
model = AutoModel.from_pretrained(model_path)
model.eval()

final_averaged_embeddings = pd.DataFrame()
chunk_size = int(1e6)
i = 0

while (i * chunk_size < len(sequences[test_name])):
    a = i * chunk_size
    b = min((i+1) * chunk_size, len(sequences[test_name]))
    print(f"computing {a} - {b}")

    seq = sequences[test_name][a:b]

    print(f"tokenizing {test_name}")
    input_ids = tokenizer(str(seq), return_tensors="pt", return_attention_mask=False, return_token_type_ids=False)["input_ids"]
    print(input_ids.shape)


    with torch.no_grad():
        embedding = model(input_ids=input_ids).last_hidden_state
    print(embedding.shape)

    # Calculate the mean for each 100-row window
    window_size = run_cfg['window_size']
    embedding_df = pd.DataFrame(StandardScaler().fit_transform(embedding[0].numpy()))
    averaged_embeddings = embedding_df.groupby(embedding_df.index // window_size).mean()
    final_averaged_embeddings = pd.concat([final_averaged_embeddings, averaged_embeddings], ignore_index=True)

    i+=1

pfx = sys.argv[1][0:-4]
averaged_embeddings.to_csv(f'{run_cfg["output_dir"]}/{pfx}_avg_embeddings.tsv', sep='\t', index=False)
    