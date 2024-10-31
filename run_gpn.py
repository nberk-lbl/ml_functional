
# run from inside gpn/gpn
import gpn.model
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import torch
from transformers import AutoModel, AutoModelForMaskedLM, AutoTokenizer

from Bio import SeqIO
import gzip
import os
os.environ['HF_HOME'] = '/pscratch/sd/n/nberk/gpn/gpn/tmp/'

print("loading model")
model_path = "songlab/gpn-brassicales"
model = AutoModelForMaskedLM.from_pretrained(model_path)

data_path = "../Phytozome"
output_dir = "../output"

print("loading fasta")
sequences = {}
with gzip.open(f'{data_path}/PhytozomeV13/Ptrichocarpa/v4.1/assembly/Ptrichocarpa_533_v4.0.fa.gz','rt') as genome_fa:
    for record in SeqIO.parse(genome_fa, "fasta"):
        sequences[record.id] = record.seq

test_name = 'Chr19'
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

seq = sequences[test_name][0:1000000]

print("setting up tokenizer")
tokenizer = AutoTokenizer.from_pretrained(model_path)
tokenizer.get_vocab()

print(f"tokenizing {test_name}")
input_ids = tokenizer(str(seq), return_tensors="pt", return_attention_mask=False, return_token_type_ids=False)["input_ids"]
print(input_ids.shape)
input_ids[0, :10]

print("embedding")
model = AutoModel.from_pretrained(model_path)
model.eval()


with torch.no_grad():
    embedding = model(input_ids=input_ids).last_hidden_state
print(embedding.shape)

print("creating visualization")
embedding_df = pd.DataFrame(StandardScaler().fit_transform(embedding[0].numpy()))
embedding_df.index.name = "Position"
embedding_df.columns.name = "Embedding dimension"

print("embedding df")
print(embedding_df)

plt.figure(figsize=(10, 6))
sns.heatmap(embedding_df.T, center=0, vmin=-3, vmax=3, cmap="coolwarm", square=True, xticklabels=100, yticklabels=100)
plt.savefig(f'{output_dir}/{test_name}.png')
