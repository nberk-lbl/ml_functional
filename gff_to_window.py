import sys
import gzip

window_size = int(sys.argv[1])
infile = "../gpn/Phytozome/PhytozomeV13/Ptrichocarpa/v4.1/annotation/Ptrichocarpa_533_v4.1.gene_exons.gff3.gz" # sys.argv[2]
target_chrom = "Chr19"
approx_chrom_len = int(16e6 / window_size)
chrom_bins = {}

with gzip.open(infile, 'rb') as gff:
    for gff_line in gff:
        
        fields = gff_line.decode("utf-8").rstrip().split("\t")
        chrom = fields[0]     

        if chrom == target_chrom:
            label = fields[2]
            if not label in chrom_bins:
                chrom_bins[label]= {
                    '+' : [0] * approx_chrom_len,
                    '-' : [0] * approx_chrom_len
                }
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            
            for nt in range(start, end):
                bin_num = int(nt / window_size)
                if bin_num >= len(chrom_bins[label][strand]):
                    print(fields)
                chrom_bins[label][strand][bin_num] += 1

header = []
for label in chrom_bins.keys():
    for strand in chrom_bins[label].keys():
        header.append(f"{label}({strand})")
print("\t".join(header))

for i in range(approx_chrom_len):
    fields = []
    for label in chrom_bins.keys():
        for strand in chrom_bins[label].keys():
            fields.append(f"{chrom_bins[label][strand][i]}")
    print("\t".join(fields))
 

            

