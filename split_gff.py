import sys
import gzip

"""
Extract regions from gff and create individual bed files

constituative exons
alternative exons
transcripts

"""


input_gff = sys.argv[1]


def parse_description(d):
    content = {}
    fields = d.rstrip().split(";")
    for f in fields:
        n, v = f.split("=")
        content[n] = v
    return(content)

pfx = input_gff[0:-8]
print(pfx)

with gzip.open(input_gff, 'rb') as gff, \
     open(f"{pfx}.constituative_exon.gff", 'w') as con, \
     open(f"{pfx}.alternative_exon.gff", 'w') as alt, \
     open(f"{pfx}.transcript.gff", 'w') as txs:
    for g_line in gff:
        g_line = g_line.decode('utf-8')
        if (g_line[0] != '#'):
            g_interval = g_line.rstrip().split("\t")
            chromosome = g_interval[0]
            region = g_interval[2]
            start = int(g_interval[3])
            end = int(g_interval[4])
            anno = g_interval[8]
            anno_fields = parse_description(anno)
            suffix = ''
            desc_dict = parse_description(anno)
            #print(desc_dict)
            if region == "exon":
                if desc_dict['constitutive'] == '1':
                    print(g_line.rstrip(), file=con)
                elif desc_dict['constitutive'] == '0':
                    print(g_line.rstrip(), file=alt)
            elif region == "mRNA":
                print(g_line.rstrip(), file=txs)