import sys
import gzip
import json

run_cfg = json.loads(open(sys.argv[1], 'r').read())
window_size = run_cfg['window_size']
threshold = run_cfg['threshold']
outfile = sys.argv[1][0:-4] + "_gff_windows.tsv"

main_gff = run_cfg['gff_path']
intron_gff = run_cfg['intron_gff']

target_chrom = run_cfg["test_seq"]

windows = {} # windows by chromosome and start position
tags = run_cfg["gff_tags"]

def handle_fields(g_line):
    global max_window
    chrom = fields[0]
    start = int(fields[3]) - 1
    end = int(fields[4]) - 1

    # if a specific chrom/contig is specified, only get those features. otherwise use all chromosomes
    if chrom == target_chrom or target_chrom == "":
        label = fields[2]
        if label in tags:
            if not chrom in windows:
                windows[chrom] = {}
            
            for nt in range(start, end):
                
                window_num = int(nt / window_size)
                if not window_num in windows[chrom]:
                    windows[chrom][window_num] = {}
                if not label in windows[chrom][window_num]:
                    windows[chrom][window_num][label] = [0] * window_size
                
                pos = nt % window_size
                windows[chrom][window_num][label][pos] = 1

## parse the tagged features from the main gff
# the main gff does not contain introns or repeats

with gzip.open(main_gff, 'rb') as gff:
    for gff_line in gff:
        g_line = gff_line.decode("utf-8").rstrip()      
        if g_line[0] != '#':
            fields = g_line.split("\t")
            handle_fields(fields)

## add introns
tags = ['intron']
with open(intron_gff) as i_gff:
    for g in i_gff:
        g_line = g.rstrip()
        fields = g_line.split("\t")
        fields[2] = 'intron'
        handle_fields(fields)

## add repeats
tags = ['repeat']
with open(run_cfg["repeat_bed"]) as bed:
    for bed_line in bed:
        chrom, start, end, name, score, strand = bed_line.rstrip().split("\t")
        fields = [chrom[3:], '', 'repeat', start, end]
        handle_fields(fields)

tags = run_cfg["gff_tags"] + ['intron', 'repeat']
for chrom in windows:
    if (chrom in run_cfg['chrom_size']):
        max_window = int(run_cfg['chrom_size'][chrom] / window_size) + 1
        for b in range(max_window):
            assigned = False
            if b in windows[chrom]:
                p_fields = set()
                p_sum = 0
                for t in tags:
                    if t in (windows[chrom][b].keys()):
                        s = sum(windows[chrom][b][t])
                        p_sum + s
                        if (s > run_cfg['threshold']):
                            p_fields.add(t)

                assigned = True
                if len(p_fields) == 1 and p_sum < 110:
                    print(f"{chrom},{b},{p_fields.pop()}")
                else:
                    print(f"{chrom},{b},ambiguous")

            if not assigned:
                print(f"{chrom},{b},intergenic")


