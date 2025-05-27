#!/usr/bin/env python3

import subprocess
import tempfile

def ispcr(primer_file: str, assembly_file: str) -> str:
    blast= ["blastn", "-query", primer_file, "-subject", assembly_file,
            "-task", "blastn-short", "-outfmt",  "6 std qlen", "-word_size", "6",
            "-penalty", "-2"]
    blast_str=subprocess.run(blast, capture_output=True, text=True)
    lines=blast_str.stdout.strip().split('\n')
    blast_list=[]
    for line in lines:
        col=line.strip().split('\t')
        blast_list.append(col)

    filtered_list=[]
    query_length=[]
    with open(primer_file) as f:
        for line in f:
            if line[0]!=">":
                query_length.append(len(line.strip()))
    for line in blast_list:
        if float(line[2])>=80 and int(line[3]) in query_length:
            filtered_list.append(line)

    sorted_good_hits = sorted(filtered_list, key=lambda x: int(x[8]))


    to_right=[]
    to_left=[]
    split_list=[]
    for line in sorted_good_hits:
        if int(line[8]) < int(line[9]):
            to_right.append(line)
        else:
            to_left.append(line)
    split_list=[to_right,to_left]


    to_right=split_list[0]
    to_left=split_list[1]
    pairedtup_list=[]
    for line_r in to_right:
        for line_l in to_left:
            if int(line_r[8])<int(line_l[9]) :
                atuple=(line_r,line_l)
                pairedtup_list.append(atuple)

    bed_str=""
    for line in pairedtup_list:
        contig=line[0][1]
        primer1=line[0]
        primer2=line[1]
        start=int(primer1[9])
        stop=int(primer2[9])-1
        bed_str += f"{contig}\t{start}\t{stop}\n"

    with tempfile.NamedTemporaryFile(mode='w+') as bed_tmp:
        bed_tmp.write(bed_str)
        bed_tmp.seek(0)
        extract=["seqtk", "subseq", assembly_file, bed_tmp.name]
        results=subprocess.run(extract, capture_output=True, text=True)

    return results.stdout
