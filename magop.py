#!/usr/bin/env python3

import argparse
import tempfile
from magnumopus import map_reads_to_ref
import magnumopus

def parse_args():
    parser = argparse.ArgumentParser(description="Process assembly, primer, read, and reference files")
    
    # Optional argument for assembly files
    parser.add_argument(
        '-a', '--assembly', 
        nargs='*',  
        help='List of assembly files'
    )
    
    # Required argument for primer file
    parser.add_argument(
        '-p', '--primers', 
        required=True, 
        help='Primer file'
    )
    
    # Optional argument for read files
    parser.add_argument(
        '-r', '--reads', 
        nargs='*',  
        help='List of read files'
    )
    
    # Optional argument for reference sequences file
    parser.add_argument(
        '-s', '--ref_seqs', 
        help='Reference sequences file'
    )
    
    return parser.parse_args()




#do reads mapping for the reads and return a temporary file
def reads_mapping(forward_read, reverse_read, ref_seqs):
    sam1 = map_reads_to_ref(
        ref=ref_seqs,
        r1=forward_read, 
        r2=reverse_read   
    )
    mapped_seq = sam1.best_consensus(fasta=True)
    
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_file:
        temp_file.write(mapped_seq.encode("utf-8"))
        mapped_seq_temp_file = temp_file.name
    
    return mapped_seq_temp_file



#for reads, do isPCR to find amplicons
def is_reads():
    args = parse_args()
    reads_dict={}

    #assign forward and reverse reads for each _1 and _2 files.
    args.reads = sorted(args.reads)
    for i in range(0, len(args.reads), 2):
        forward_read = args.reads[i]  
        reverse_read = args.reads[i + 1] 

        #read mapping
        mapped_seq_temp_file = reads_mapping(forward_read, reverse_read, args.ref_seqs)

        #run isPCR to get 16S amplicon
        amplicons = magnumopus.ispcr(args.primers, mapped_seq_temp_file)

        #store species name and amplicon seq in dict
        contents = amplicons.split('\n')
        sequence=contents[1]    

        with open(forward_read, "r") as forward_read_file:
            for line in forward_read_file:
                if line.startswith("@"):
                    species_classification=line.split('.')
                    header=species_classification[0].lstrip("@")
    
        reads_dict[header]=sequence

    return reads_dict

#for assembly, do isPCR to find amplicons
def is_assembly():
    args = parse_args()
    assembly_dict={}

    #run isPCR to get 16S amplicon
    for assembly in args.assembly:

        amplicons = magnumopus.ispcr(args.primers, assembly)

        #get just the first amplicon in output since there are multiple
        amplicon_list = amplicons.split(">")
        first_amplicon = ">" + amplicon_list[1] 

        #store species name and amplicon seq in dict
        contents = first_amplicon.split('\n')
        sequence=contents[1]        

        species_classification=contents[0].split(':')
        header=species_classification[0].lstrip(">")
        assembly_dict[header]=sequence

    return assembly_dict


#for reference, do isPCR to find amplicons
def is_ref():
    args = parse_args()
    ref_dict={}

    with open(args.ref_seqs, 'r') as ref_file:
            ref_contents = ref_file.read()
            entries = ref_contents.strip().split('>')
            entries = [entry for entry in entries if entry]
            

            for entry in entries:
                #make temp file for each species with header and sequence
                lines = entry.strip().split('\n')
                header = lines[0] 
                sequence = '\n'.join(lines[1:])  

                temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fna')
                temp_file.write(f'>{header}\n{sequence}\n')
                temp_file.close() 

                ref_temp_file=temp_file.name

                #run isPCR to get 16S amplicon
                amplicons = magnumopus.ispcr(args.primers, ref_temp_file)


                #store species name and amplicon seq in dict
                contents = amplicons.split('\n')
                sequence=contents[1]        

                species_classification=contents[0].split(':')
                header=species_classification[0].lstrip(">")
                ref_dict[header]=sequence

    return ref_dict


#do needleman-wunsch to determine best orientation 
def do_newick():
    #combine all dict into one 
    ref_amplicon_dict=is_ref()
    reads_amplicon_dict=is_reads()
    assembly_amplicon_dict=is_assembly()
    all_amplicon_dict = {**ref_amplicon_dict, **reads_amplicon_dict, **assembly_amplicon_dict}
    aligned_amplicon_dict={}

    #make the very first sequence listed in the dict as anchor sequence
    anchor_amplicon=next(iter(all_amplicon_dict.values()))

    for species,seq in all_amplicon_dict.items():
        forward_amplicon=seq

        #make reverse compliment sequence
        reverese_amplicon=seq[::-1]
        reverese_comp_amplicon=""
        for base in reverese_amplicon:
            if base=="A":
                reverese_comp_amplicon+="T"
            elif base=="T":
                reverese_comp_amplicon+="A"
            elif base=="C":
                reverese_comp_amplicon+="G"
            elif base=="G":
                reverese_comp_amplicon+="C"
        
        # #use needleman wunsch package
        aln_forward, score_forward = magnumopus.needleman_wunsch(anchor_amplicon, forward_amplicon, 1, -1, -1)
        aln_reverse_comp, score_reverse_comp = magnumopus.needleman_wunsch(anchor_amplicon, reverese_comp_amplicon, 1, -1, -1)

        #add best orientation seq to aligned_amplicon_dict
        if score_forward>score_reverse_comp:
            aligned_amplicon_dict[species]=forward_amplicon
            
        else:
            aligned_amplicon_dict[species]=reverese_comp_amplicon
    return aligned_amplicon_dict


#using the updated dictionary, write a fasta format file with species and amplicon sequences
def amplicon_to_fasta():
    output_file = 'amplicons.fasta'
    aligned_amplicon_dict=do_newick()
    #write a fasta format file
    with open(output_file, 'w') as fasta_file:
        for species, amplicon in aligned_amplicon_dict.items():
            fasta_file.write(f'>{species}\n')
            fasta_file.write(f'{amplicon}\n')

    #open the new fasta format file and print contents
    with open(output_file, 'r') as file:
        contents = file.read()
        print(contents)



if __name__ == "__main__":
    amplicon_to_fasta()
