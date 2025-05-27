# Magnumopus - 16S rRNA V4 Region Extractor and Aligner

`magop.py` is a command-line tool for extracting and aligning the **V4 region** of the **16S rRNA gene** from assemblies and/or Illumina paired-end read files. It combines PCR simulation (`isPCR.py`), read mapping (`mapping.py`), and Needleman-Wunsch alignment (`nw.py`) into a single streamlined pipeline.

## 🔍 Features

- Extract 16S V4 region amplicons from assembly files using user-provided primers.
- Map paired-end reads to reference sequences and extract consensus V4 regions.
- Align and consistently orient extracted amplicons using Needleman-Wunsch scoring.
- Output final amplicons in FASTA format to `stdout`.

---

## 📦 Requirements

- Python 3.x
- BLAST+ (`blastn`)
- `isPCR` tool (e.g., from UCSC)
- Python libraries:  
  - `argparse`
  - `subprocess`
  - `Bio` (Biopython)

Install Biopython (if not already):

```bash
pip install biopython
```

## 🚀 Usage
```bash
./magop.py [-a ASSEMBLY [ASSEMBLY ...]] -p PRIMERS [-r READS [READS ...]] [-s REF_SEQS]
```

### Arguments
| Argument | Required | Description                                                                                            |
| -------- | -------- | ------------------------------------------------------------------------------------------------------ |
| `-a`     | No       | Space-separated list or glob of assembly `.fasta` files                                                |
| `-p`     | Yes      | Path to primer file in `.fna` format                                                                   |
| `-r`     | No       | Space-separated list or glob of paired-end read files (must follow `*_1.fastq` and `*_2.fastq` naming) |
| `-s`     | No       | Reference `.fna` file for mapping reads                                                                |


## 📥 Input File Format

- Assemblies (-a): FASTA files of genome assemblies.
- Reads (-r): Paired-end FASTQ files with filenames like sampleID_1.fastq and sampleID_2.fastq.
- Primers (-p): Primer sequences in FASTA format (e.g., 515F and 806R primers).
- Reference Sequences (-s): Reference FASTA file with full-length or partial 16S sequences.

## 📤 Output

- Prints aligned, consistently oriented V4 region amplicons to stdout in FASTA format.


## 🧪 Example Usage
./magop.py \
  -a data/assemblies/* \
  -p data/primers/general_16S_515f_806r.fna \
  -r data/reads/Ecoli_1.fastq data/reads/Ecoli_2.fastq \
  -s data/refs/V4.fna


## 📁 Directory Structure
```bash
├── magnumopus
│   ├── __init__.py
│   ├── ispcr.py
│   ├── mapping.py
│   ├── nw.py
│   ├── run_external.py
│   ├── sam.py
├── magop.py
├── data/
│   ├── assemblies/
│   │   └── Escherichia_coli_K12.fna
│   │   └── Pseudomonas_aeruginosa_PAO1.fna
│   │   └── ...
│   ├── reads/
│   │   ├── ERR11767307_1.fastq
│   │   └── ERR11767307_2.fastq
│   │   └── ...
│   ├── primers/
│   │   └── general_16S_515f_806r.fna
│   └── refs/
│       └── V4.fna
└── README.md
```


