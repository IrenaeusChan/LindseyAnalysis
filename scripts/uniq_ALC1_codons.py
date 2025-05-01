from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import csv
import argparse
import pandas as pd
import numpy as np
import re

# Handle the input arguments
parser = argparse.ArgumentParser(description='Filter merged FASTQ reads and count codon variants directly.')
parser.add_argument("-s", "--string", type=str, required=True, help="Prefix for input and output files")
parser.add_argument("-l", "--library", type=str, required=True, help="Library identifier (A, B, C, or D)")
parser.add_argument("--sequence", type=str, required=True, help="Library DNA sequences")
parser.add_argument("--syn", type=str, required=True, help="Synonymous mutation codon")
parser.add_argument("--wt", type=str, required=True, help="Wildtype codon")
args = parser.parse_args()

lengths = defaultdict(int)

# Normalizes the variant counts to counts per million
def formater(table, length):
    for x in list(table.columns.values[1::]):
        table[x] = pd.to_numeric(table[x], errors='coerce').fillna(0).astype(np.int64)
        table[x] = (table[x])/(length/1000000)
    return table

# Translates incoming codons to the amino acid
def translate_codons(codons):
    amino_acids = []
    for codon in codons:
        codon_seq = Seq(codon)
        amino_acids.append(str(codon_seq.translate()))
    return amino_acids

# Library DNA sequences (synonymously recoded WT)
# Enter information as shown, 1) codon position with syn. mutation, 2) syn.
# mutation, 3) wildtype codon, 4) entire sequence with syn. mutation
syn_codon_pos = 25
syn_codon = args.syn
wt_codon = args.wt
Library = args.sequence
# if args.library == "Library_1": 
#     syn_codon_pos = 25
#     syn_codon = "GCT"
#     wt_codon = "GCC"
#     Library = "GAGCGGGCTGGCGCGACAAGCAGAGGAGGGCAGGCCCCAGGTTTCTTGTTGCGACTCCACACGGAAGGAAGAGCTGAAGCGGCTCGAGTGCAAGAACAAGACCTGAGACAATGGGGCCTGACCGGGATTCACCTCAGGTCCTATCAACTG"
# elif args.library == "Library_2":
#     syn_codon_pos = 25
#     syn_codon = "GGA"
#     wt_codon = "GGC"
#     Library = "GAGGGCGTGAACTGGTTGGCCCAGCGATTCCATTGTCAGAACGGCTGCATCCTCGGTGACGAAATGGGGCTGGGAAAGACATGTCAAACGATCGCGCTTTTTATCTACTTGGCTGGCCGGCTGAATGATGAAGGACCTTTTCTGATTCTT"
# elif args.library == "Library_3":
#     syn_codon_pos = 25
#     syn_codon = "TAT"
#     wt_codon = "TAC"
#     Library = "TGCCCGCTGAGCGTTCTGTCTAATTGGAAAGAGGAGATGCAGCGATTTGCGCCTGGATTGTCATGTGTGACATATGCGGGCGATAAAGAAGAAAGGGCATGTTTGCAGCAAGATCTGAAACAAGAATCCCGGTTCCATGTTTTGCTTACT"
# elif args.library == "Library_4":
#     syn_codon_pos = 25
#     syn_codon = "GCC"
#     wt_codon = "GCG"
#     Library = "ACTTACGAAATTTGCTTGAAAGATGCGAGTTTCCTTAAATCCTTCCCCTGGTCCGTTTTGGTGGTAGACGAAGCCCACAGGCTTAAGAACCAGTCTAGTCTCCTGCATAAAACGCTTAGCGAATTTAGTGTCGTTTTTTCACTTCTCCTT"

# Use re.findall to split the sequence into groups of 3 bases.
# The pattern '.{1,3}' will match 1 to 3 characters, so if the sequence length
# isn't a multiple of 3, the final group might be shorter.
codons = re.findall(r'.{1,3}', Library)

# Format into a dictionary with key "Sequence"
# WT codon sequences for each library (synonymously recoded wildtype)
WT_dict = {args.library: codons}
WT = WT_dict.get(args.library, [])
#print(WT)

if not WT:
    raise ValueError("Invalid library selection. Choose A, B, C, or D.")

# Generate codon counts dictionary dynamically
codons = {}
for i in range(len(WT)):
    codons[f"codon{i + 1}"] = defaultdict(int)

#print(codons)

# Filtering sequences: Only keep sequences with correct length (number of codons * 3)
print('Filtering sequences...')
with open(f"{args.string}.filtered.fasta", "w") as handle:
    for seq_record in SeqIO.parse(f"{args.string}.fasta", "fasta"):
        lengths[len(seq_record)] += 1
        if len(seq_record) == len(WT) * 3:
            SeqIO.write(seq_record, handle, "fasta")


with open(f'{args.string}.lengths.csv', 'w') as f:
    writer = csv.writer(f)
    for row in lengths.items():
        writer.writerow(row)

# Correct reading frame and orient sequences based on the first two codons
sequences = defaultdict(int)
print('Orienting sequences...')
with open(f"{args.string}.reversed.fasta", "w") as handle:
    for seq_record in SeqIO.parse(f"{args.string}.filtered.fasta", "fasta"):
        first_codon = str(seq_record.seq[0:3])
        second_codon = str(seq_record.seq[3:6])
        # If the first two codons do not match the WT, assume wrong orientation and take reverse complement.
        if first_codon != WT[0] and second_codon != WT[1]:
            rc_record = seq_record.reverse_complement(id=True, name=True, description=True)
            SeqIO.write(rc_record, handle, "fasta")
            sequences[str(rc_record.seq)] += 1
        else:
            SeqIO.write(seq_record, handle, "fasta")
            sequences[str(seq_record.seq)] += 1

with open(f'{args.string}.csv', 'w') as f:
    writer = csv.writer(f)
    for row in sequences.items():
        writer.writerow(row)

# Count codon variants directly from nucleotide sequences (without translation)
print('Counting codon variants directly...')

total_count = 0
count_wt_spike = 0
count_syn_mut = 0
count_wt = 0
count_variant = 0

# Parse fasta file sequence-by-sequence
for seq_record in SeqIO.parse(f"{args.string}.reversed.fasta", "fasta"):
    total_count += 1
    s = str(seq_record.seq)

    # Split sequence into codons (3 nucleotides each)
    seq_codons = [s[i:i+3] for i in range(0, len(s), 3)]

    # Count the number of codons that differ from WT
    diff_positions = [i for i in range(len(WT)) if seq_codons[i] != WT[i]]
    num_differences = len(diff_positions)
    if num_differences == 0:
        # Case 4: No differences - count all WT codons
        for i, codon in enumerate(WT):
            codons[f"codon{i + 1}"][codon] += 1
        count_wt_spike += 1

    elif num_differences == 1:
        # Case 1: Only one codon is different
        if diff_positions[0] == syn_codon_pos - 1:  # Adjust for 0-based index
            if seq_codons[syn_codon_pos - 1] == wt_codon:
                count_wt += 1
                continue  # Disregard this sequence
            else:
                codons[f"codon{syn_codon_pos}"][seq_codons[syn_codon_pos - 1]] += 1
                count_syn_mut +=1

    elif num_differences == 2:
        # Case 2: Two codons are different
        if syn_codon_pos - 1 in diff_positions:
            other_codon_pos = [pos for pos in diff_positions if pos != syn_codon_pos - 1][0]
            if seq_codons[syn_codon_pos - 1] == wt_codon:
                codons[f"codon{other_codon_pos + 1}"][seq_codons[other_codon_pos]] += 1
                count_variant += 1

print(f"total_count {total_count}")
print(f"count_wt_spike {count_wt_spike}")
print(f"count_syn_mut {count_syn_mut}")
print(f"count_wt {count_wt}")
print(f"count_variant {count_variant}")

    # Case 3: More than two codons are different - disregard
    # Implicitly handled by not processing further if num_differences > 2

# Save raw codon counts to a CSV file
df = pd.DataFrame(codons).fillna(0).astype(int)
df.rename(columns={0: 'Codon'}, inplace=True)
df.to_csv(f"{args.string}.raw_counts.csv", index=True)

# Replace "WT counts" to WT spike-in counts
WT_spike = df.at[syn_codon, f"codon{syn_codon_pos}"]
codon_count = 0
for x in WT:
    codon_count += 1
    df.at[x, f"codon{codon_count}"] = WT_spike
    if codon_count == 25:
        df.at["GGC", f"codon{codon_count}"] = 0


# Convert row names to AA and merge
Rownames = df.index.tolist()
aa_Rownames = translate_codons(Rownames)
df["Amino_Acid"] = aa_Rownames
df_grouped = df.groupby("Amino_Acid").sum().reset_index()
df_grouped.to_csv(f"{args.string}.WT_spike.csv", index=True)

# normalize counts and modifiy WT counts to spike-in counts
Length = len([1 for line in open(f"{args.string}.reversed.fasta") if line.startswith(">")])
table1_num = formater(df_grouped, Length)
table1_num.to_csv(f"{args.string}.WT_spike.RCPM.csv", index=True)

print("Processing complete!")
