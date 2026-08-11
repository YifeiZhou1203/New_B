from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import torch
import pickle
import re
from pathlib import Path
from collections import Counter

# ==========================================================
# Input files for the new NbT2T dataset
# Please Please Please Put this script in the same folder as these .gz files, otherwise, it wont work
# ==========================================================

AA_FASTA = "NbT2T.pep.fa.gz"
CDS_FASTA = "NbT2T.cds.fa.gz"

OUTPUT_PAIR_TXT = "clean_aa_codon_pairs_NbT2T.txt"
OUTPUT_TENSORS = "encoded_tensors_NbT2T.pt"
OUTPUT_VOCAB = "vocab_NbT2T.pkl"
OUTPUT_STATS = "preprocess_stats_NbT2T.txt"

# Optional filters; you can try more if you want
MIN_AA_LEN = 30
MAX_AA_LEN = 3000

STOP_CODONS = {"TAA", "TAG", "TGA"}
VALID_BASES = set("ATCG")
VALID_AAS = set("ARNDCQEGHILKMFPSTWYV")


def open_fasta(path):
    """Open normal FASTA or gzipped FASTA."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalize_id(text):
    """Make FASTA IDs easier to match between CDS and protein files"""
    text = str(text).strip()
    token = text.split()[0]

    # Try useful header attributes first.
    for pattern in [
        r"(?:transcript_id|transcript|Parent|ID|protein_id|gene)[:=]([^;\s]+)",
    ]:
        match = re.search(pattern, text)
        if match:
            token = match.group(1)
            break

    # Keep the most informative field if pipe-delimited.
    if "|" in token:
        parts = [p for p in token.split("|") if p]
        token = parts[-1]

    # Remove common CDS/protein suffixes.
    token = re.sub(r"\.(cds|pep|protein)$", "", token, flags=re.IGNORECASE)
    token = re.sub(r"(_cds|_pep|_protein)$", "", token, flags=re.IGNORECASE)

    # Some datasets use transcript.1 / transcript.p1 style.
    token = re.sub(r"\.p\d+$", "", token, flags=re.IGNORECASE)

    return token


def read_fasta_dict(path, seq_type):
    """Read FASTA into dictionary after ID normalization. If duplicate IDs exist, keep the longest sequence"""
    records = {}
    duplicate_ids = 0

    with open_fasta(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            rid = normalize_id(record.description)
            seq = str(record.seq).upper().replace("U", "T")
            seq = seq.replace(" ", "").replace("\n", "")

            if seq_type == "protein":
                seq = remove_terminal_stop_from_protein(seq)

            if rid in records:
                duplicate_ids += 1
                if len(seq) > len(records[rid]):
                    records[rid] = seq
            else:
                records[rid] = seq

    return records, duplicate_ids


def remove_terminal_stop_from_cds(dna_seq):
    """Remove terminal stop codon if present"""
    if len(dna_seq) >= 3 and dna_seq[-3:] in STOP_CODONS:
        return dna_seq[:-3]
    return dna_seq


#def remove_terminal_stop_from_protein(aa_seq):

def remove_terminal_stop_from_protein(aa_seq):
    """Remove terminal stop symbol from protein sequence if present.NbT2T uses '.' instead of '*'"""
    if aa_seq.endswith("*") or aa_seq.endswith("."):
        return aa_seq[:-1]
    return aa_seq


def split_codons(dna_seq):
    return [dna_seq[i:i + 3] for i in range(0, len(dna_seq), 3)]


def is_valid_dna(seq):
    return all(base in VALID_BASES for base in seq)


def is_valid_aa(seq):
    return all(aa in VALID_AAS for aa in seq)


def translate_cds_without_terminal_stop(dna_seq):
    """Translate after terminal stop codon has already been removed"""
    return str(Seq(dna_seq).translate())


def gc_content(seq):
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


# Check input files


for path in [AA_FASTA, CDS_FASTA]:
    if not Path(path).exists():
        raise FileNotFoundError(
            f"Cannot find {path}. Put this script in the same folder as the NbT2T files "
            f"or edit AA_FASTA/CDS_FASTA at the top of the script."
        )

# Read files


protein_dict, protein_duplicates = read_fasta_dict(AA_FASTA, "protein")
cds_dict, cds_duplicates = read_fasta_dict(CDS_FASTA, "cds")

print(f"Protein sequences loaded: {len(protein_dict)}")
print(f"CDS sequences loaded:     {len(cds_dict)}")
print(f"Duplicate protein IDs handled: {protein_duplicates}")
print(f"Duplicate CDS IDs handled:     {cds_duplicates}")

# Pair and clean sequences


clean_pairs = []

stats = Counter()
matched_ids = sorted(set(protein_dict) & set(cds_dict))
stats["matched_aa_cds_ids"] = len(matched_ids)
stats["protein_only_ids"] = len(set(protein_dict) - set(cds_dict))
stats["cds_only_ids"] = len(set(cds_dict) - set(protein_dict))

for gene_id in matched_ids:
    aa_seq = protein_dict[gene_id]
    dna_seq = remove_terminal_stop_from_cds(cds_dict[gene_id])

    if len(dna_seq) == 0 or len(aa_seq) == 0:
        stats["removed_empty_sequence"] += 1
        continue

    if len(dna_seq) % 3 != 0:
        stats["removed_length_not_multiple_of_3"] += 1
        continue

    if not is_valid_dna(dna_seq):
        stats["removed_invalid_dna"] += 1
        continue

    if not is_valid_aa(aa_seq):
        stats["removed_invalid_aa"] += 1
        continue

    if MIN_AA_LEN is not None and len(aa_seq) < MIN_AA_LEN:
        stats["removed_too_short"] += 1
        continue

    if MAX_AA_LEN is not None and len(aa_seq) > MAX_AA_LEN:
        stats["removed_too_long"] += 1
        continue

    translated_aa = translate_cds_without_terminal_stop(dna_seq)

    if "*" in translated_aa:
        stats["removed_internal_stop"] += 1
        continue

    if translated_aa != aa_seq:
        stats["removed_translation_mismatch"] += 1
        continue

    codons = split_codons(dna_seq)

    if len(codons) != len(aa_seq):
        stats["removed_aa_codon_length_mismatch"] += 1
        continue

    clean_pairs.append((gene_id, list(aa_seq), codons))

stats["clean_valid_pairs"] = len(clean_pairs)

print("\n=== Cleaning Summary ===")
for key in [
    "matched_aa_cds_ids",
    "clean_valid_pairs",
    "protein_only_ids",
    "cds_only_ids",
    "removed_empty_sequence",
    "removed_length_not_multiple_of_3",
    "removed_invalid_dna",
    "removed_invalid_aa",
    "removed_too_short",
    "removed_too_long",
    "removed_internal_stop",
    "removed_translation_mismatch",
    "removed_aa_codon_length_mismatch",
]:
    print(f"{key:36s}: {stats[key]}")

if not clean_pairs:
    raise ValueError(
        "No clean AA-CDS pairs were found. The most likely cause is ID mismatch between "
        "NbT2T.pep.fa.gz and NbT2T.cds.fa.gz. Print a few FASTA headers and adjust normalize_id()."
    )

# Save clean text pairs

with open(OUTPUT_PAIR_TXT, "w") as f:
    for gene_id, aa_list, codon_list in clean_pairs:
        aa_seq = "".join(aa_list)
        codon_seq = " ".join(codon_list)
        f.write(f"{gene_id}\t{aa_seq}\t{codon_seq}\n")

print(f"\nSaved clean pairs to: {OUTPUT_PAIR_TXT}")

# Build vocabularies


aa_vocab = {
    "<PAD>": 0, "<SOS>": 1, "<EOS>": 2,
    "A": 3, "R": 4, "N": 5, "D": 6, "C": 7,
    "Q": 8, "E": 9, "G": 10, "H": 11, "I": 12,
    "L": 13, "K": 14, "M": 15, "F": 16, "P": 17,
    "S": 18, "T": 19, "W": 20, "Y": 21, "V": 22,
}

bases = ["T", "C", "A", "G"]
all_codons = [a + b + c for a in bases for b in bases for c in bases]

codon_vocab = {"<PAD>": 0, "<SOS>": 1, "<EOS>": 2}
for i, codon in enumerate(all_codons, start=3):
    codon_vocab[codon] = i

# If your model should never generate stop codons, keep this mask for decoding/training.
sense_codons = [c for c in all_codons if c not in STOP_CODONS]
stop_codon_ids = [codon_vocab[c] for c in STOP_CODONS]
sense_codon_ids = [codon_vocab[c] for c in sense_codons]


# Encoding functions: not necessary, can do it further in the test.py later

def encode_aa(seq):
    return [aa_vocab["<SOS>"]] + [aa_vocab[aa] for aa in seq] + [aa_vocab["<EOS>"]]


def encode_codons(codons):
    return [codon_vocab["<SOS>"]] + [codon_vocab[codon] for codon in codons] + [codon_vocab["<EOS>"]]


def pad_sequences(seqs, pad_value):
    max_len = max(len(seq) for seq in seqs)
    return torch.tensor([seq + [pad_value] * (max_len - len(seq)) for seq in seqs], dtype=torch.long)


# Encode clean pairs

src_sequences = []
tgt_sequences = []
gene_ids = []
lengths_aa = []
gc_values = []

for gene_id, aa_list, codon_list in clean_pairs:
    src_sequences.append(encode_aa(aa_list))
    tgt_sequences.append(encode_codons(codon_list))
    gene_ids.append(gene_id)
    lengths_aa.append(len(aa_list))
    gc_values.append(gc_content("".join(codon_list)))

src_tensor = pad_sequences(src_sequences, aa_vocab["<PAD>"])
tgt_tensor = pad_sequences(tgt_sequences, codon_vocab["<PAD>"])

print("\n=== Tensor Shapes ===")
print(f"Source tensor shape: {tuple(src_tensor.shape)}")
print(f"Target tensor shape: {tuple(tgt_tensor.shape)}")
print(f"Average AA length:   {sum(lengths_aa) / len(lengths_aa):.2f}")
print(f"Max AA length:       {max(lengths_aa)}")
print(f"Average CDS GC:      {sum(gc_values) / len(gc_values):.4f}")


torch.save(
    {
        "gene_ids": gene_ids,
        "src_tensor": src_tensor,
        "tgt_tensor": tgt_tensor,
        "lengths_aa": lengths_aa,
        "gc_values": gc_values,
    },
    OUTPUT_TENSORS,
)

with open(OUTPUT_VOCAB, "wb") as f:
    pickle.dump(
        {
            "aa_vocab": aa_vocab,
            "codon_vocab": codon_vocab,
            "sense_codon_ids": sense_codon_ids,
            "stop_codon_ids": stop_codon_ids,
        },
        f,
    )

with open(OUTPUT_STATS, "w") as f:
    f.write("=== NbT2T preprocessing stats ===\n")
    for key, value in stats.items():
        f.write(f"{key}: {value}\n")
    f.write(f"average_aa_length: {sum(lengths_aa) / len(lengths_aa):.2f}\n")
    f.write(f"max_aa_length: {max(lengths_aa)}\n")
    f.write(f"average_cds_gc: {sum(gc_values) / len(gc_values):.4f}\n")

print(f"\nSaved tensors to: {OUTPUT_TENSORS}")
print(f"Saved vocab to:   {OUTPUT_VOCAB}")
print(f"Saved stats to:   {OUTPUT_STATS}")

# Preview first example

first_gene, first_aa, first_codons = clean_pairs[0]
print("\n=== First Clean Example ===")
print("Gene ID:", first_gene)
print("AA sequence first 20:", "".join(first_aa[:20]))
print("Codons first 20:", first_codons[:20])
