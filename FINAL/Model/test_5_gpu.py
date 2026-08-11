from Bio import SeqIO
import torch
from torch.utils.data import Dataset, DataLoader, random_split
from torch.nn.utils.rnn import pad_sequence
import math
import torch.nn as nn
import random
import warnings
from collections import defaultdict, Counter
import numpy as np
import matplotlib.pyplot as plt
import pickle
import subprocess
import shutil
import re


warnings.filterwarnings(
    "ignore",
    message="The PyTorch API of nested tensors is in prototype stage"
)


CLEAN_PAIR_FILE = "clean_aa_codon_pairs_NbT2T.txt"

MAX_AA_LENGTH = 480
MAX_SAMPLES = 20000
BATCH_SIZE = 16
NUM_EPOCHS = 25
LEARNING_RATE = 0.0005
RANDOM_SEED = 1234

TARGET_GC = 0.43
LAMBDA_GC = 3


TRNA_SCAN_FILE = "nbenthamiana_NbT2T_trna.out"
USE_TRNASCAN_TAI = True
LAMBDA_TAI = 0.1
TARGET_TAI = 0.373
TARGET_MFE_PROXY = 0.48
LAMBDA_MFE = 0.1

GC_DECODING_BONUS = 1.5
TAI_DECODING_BONUS = 0
USE_POST_MFE_SELECTION = False
POST_TOP_K = 4
POST_NUM_CANDIDATES = 12
POST_MFE_MAX_RNAFOLD_LEN = 3000
TARGET_RNAFOLD_MFE_PER_NT = -0.30
POST_GC_WEIGHT = 1
POST_TAI_WEIGHT = 0
POST_MFE_WEIGHT = 0.3


RUN_EXACT_RNAFOLD_ON_TEST_SET = False


FORBIDDEN_RESTRICTION_SITES = {
    "BamHI": "GGATCC",
    "BglII": "AGATCT",
    "KasI": "GGCGCC",
    "SacI": "GAGCTC",
}
POLYA_SIGNALS = {
    "canonical_polyA": "AATAAA",
    "variant_polyA": "ATTAAA",
}
MOTIF_PENALTY = 10.0
POLYA_PENALTY = 10.0

GENERATE_RNAFOLD_CANDIDATES = True
NUM_SEQUENCES_FOR_CANDIDATES = 10
CANDIDATES_PER_SEQUENCE = 12
CANDIDATE_TOP_K = 4
RNAFOLD_CANDIDATE_FASTA = "final_candidates_for_rnafold_NbT2T.fasta"
RNAFOLD_CANDIDATE_TSV = "final_candidates_for_rnafold_summary_NbT2T.tsv"

random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)

print("Dataset file:", CLEAN_PAIR_FILE)
print("Recommended for adalimumab HC: MAX_AA_LENGTH=480 covers the 473-aa heavy chain while staying faster than 500/512.")
print("Windows CUDA version: put this script in the same folder as clean_aa_codon_pairs_NbT2T.txt and nbenthamiana_NbT2T_trna.out.")


genetic_code = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

aa_to_synonymous_codons = defaultdict(list)
for codon, aa in genetic_code.items():
    if aa != "*":
        aa_to_synonymous_codons[aa].append(codon)


def build_valid_codon_mask(aa_vocab, codon_vocab, aa_to_synonymous_codons, device):
    valid_mask = torch.zeros(len(aa_vocab), len(codon_vocab), dtype=torch.bool, device=device)

    for aa, codons in aa_to_synonymous_codons.items():
        aa_idx = aa_vocab[aa]
        for codon in codons:
            codon_idx = codon_vocab[codon]
            valid_mask[aa_idx, codon_idx] = True

    return valid_mask


def valid_aa(seq):
    aa_set = set("ARNDCQEGHILKMFPSTWYV")
    return all(s in aa_set for s in seq)


def valid_codon_list(codons):
    return all(c in genetic_code and genetic_code[c] != "*" for c in codons)


def load_clean_pairs(path):
    loaded_pairs = []
    skipped_bad_format = 0
    skipped_bad_aa = 0
    skipped_bad_codons = 0
    skipped_len_mismatch = 0

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")

            if len(parts) == 3:
                gene_id, aa_seq, codon_text = parts
            elif len(parts) == 2:
                gene_id = None
                aa_seq, codon_text = parts
            else:
                skipped_bad_format += 1
                continue

            aa_seq = aa_seq.strip().upper()
            codons = [c.strip().upper() for c in codon_text.split() if c.strip()]

            if not valid_aa(aa_seq):
                skipped_bad_aa += 1
                continue

            if not valid_codon_list(codons):
                skipped_bad_codons += 1
                continue

            if len(aa_seq) != len(codons):
                skipped_len_mismatch += 1
                continue

            loaded_pairs.append((list(aa_seq), codons))

    print("Loaded clean pairs:", len(loaded_pairs))
    print("Skipped bad format:", skipped_bad_format)
    print("Skipped bad AA:", skipped_bad_aa)
    print("Skipped bad codons:", skipped_bad_codons)
    print("Skipped length mismatch:", skipped_len_mismatch)

    return loaded_pairs


clean_pairs = load_clean_pairs(CLEAN_PAIR_FILE)


short_pairs = []
for aa_seq, codons in clean_pairs:
    if len(aa_seq) <= MAX_AA_LENGTH:
        short_pairs.append((aa_seq, codons))

print("After length filter:", len(short_pairs))

if MAX_SAMPLES is not None and len(short_pairs) > MAX_SAMPLES:
    clean_pairs = random.sample(short_pairs, MAX_SAMPLES)
else:
    clean_pairs = short_pairs

print("Final dataset size:", len(clean_pairs))

if len(clean_pairs) == 0:
    raise ValueError("No usable sequence pairs left. Check CLEAN_PAIR_FILE and MAX_AA_LENGTH.")

with open("clean_pairs_used_for_training_NbT2T.txt", "w") as w:
    for aa_seq, codons in clean_pairs:
        w.write(f"{''.join(aa_seq)}\t{' '.join(codons)}\n")


bases = ["T", "C", "A", "G"]
all_codons = []

for a in bases:
    for b in bases:
        for c in bases:
            all_codons.append(a + b + c)

codon_vocab = {
    "<PAD>": 0,
    "<SOS>": 1,
    "<EOS>": 2
}

aa_vocab = {
    "<PAD>": 0,
    "<SOS>": 1,
    "<EOS>": 2,
    "A": 3, "R": 4, "N": 5, "D": 6, "C": 7,
    "Q": 8, "E": 9, "G": 10, "H": 11, "I": 12,
    "L": 13, "K": 14, "M": 15, "F": 16, "P": 17,
    "S": 18, "T": 19, "W": 20, "Y": 21, "V": 22
}

index = 3
for codon in all_codons:
    codon_vocab[codon] = index
    index += 1

idx_to_codon = {v: k for k, v in codon_vocab.items()}
idx_to_aa = {v: k for k, v in aa_vocab.items()}

def encode_aa(seq):
    return [aa_vocab["<SOS>"]] + [aa_vocab[s] for s in seq] + [aa_vocab["<EOS>"]]


def encode_codon(seq):
    return [codon_vocab["<SOS>"]] + [codon_vocab[s] for s in seq] + [codon_vocab["<EOS>"]]


src_sequences = [encode_aa(aa) for aa, codons in clean_pairs]
tgt_sequences = [encode_codon(codons) for aa, codons in clean_pairs]


class CodonDataset(Dataset):
    def __init__(self, src_sequences, tgt_sequences, gc_target):
        self.src_sequences = src_sequences
        self.tgt_sequences = tgt_sequences
        self.gc_target = gc_target

    def __len__(self):
        return len(self.src_sequences)

    def __getitem__(self, idx):
        src = torch.tensor(self.src_sequences[idx], dtype=torch.long)
        tgt = torch.tensor(self.tgt_sequences[idx], dtype=torch.long)
        gc = torch.tensor([self.gc_target], dtype=torch.float32)
        return src, tgt, gc

def collate_fn(batch):
    src_batch = [item[0] for item in batch]
    tgt_batch = [item[1] for item in batch]
    gc_batch = torch.stack([item[2] for item in batch])

    src_padded = pad_sequence(
        src_batch,
        batch_first=True,
        padding_value=aa_vocab["<PAD>"]
    )
    tgt_padded = pad_sequence(
        tgt_batch,
        batch_first=True,
        padding_value=codon_vocab["<PAD>"]
    )
    return src_padded, tgt_padded, gc_batch


full_dataset = CodonDataset(src_sequences, tgt_sequences, TARGET_GC)

total_size = len(full_dataset)
train_size = int(0.8 * total_size)
val_size = int(0.1 * total_size)
test_size = total_size - train_size - val_size

split_generator = torch.Generator().manual_seed(RANDOM_SEED)

train_dataset, val_dataset, test_dataset = random_split(
    full_dataset,
    [train_size, val_size, test_size],
    generator=split_generator
)

train_loader = DataLoader(
    train_dataset,
    batch_size=BATCH_SIZE,
    shuffle=True,
    collate_fn=collate_fn
)
val_loader = DataLoader(
    val_dataset,
    batch_size=BATCH_SIZE,
    shuffle=False,
    collate_fn=collate_fn
)
test_loader = DataLoader(
    test_dataset,
    batch_size=BATCH_SIZE,
    shuffle=False,
    collate_fn=collate_fn
)


def build_empirical_codon_distribution(train_dataset, aa_vocab, codon_vocab, idx_to_aa, idx_to_codon):
    aa_to_codon_counts = defaultdict(Counter)

    for i in range(len(train_dataset)):
        src_tensor, tgt_tensor, _ = train_dataset[i]

        aa_seq = []
        for idx in src_tensor.tolist():
            aa = idx_to_aa[idx]
            if aa not in ["<SOS>", "<EOS>", "<PAD>"]:
                aa_seq.append(aa)

        codon_seq = []
        for idx in tgt_tensor.tolist():
            codon = idx_to_codon[idx]
            if codon not in ["<SOS>", "<EOS>", "<PAD>"]:
                codon_seq.append(codon)

        for aa, codon in zip(aa_seq, codon_seq):
            aa_to_codon_counts[aa][codon] += 1

    q_table = torch.zeros(len(aa_vocab), len(codon_vocab))

    for aa, counter in aa_to_codon_counts.items():
        aa_idx = aa_vocab[aa]
        total = sum(counter.values())

        for codon, count in counter.items():
            codon_idx = codon_vocab[codon]
            q_table[aa_idx, codon_idx] = count / total

    return q_table


def compute_entropy_from_distribution(q_table, eps=1e-8):
    q_safe = q_table.clamp(min=eps)
    entropy = -(q_safe * torch.log(q_safe)).sum(dim=-1)
    return entropy


class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=6000):
        super().__init__()

        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)

        div_term = torch.exp(
            torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model)
        )

        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)

        pe = pe.unsqueeze(0)
        self.register_buffer("pe", pe)

    def forward(self, x):
        seq_len = x.size(1)
        return x + self.pe[:, :seq_len, :]


class EncoderOnlyCodonTransformer(nn.Module):
    def __init__(
        self,
        src_vocab_size,
        tgt_vocab_size,
        d_model=128,
        nhead=4,
        num_encoder_layers=4,
        dim_feedforward=128,
        dropout=0.1
    ):
        super().__init__()

        self.d_model = d_model
        self.src_embedding = nn.Embedding(src_vocab_size, d_model)
        self.gc_projection = nn.Linear(1, d_model)
        self.positional_encoding = PositionalEncoding(d_model)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first=True
        )

        self.encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers=num_encoder_layers,
            enable_nested_tensor=False
        )

        self.fc_out = nn.Linear(d_model, tgt_vocab_size)
    def forward(self, src, gc_target, src_key_padding_mask=None):
        x = self.src_embedding(src) * math.sqrt(self.d_model)
        gc_embed = self.gc_projection(gc_target)
        gc_embed = gc_embed.unsqueeze(1)
        gc_embed = gc_embed.expand(-1, src.size(1), -1)

        x = x + gc_embed
        x = self.positional_encoding(x)
        x = self.encoder(x, src_key_padding_mask=src_key_padding_mask)
        output = self.fc_out(x)
        return output


def build_aligned_targets(src_batch, tgt_batch):
    batch_size, src_len = src_batch.shape

    aligned = torch.full(
        (batch_size, src_len),
        codon_vocab["<PAD>"],
        dtype=torch.long,
        device=src_batch.device
    )

    for i in range(batch_size):
        tgt_seq = tgt_batch[i]
        tgt_real_len = int((tgt_seq != codon_vocab["<PAD>"]).sum().item())
        codon_tokens = tgt_seq[1:tgt_real_len - 1]
        aligned[i, 1:1 + len(codon_tokens)] = codon_tokens

    return aligned


def build_gc_weight_vector(codon_vocab, device):
    gc_weights = torch.zeros(len(codon_vocab), device=device)

    for codon, idx in codon_vocab.items():
        if codon in ["<PAD>", "<SOS>", "<EOS>"]:
            gc_weights[idx] = 0.0
        else:
            gc_weights[idx] = sum(1 for b in codon if b in ["G", "C"]) / 3.0

    return gc_weights

def expected_gc_from_logits(logits, gc_weights):
    probs = torch.softmax(logits, dim=-1)
    expected_gc = (probs * gc_weights.view(1, 1, -1)).sum(dim=-1)
    return expected_gc


def compute_gc_loss(logits, aligned_target, gc_batch, gc_weights):
    expected_gc = expected_gc_from_logits(logits, gc_weights)

    valid_mask = (aligned_target != codon_vocab["<PAD>"]).float()

    gc_sum = (expected_gc * valid_mask).sum(dim=1)
    valid_count = valid_mask.sum(dim=1).clamp(min=1.0)
    avg_gc_per_seq = gc_sum / valid_count

    target_gc = gc_batch.squeeze(1)

    gc_loss = ((avg_gc_per_seq - target_gc) ** 2).mean()
    return gc_loss, avg_gc_per_seq.mean().item()


def reverse_complement(seq):
    comp = str.maketrans("ATCG", "TAGC")
    return seq.upper().translate(comp)[::-1]


aa1_to_aa3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"
}
aa3_to_aa1 = {v: k for k, v in aa1_to_aa3.items()}
aa3_to_aa1["iMet"] = "M"


def parse_trnascan_to_codon_counts(trna_file, genetic_code):
    codon_counts = Counter()
    anticodon_counts = defaultdict(Counter)
    parsed_rows = 0
    used_rows = 0
    skipped_rows = 0

    standard_types = set(aa3_to_aa1.keys())

    try:
        f = open(trna_file, "r")
    except FileNotFoundError:
        print(f"WARNING: {trna_file} not found. Falling back to training codon-usage tAI proxy.")
        return codon_counts, anticodon_counts, parsed_rows, used_rows, skipped_rows

    with f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Sequence") or line.startswith("Name") or line.startswith("-"):
                continue

            parts = line.split()
            if len(parts) < 6:
                continue

            parsed_rows += 1

            trna_type = parts[4]
            anticodon = parts[5].upper()
            lower_line = line.lower()

            if "pseudo" in lower_line:
                skipped_rows += 1
                continue
            if trna_type not in standard_types:
                skipped_rows += 1
                continue
            if len(anticodon) != 3 or any(b not in "ATCG" for b in anticodon):
                skipped_rows += 1
                continue

            aa1 = aa3_to_aa1[trna_type]
            codon = reverse_complement(anticodon)

            if genetic_code.get(codon) != aa1:
                skipped_rows += 1
                continue

            anticodon_counts[aa1][anticodon] += 1
            codon_counts[codon] += 1
            used_rows += 1

    return codon_counts, anticodon_counts, parsed_rows, used_rows, skipped_rows


def build_tai_weights_from_trnascan(trna_file, aa_to_synonymous_codons, codon_vocab, device):
    codon_counts, anticodon_counts, parsed_rows, used_rows, skipped_rows = parse_trnascan_to_codon_counts(
        trna_file,
        genetic_code
    )

    if used_rows == 0:
        return None, codon_counts, anticodon_counts

    pseudocount = 0.5
    weights = torch.zeros(len(codon_vocab), device=device)

    for aa, codons in aa_to_synonymous_codons.items():
        raw_counts = [codon_counts[c] for c in codons]
        max_count = max(raw_counts)

        if max_count == 0:
            for codon in codons:
                weights[codon_vocab[codon]] = 1.0
            continue

        for codon in codons:
            count = codon_counts[codon]
            weights[codon_vocab[codon]] = (count + pseudocount) / (max_count + pseudocount)

    print("\n=== tRNAscan-SE tAI-like weights ===")
    print(f"tRNAscan file: {trna_file}")
    print(f"Parsed tRNA rows: {parsed_rows}")
    print(f"Used standard tRNA rows: {used_rows}")
    print(f"Skipped rows: {skipped_rows}")
    print("Interpretation: tRNA gene-copy-based tAI-like proxy, not experimental tRNA abundance.")

    with open("trnascan_codon_gene_counts_NbT2T.csv", "w") as out:
        out.write("amino_acid,codon,trna_gene_count,tai_gene_copy_weight\n")
        for aa in sorted(aa_to_synonymous_codons.keys()):
            for codon in sorted(aa_to_synonymous_codons[aa]):
                out.write(
                    f"{aa},{codon},{codon_counts[codon]},"
                    f"{float(weights[codon_vocab[codon]].detach().cpu().item()):.6f}\n"
                )

    print("Saved: trnascan_codon_gene_counts_NbT2T.csv")
    return weights, codon_counts, anticodon_counts


def build_codon_adaptation_weight_vector(q_codon_given_aa, aa_to_synonymous_codons, aa_vocab, codon_vocab, device):
    weights = torch.zeros(len(codon_vocab), device=device)

    for aa, codons in aa_to_synonymous_codons.items():
        aa_idx = aa_vocab[aa]
        codon_indices = [codon_vocab[c] for c in codons]
        freqs = q_codon_given_aa[aa_idx, codon_indices]
        max_freq = freqs.max().clamp(min=1e-8)

        for codon in codons:
            codon_idx = codon_vocab[codon]
            weights[codon_idx] = q_codon_given_aa[aa_idx, codon_idx] / max_freq

    print("\nWARNING: Using training codon-usage tAI proxy, not tRNAscan-SE tRNA gene-copy weights.")
    return weights

def compute_tai_loss(logits, aligned_target, tai_weights, target_tai=0.85, eps=1e-8):
    probs = torch.softmax(logits, dim=-1)
    valid_mask = (aligned_target != codon_vocab["<PAD>"]).float()


    expected_tai_per_pos = (probs * tai_weights.view(1, 1, -1)).sum(dim=-1)

    tai_sum = (expected_tai_per_pos * valid_mask).sum(dim=1)
    valid_count = valid_mask.sum(dim=1).clamp(min=1.0)
    avg_tai_per_seq = tai_sum / valid_count

    target = torch.full_like(avg_tai_per_seq, target_tai)
    tai_loss = ((avg_tai_per_seq - target) ** 2).mean()

    return tai_loss, avg_tai_per_seq.mean().item()


def build_codon_pairing_weight_vector(codon_vocab, device):
    weights = torch.zeros(len(codon_vocab), device=device)

    for codon, idx in codon_vocab.items():
        if codon in ["<PAD>", "<SOS>", "<EOS>"]:
            weights[idx] = 0.0
        else:
            rna = codon.replace("T", "U")
            gc = sum(1 for b in rna if b in ["G", "C"]) / 3.0
            au = sum(1 for b in rna if b in ["A", "U"]) / 3.0
            weights[idx] = 0.7 * gc + 0.3 * au

    return weights


def compute_mfe_proxy_loss(logits, aligned_target, mfe_proxy_weights, target_mfe_proxy=0.18):
    probs = torch.softmax(logits, dim=-1)
    valid_mask = (aligned_target != codon_vocab["<PAD>"]).float()

    expected_pairing_per_pos = (probs * mfe_proxy_weights.view(1, 1, -1)).sum(dim=-1)

    pairing_sum = (expected_pairing_per_pos * valid_mask).sum(dim=1)
    valid_count = valid_mask.sum(dim=1).clamp(min=1.0)
    avg_pairing_per_seq = pairing_sum / valid_count

    target = torch.full_like(avg_pairing_per_seq, target_mfe_proxy)
    mfe_proxy_loss = ((avg_pairing_per_seq - target) ** 2).mean()

    return mfe_proxy_loss, avg_pairing_per_seq.mean().item()


def tai_of_codon_sequence(codon_list, tai_weight_dict):
    weights = []
    for codon in codon_list:
        w = tai_weight_dict.get(codon, 0.0)
        if w > 0:
            weights.append(w)
    if len(weights) == 0:
        return 0.0
    return float(np.exp(np.mean(np.log(weights))))


def summarize_dataset_tai(clean_pairs, tai_weight_dict, output_path="dataset_tai_summary_NbT2T.txt"):
    dataset_tais = [tai_of_codon_sequence(codons, tai_weight_dict) for _, codons in clean_pairs]
    dataset_tais = np.array(dataset_tais, dtype=float)

    summary = {
        "mean": float(np.mean(dataset_tais)),
        "median": float(np.median(dataset_tais)),
        "std": float(np.std(dataset_tais)),
        "min": float(np.min(dataset_tais)),
        "max": float(np.max(dataset_tais)),
        "n": int(len(dataset_tais)),
    }

    print("\n=== Dataset native tAI-like summary ===")
    print(f"Dataset average tAI-like: {summary['mean']:.4f}")
    print(f"Dataset median tAI-like:  {summary['median']:.4f}")
    print(f"Dataset tAI-like std:     {summary['std']:.4f}")
    print(f"Dataset tAI-like range:   {summary['min']:.4f} - {summary['max']:.4f}")
    print("Interpretation: this is a tRNAscan-SE gene-copy-based tAI-like proxy, not true experimental tAI.")

    with open(output_path, "w") as out:
        for key, value in summary.items():
            out.write(f"{key}\t{value}\n")

    print(f"Saved: {output_path}")
    return summary


def mfe_proxy_of_codon_sequence(codon_list):
    if len(codon_list) == 0:
        return 0.0
    scores = []
    for codon in codon_list:
        rna = codon.replace("T", "U")
        gc = sum(1 for b in rna if b in ["G", "C"]) / 3.0
        au = sum(1 for b in rna if b in ["A", "U"]) / 3.0
        scores.append(0.7 * gc + 0.3 * au)
    return float(np.mean(scores))


def run_rnafold_mfe(rna_sequence, max_len=5000):
    if shutil.which("RNAfold") is None or len(rna_sequence) > max_len:
        return None

    result = subprocess.run(
        ["RNAfold", "--noPS"],
        input=rna_sequence + "\n",
        text=True,
        capture_output=True,
        check=False
    )

    if result.returncode != 0:
        return None

    match = re.search(r"\(([-+]?[0-9]*\.?[0-9]+)\)", result.stdout)
    if match is None:
        return None

    return float(match.group(1))

def count_motif_hits(dna_sequence, motif_dict):
    hits = {}
    total = 0
    seq = dna_sequence.upper()
    for name, motif in motif_dict.items():
        motif = motif.upper()
        count = seq.count(motif)
        hits[name] = count
        total += count
    return hits, total


def motif_penalty_for_codons(codon_list):
    dna = "".join(codon_list)
    restriction_hits, restriction_total = count_motif_hits(dna, FORBIDDEN_RESTRICTION_SITES)
    polya_hits, polya_total = count_motif_hits(dna, POLYA_SIGNALS)
    penalty = MOTIF_PENALTY * restriction_total + POLYA_PENALTY * polya_total
    return penalty, restriction_hits, polya_hits, restriction_total, polya_total


def score_sequence_for_post_selection(codon_list, tai_weight_dict, target_gc, target_mfe_proxy):
    if len(codon_list) == 0:
        return float("inf"), {"gc": 0.0, "tai": 0.0, "mfe_proxy": 0.0, "rnafold_mfe": None}

    seq = "".join(codon_list)
    gc = (sum(1 for b in seq if b in ["G", "C"]) / len(seq)) if len(seq) > 0 else 0.0
    tai = tai_of_codon_sequence(codon_list, tai_weight_dict)
    mfe_proxy = mfe_proxy_of_codon_sequence(codon_list)

    rna_sequence = "".join(codon_list).replace("T", "U")
    rnafold_mfe = run_rnafold_mfe(rna_sequence, max_len=POST_MFE_MAX_RNAFOLD_LEN)

    gc_penalty = abs(gc - target_gc)
    tai_bonus = tai

    if rnafold_mfe is not None and len(rna_sequence) > 0:
        mfe_per_nt = rnafold_mfe / len(rna_sequence)
        mfe_penalty = abs(mfe_per_nt - TARGET_RNAFOLD_MFE_PER_NT)
    else:
        mfe_penalty = abs(mfe_proxy - target_mfe_proxy)

    motif_penalty, restriction_hits, polya_hits, restriction_total, polya_total = motif_penalty_for_codons(codon_list)

    score = (
        POST_GC_WEIGHT * gc_penalty
        - POST_TAI_WEIGHT * tai_bonus
        + POST_MFE_WEIGHT * mfe_penalty
        + motif_penalty
    )

    details = {
        "gc": gc,
        "tai": tai,
        "mfe_proxy": mfe_proxy,
        "rnafold_mfe": rnafold_mfe,
        "restriction_site_total": restriction_total,
        "polyA_signal_total": polya_total,
        "restriction_hits": restriction_hits,
        "polyA_hits": polya_hits
    }
    return score, details


def post_mfe_select_codons_for_sequence(
    logits_seq,
    src_seq,
    valid_positions,
    gc_weights,
    tai_weights,
    tai_weight_dict,
    top_k=3,
    num_candidates=8
):
    adjusted_logits = (
        logits_seq
        + GC_DECODING_BONUS * gc_weights.view(1, -1)
        + TAI_DECODING_BONUS * tai_weights.view(1, -1)
    )

    position_choices = []

    for pos in valid_positions.tolist():
        scores = adjusted_logits[pos]
        finite_idx = torch.where(scores > -1e8)[0]

        if finite_idx.numel() == 0:
            position_choices.append([])
            continue

        k = min(top_k, finite_idx.numel())
        top_values, top_indices_local = torch.topk(scores[finite_idx], k=k)
        top_indices = finite_idx[top_indices_local]
        position_choices.append(top_indices.detach().cpu().tolist())

    candidates = []

    greedy = []
    for choices in position_choices:
        if len(choices) > 0:
            greedy.append(idx_to_codon[choices[0]])
    candidates.append(greedy)

    for _ in range(max(0, num_candidates - 1)):
        cand = []
        for choices in position_choices:
            if len(choices) == 0:
                continue
            if len(choices) == 1:
                chosen_idx = choices[0]
            else:
                chosen_idx = random.choice(choices)
            cand.append(idx_to_codon[chosen_idx])
        candidates.append(cand)

    best_score = float("inf")
    best_candidate = candidates[0]
    best_details = None

    for cand in candidates:
        score, details = score_sequence_for_post_selection(
            cand,
            tai_weight_dict,
            target_gc=TARGET_GC,
            target_mfe_proxy=TARGET_MFE_PROXY
        )
        if score < best_score:
            best_score = score
            best_candidate = cand
            best_details = details

    return best_candidate, best_score, best_details


def compute_synonymous_entropy_loss(
    logits,
    src_batch,
    aligned_target,
    q_codon_given_aa,
    target_entropy_by_aa,
    alpha=0.7,
    lambda_entropy=0.2,
    eps=1e-8
):

    valid_mask = aligned_target != codon_vocab["<PAD>"]

    log_probs = torch.log_softmax(logits, dim=-1)
    probs = torch.softmax(logits, dim=-1)

    q_soft = q_codon_given_aa[src_batch]

    one_hot = torch.zeros_like(logits)
    safe_target = aligned_target.clone()
    safe_target[safe_target == codon_vocab["<PAD>"]] = 0
    one_hot.scatter_(-1, safe_target.unsqueeze(-1), 1.0)

    mixed_target = (1 - alpha) * one_hot + alpha * q_soft

    soft_ce_per_pos = -(mixed_target * log_probs).sum(dim=-1)
    soft_ce_loss = soft_ce_per_pos[valid_mask].mean()

    pred_entropy = -(probs.clamp(min=eps) * torch.log(probs.clamp(min=eps))).sum(dim=-1)

    target_entropy = target_entropy_by_aa[src_batch]

    entropy_loss = ((pred_entropy - target_entropy) ** 2)[valid_mask].mean()

    total_loss = soft_ce_loss + lambda_entropy * entropy_loss

    return total_loss, soft_ce_loss, entropy_loss


if torch.cuda.is_available():
    device = torch.device("cuda")
    print("Using device:", device)
    print("CUDA GPU:", torch.cuda.get_device_name(0))
    print("CUDA version used by PyTorch:", torch.version.cuda)
    torch.backends.cudnn.benchmark = True
elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
    device = torch.device("mps")
    print("Using device:", device)
    print("MPS available. This is slower than NVIDIA CUDA for this model.")
else:
    device = torch.device("cpu")
    print("Using device:", device)
    print("WARNING: No CUDA/MPS GPU detected. Training will be slow.")

print("MAX_AA_LENGTH:", MAX_AA_LENGTH)
print("BATCH_SIZE:", BATCH_SIZE)
print("If CUDA runs out of memory, set BATCH_SIZE = 8. If stable, try 24 or 32.")

valid_codon_mask = build_valid_codon_mask(
    aa_vocab,
    codon_vocab,
    aa_to_synonymous_codons,
    device
)

print("GC target conditioning value:", TARGET_GC)
print("GC penalty lambda:", LAMBDA_GC)


model = EncoderOnlyCodonTransformer(
    src_vocab_size=len(aa_vocab),
    tgt_vocab_size=len(codon_vocab)
).to(device)

print(model)


optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)

gc_weights = build_gc_weight_vector(codon_vocab, device)

q_codon_given_aa = build_empirical_codon_distribution(
    train_dataset,
    aa_vocab,
    codon_vocab,
    idx_to_aa,
    idx_to_codon
).to(device)

target_entropy_by_aa = compute_entropy_from_distribution(
    q_codon_given_aa
).to(device)

if USE_TRNASCAN_TAI:
    tai_weights, trnascan_codon_counts, trnascan_anticodon_counts = build_tai_weights_from_trnascan(
        TRNA_SCAN_FILE,
        aa_to_synonymous_codons,
        codon_vocab,
        device
    )
else:
    tai_weights = None

if tai_weights is None:
    tai_weights = build_codon_adaptation_weight_vector(
        q_codon_given_aa,
        aa_to_synonymous_codons,
        aa_vocab,
        codon_vocab,
        device
    )

mfe_proxy_weights = build_codon_pairing_weight_vector(codon_vocab, device)

tai_weight_dict_for_selection = {
    codon: float(tai_weights[idx].detach().cpu().item())
    for codon, idx in codon_vocab.items()
    if codon not in ["<PAD>", "<SOS>", "<EOS>"]
}

dataset_tai_summary = summarize_dataset_tai(
    clean_pairs,
    tai_weight_dict_for_selection,
    output_path="dataset_tai_summary_NbT2T.txt"
)

ALPHA_SOFT_TARGET = 0.7
LAMBDA_ENTROPY = 0.2

print("MAX_SAMPLES:", MAX_SAMPLES)
print("NUM_EPOCHS:", NUM_EPOCHS)
print("BATCH_SIZE:", BATCH_SIZE)
print("ALPHA_SOFT_TARGET:", ALPHA_SOFT_TARGET)
print("LAMBDA_ENTROPY:", LAMBDA_ENTROPY)
print("LAMBDA_TAI:", LAMBDA_TAI)
print("LAMBDA_MFE:", LAMBDA_MFE)
print("TARGET_TAI:", TARGET_TAI)
print("TARGET_MFE_PROXY:", TARGET_MFE_PROXY)
print("RNAfold available:", shutil.which("RNAfold") is not None)
print("USE_TRNASCAN_TAI:", USE_TRNASCAN_TAI)
print("TRNA_SCAN_FILE:", TRNA_SCAN_FILE)
print("USE_POST_MFE_SELECTION:", USE_POST_MFE_SELECTION)
print("POST_TOP_K:", POST_TOP_K)
print("POST_NUM_CANDIDATES:", POST_NUM_CANDIDATES)
print("TAI_DECODING_BONUS:", TAI_DECODING_BONUS)
print("FORBIDDEN_RESTRICTION_SITES:", FORBIDDEN_RESTRICTION_SITES)
print("POLYA_SIGNALS:", POLYA_SIGNALS)
print("GENERATE_RNAFOLD_CANDIDATES:", GENERATE_RNAFOLD_CANDIDATES)
print("RUN_EXACT_RNAFOLD_ON_TEST_SET:", RUN_EXACT_RNAFOLD_ON_TEST_SET)


train_losses = []
val_losses = []
train_ce_losses = []
val_ce_losses = []
train_gc_losses = []
val_gc_losses = []
train_expected_gcs = []
val_expected_gcs = []
train_entropy_losses = []
val_entropy_losses = []
train_tai_losses = []
val_tai_losses = []
train_mfe_proxy_losses = []
val_mfe_proxy_losses = []
train_expected_tais = []
val_expected_tais = []
train_expected_mfe_proxies = []
val_expected_mfe_proxies = []


for epoch in range(NUM_EPOCHS):
    model.train()
    total_train_loss = 0.0
    total_train_ce = 0.0
    total_train_gc = 0.0
    total_train_expected_gc = 0.0
    total_train_entropy = 0.0
    total_train_tai = 0.0
    total_train_mfe_proxy = 0.0
    total_train_expected_tai = 0.0
    total_train_expected_mfe_proxy = 0.0

    for src_batch, tgt_batch, gc_batch in train_loader:
        src_batch = src_batch.to(device)
        tgt_batch = tgt_batch.to(device)
        gc_batch = gc_batch.to(device)

        src_pad_mask = (src_batch == aa_vocab["<PAD>"])
        aligned_target = build_aligned_targets(src_batch, tgt_batch)

        output = model(
            src=src_batch,
            gc_target=gc_batch,
            src_key_padding_mask=src_pad_mask
        )


        mask_for_batch = valid_codon_mask[src_batch]
        masked_output = output.masked_fill(~mask_for_batch, -1e9)

        syn_loss, soft_ce_loss, entropy_loss = compute_synonymous_entropy_loss(
            logits=masked_output,
            src_batch=src_batch,
            aligned_target=aligned_target,
            q_codon_given_aa=q_codon_given_aa,
            target_entropy_by_aa=target_entropy_by_aa,
            alpha=ALPHA_SOFT_TARGET,
            lambda_entropy=LAMBDA_ENTROPY
        )

        gc_loss, avg_expected_gc = compute_gc_loss(
            masked_output,
            aligned_target,
            gc_batch,
            gc_weights
        )

        tai_loss, avg_expected_tai = compute_tai_loss(
            masked_output,
            aligned_target,
            tai_weights,
            target_tai=TARGET_TAI
        )

        mfe_proxy_loss, avg_expected_mfe_proxy = compute_mfe_proxy_loss(
            masked_output,
            aligned_target,
            mfe_proxy_weights,
            target_mfe_proxy=TARGET_MFE_PROXY
        )

        loss = (
            syn_loss
            + LAMBDA_GC * gc_loss
            + LAMBDA_TAI * tai_loss
            + LAMBDA_MFE * mfe_proxy_loss
        )


        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        total_train_loss += loss.item()
        total_train_ce += soft_ce_loss.item()
        total_train_gc += gc_loss.item()
        total_train_expected_gc += avg_expected_gc
        total_train_entropy += entropy_loss.item()
        total_train_tai += tai_loss.item()
        total_train_mfe_proxy += mfe_proxy_loss.item()
        total_train_expected_tai += avg_expected_tai
        total_train_expected_mfe_proxy += avg_expected_mfe_proxy

    avg_train_loss = total_train_loss / len(train_loader)
    avg_train_ce = total_train_ce / len(train_loader)
    avg_train_gc = total_train_gc / len(train_loader)
    avg_train_expected_gc = total_train_expected_gc / len(train_loader)
    avg_train_entropy = total_train_entropy / len(train_loader)
    avg_train_tai = total_train_tai / len(train_loader)
    avg_train_mfe_proxy = total_train_mfe_proxy / len(train_loader)
    avg_train_expected_tai = total_train_expected_tai / len(train_loader)
    avg_train_expected_mfe_proxy = total_train_expected_mfe_proxy / len(train_loader)

    train_losses.append(avg_train_loss)
    train_ce_losses.append(avg_train_ce)
    train_gc_losses.append(avg_train_gc)
    train_expected_gcs.append(avg_train_expected_gc)
    train_entropy_losses.append(avg_train_entropy)
    train_tai_losses.append(avg_train_tai)
    train_mfe_proxy_losses.append(avg_train_mfe_proxy)
    train_expected_tais.append(avg_train_expected_tai)
    train_expected_mfe_proxies.append(avg_train_expected_mfe_proxy)

    model.eval()
    total_val_loss = 0.0
    total_val_ce = 0.0
    total_val_gc = 0.0
    total_val_expected_gc = 0.0
    total_val_entropy = 0.0
    total_val_tai = 0.0
    total_val_mfe_proxy = 0.0
    total_val_expected_tai = 0.0
    total_val_expected_mfe_proxy = 0.0

    with torch.no_grad():
        for src_batch, tgt_batch, gc_batch in val_loader:
            src_batch = src_batch.to(device)
            tgt_batch = tgt_batch.to(device)
            gc_batch = gc_batch.to(device)

            src_pad_mask = (src_batch == aa_vocab["<PAD>"])
            aligned_target = build_aligned_targets(src_batch, tgt_batch)

            output = model(
                src=src_batch,
                gc_target=gc_batch,
                src_key_padding_mask=src_pad_mask
            )

            mask_for_batch = valid_codon_mask[src_batch]
            masked_output = output.masked_fill(~mask_for_batch, -1e9)


            syn_loss, soft_ce_loss, entropy_loss = compute_synonymous_entropy_loss(
                logits=masked_output,
                src_batch=src_batch,
                aligned_target=aligned_target,
                q_codon_given_aa=q_codon_given_aa,
                target_entropy_by_aa=target_entropy_by_aa,
                alpha=ALPHA_SOFT_TARGET,
                lambda_entropy=LAMBDA_ENTROPY
            )


            gc_loss, avg_expected_gc = compute_gc_loss(
                masked_output,
                aligned_target,
                gc_batch,
                gc_weights
            )

            tai_loss, avg_expected_tai = compute_tai_loss(
                masked_output,
                aligned_target,
                tai_weights,
                target_tai=TARGET_TAI
            )

            mfe_proxy_loss, avg_expected_mfe_proxy = compute_mfe_proxy_loss(
                masked_output,
                aligned_target,
                mfe_proxy_weights,
                target_mfe_proxy=TARGET_MFE_PROXY
            )

            loss = (
                syn_loss
                + LAMBDA_GC * gc_loss
                + LAMBDA_TAI * tai_loss
                + LAMBDA_MFE * mfe_proxy_loss
            )

            total_val_loss += loss.item()
            total_val_ce += soft_ce_loss.item()
            total_val_gc += gc_loss.item()
            total_val_expected_gc += avg_expected_gc
            total_val_entropy += entropy_loss.item()
            total_val_tai += tai_loss.item()
            total_val_mfe_proxy += mfe_proxy_loss.item()
            total_val_expected_tai += avg_expected_tai
            total_val_expected_mfe_proxy += avg_expected_mfe_proxy

    avg_val_loss = total_val_loss / len(val_loader)
    avg_val_ce = total_val_ce / len(val_loader)
    avg_val_gc = total_val_gc / len(val_loader)
    avg_val_expected_gc = total_val_expected_gc / len(val_loader)
    avg_val_entropy = total_val_entropy / len(val_loader)
    avg_val_tai = total_val_tai / len(val_loader)
    avg_val_mfe_proxy = total_val_mfe_proxy / len(val_loader)
    avg_val_expected_tai = total_val_expected_tai / len(val_loader)
    avg_val_expected_mfe_proxy = total_val_expected_mfe_proxy / len(val_loader)

    val_losses.append(avg_val_loss)
    val_ce_losses.append(avg_val_ce)
    val_gc_losses.append(avg_val_gc)
    val_expected_gcs.append(avg_val_expected_gc)
    val_entropy_losses.append(avg_val_entropy)
    val_tai_losses.append(avg_val_tai)
    val_mfe_proxy_losses.append(avg_val_mfe_proxy)
    val_expected_tais.append(avg_val_expected_tai)
    val_expected_mfe_proxies.append(avg_val_expected_mfe_proxy)

    print(
        f"Epoch {epoch + 1}/{NUM_EPOCHS}, "
        f"Train Total Loss: {avg_train_loss:.4f}, "
        f"Val Total Loss: {avg_val_loss:.4f}, "
        f"Train Soft CE: {avg_train_ce:.4f}, "
        f"Val Soft CE: {avg_val_ce:.4f}, "
        f"Train Entropy: {avg_train_entropy:.4f}, "
        f"Val Entropy: {avg_val_entropy:.4f}, "
        f"Train GC Loss: {avg_train_gc:.4f}, "
        f"Val GC Loss: {avg_val_gc:.4f}, "
        f"Train tAI Loss: {avg_train_tai:.4f}, "
        f"Val tAI Loss: {avg_val_tai:.4f}, "
        f"Train MFE Proxy Loss: {avg_train_mfe_proxy:.4f}, "
        f"Val MFE Proxy Loss: {avg_val_mfe_proxy:.4f}, "
        f"Train Exp GC: {avg_train_expected_gc:.4f}, "
        f"Val Exp GC: {avg_val_expected_gc:.4f}, "
        f"Train Exp tAI: {avg_train_expected_tai:.4f}, "
        f"Val Exp tAI: {avg_val_expected_tai:.4f}, "
        f"Train Exp MFE Proxy: {avg_train_expected_mfe_proxy:.4f}, "
        f"Val Exp MFE Proxy: {avg_val_expected_mfe_proxy:.4f}"
    )


model.eval()
total_test_loss = 0.0
total_test_ce = 0.0
total_test_gc = 0.0
total_test_expected_gc = 0.0
total_test_entropy = 0.0
total_test_tai = 0.0
total_test_mfe_proxy = 0.0
total_test_expected_tai = 0.0
total_test_expected_mfe_proxy = 0.0
model_correct = 0
model_total = 0

all_true_codons = []
all_pred_codons = []
all_true_aas = []

with torch.no_grad():
    for src_batch, tgt_batch, gc_batch in test_loader:
        src_batch = src_batch.to(device)
        tgt_batch = tgt_batch.to(device)
        gc_batch = gc_batch.to(device)

        src_pad_mask = (src_batch == aa_vocab["<PAD>"])
        aligned_target = build_aligned_targets(src_batch, tgt_batch)

        output = model(
            src=src_batch,
            gc_target=gc_batch,
            src_key_padding_mask=src_pad_mask
        )

        mask_for_batch = valid_codon_mask[src_batch]
        masked_output = output.masked_fill(~mask_for_batch, -1e9)

        syn_loss, soft_ce_loss, entropy_loss = compute_synonymous_entropy_loss(
            logits=masked_output,
            src_batch=src_batch,
            aligned_target=aligned_target,
            q_codon_given_aa=q_codon_given_aa,
            target_entropy_by_aa=target_entropy_by_aa,
            alpha=ALPHA_SOFT_TARGET,
            lambda_entropy=LAMBDA_ENTROPY
        )

        gc_loss, avg_expected_gc = compute_gc_loss(
            masked_output,
            aligned_target,
            gc_batch,
            gc_weights
        )

        tai_loss, avg_expected_tai = compute_tai_loss(
            masked_output,
            aligned_target,
            tai_weights,
            target_tai=TARGET_TAI
        )

        mfe_proxy_loss, avg_expected_mfe_proxy = compute_mfe_proxy_loss(
            masked_output,
            aligned_target,
            mfe_proxy_weights,
            target_mfe_proxy=TARGET_MFE_PROXY
        )

        loss = (
            syn_loss
            + LAMBDA_GC * gc_loss
            + LAMBDA_TAI * tai_loss
            + LAMBDA_MFE * mfe_proxy_loss
        )

        total_test_loss += loss.item()
        total_test_ce += soft_ce_loss.item()
        total_test_gc += gc_loss.item()
        total_test_expected_gc += avg_expected_gc
        total_test_entropy += entropy_loss.item()
        total_test_tai += tai_loss.item()
        total_test_mfe_proxy += mfe_proxy_loss.item()
        total_test_expected_tai += avg_expected_tai
        total_test_expected_mfe_proxy += avg_expected_mfe_proxy

        mask = aligned_target != codon_vocab["<PAD>"]

        for i in range(src_batch.size(0)):
            valid_positions = torch.where(mask[i])[0]

            if USE_POST_MFE_SELECTION:
                pred_codon_list, post_score, post_details = post_mfe_select_codons_for_sequence(
                    logits_seq=masked_output[i],
                    src_seq=src_batch[i],
                    valid_positions=valid_positions,
                    gc_weights=gc_weights,
                    tai_weights=tai_weights,
                    tai_weight_dict=tai_weight_dict_for_selection,
                    top_k=POST_TOP_K,
                    num_candidates=POST_NUM_CANDIDATES
                )
            else:
                adjusted_output = (
                    masked_output[i]
                    + GC_DECODING_BONUS * gc_weights.view(1, -1)
                    + TAI_DECODING_BONUS * tai_weights.view(1, -1)
                )
                pred_indices = adjusted_output.argmax(dim=-1)
                pred_codon_list = [idx_to_codon[pred_indices[pos].item()] for pos in valid_positions]

            for local_idx, j in enumerate(valid_positions.tolist()):
                aa_idx = src_batch[i, j].item()
                true_idx = aligned_target[i, j].item()
                pred_codon = pred_codon_list[local_idx]
                true_codon = idx_to_codon[true_idx]

                all_true_aas.append(idx_to_aa[aa_idx])
                all_true_codons.append(true_codon)
                all_pred_codons.append(pred_codon)

                if pred_codon == true_codon:
                    model_correct += 1
                model_total += 1

avg_test_loss = total_test_loss / len(test_loader)
avg_test_ce = total_test_ce / len(test_loader)
avg_test_gc_loss = total_test_gc / len(test_loader)
avg_test_expected_gc = total_test_expected_gc / len(test_loader)
avg_test_entropy = total_test_entropy / len(test_loader)
avg_test_tai_loss = total_test_tai / len(test_loader)
avg_test_mfe_proxy_loss = total_test_mfe_proxy / len(test_loader)
avg_test_expected_tai = total_test_expected_tai / len(test_loader)
avg_test_expected_mfe_proxy = total_test_expected_mfe_proxy / len(test_loader)
model_accuracy = model_correct / model_total

print(f"\nGC-Controlled Transformer Test Total Loss: {avg_test_loss:.4f}")
print(f"GC-Controlled Transformer Test CE Loss:    {avg_test_ce:.4f}")
print(f"GC-Controlled Transformer Test GC Loss:    {avg_test_gc_loss:.4f}")
print(f"GC-Controlled Transformer Expected GC:     {avg_test_expected_gc:.4f}")
print(f"GC-Controlled Transformer Test Accuracy:   {model_accuracy:.4f}")


def decode_aa_indices(src_tensor):
    aa_list = []
    for idx in src_tensor.tolist():
        token = idx_to_aa[idx]
        if token not in ["<SOS>", "<EOS>", "<PAD>"]:
            aa_list.append(token)
    return aa_list


def decode_codon_indices(tgt_tensor):
    codon_list = []
    for idx in tgt_tensor.tolist():
        token = idx_to_codon[idx]
        if token not in ["<SOS>", "<EOS>", "<PAD>"]:
            codon_list.append(token)
    return codon_list


def build_most_frequent_codon_table(train_dataset):
    aa_to_codon_counts = defaultdict(Counter)

    for i in range(len(train_dataset)):
        src_tensor, tgt_tensor, _ = train_dataset[i]
        aa_seq = decode_aa_indices(src_tensor)
        codon_seq = decode_codon_indices(tgt_tensor)

        assert len(aa_seq) == len(codon_seq), f"Length mismatch at item {i}"

        for aa, codon in zip(aa_seq, codon_seq):
            aa_to_codon_counts[aa][codon] += 1

    baseline_table = {}
    for aa, counter in aa_to_codon_counts.items():
        baseline_table[aa] = counter.most_common(1)[0][0]

    return baseline_table, aa_to_codon_counts


def evaluate_most_frequent_codon_baseline(test_dataset, baseline_table):
    correct = 0
    total = 0

    true_codons = []
    pred_codons = []
    true_aas = []

    for i in range(len(test_dataset)):
        src_tensor, tgt_tensor, _ = test_dataset[i]

        aa_seq = decode_aa_indices(src_tensor)
        true_codon_seq = decode_codon_indices(tgt_tensor)

        assert len(aa_seq) == len(true_codon_seq), f"Length mismatch at item {i}"

        pred_codon_seq = [baseline_table[aa] for aa in aa_seq]

        for aa, pred, true in zip(aa_seq, pred_codon_seq, true_codon_seq):
            true_aas.append(aa)
            true_codons.append(true)
            pred_codons.append(pred)

            if pred == true:
                correct += 1
            total += 1

    accuracy = correct / total
    return accuracy, true_aas, true_codons, pred_codons


baseline_table, aa_to_codon_counts = build_most_frequent_codon_table(train_dataset)

print("\nMost-frequent-codon table:")
for aa in sorted(baseline_table.keys()):
    print(f"{aa} -> {baseline_table[aa]}")

baseline_accuracy, baseline_true_aas, baseline_true_codons, baseline_pred_codons = \
    evaluate_most_frequent_codon_baseline(test_dataset, baseline_table)

print(f"\nMost-Frequent-Codon Baseline Accuracy: {baseline_accuracy:.4f}")


def accuracy_by_amino_acid(true_aas, true_codons, pred_codons):
    aa_correct = Counter()
    aa_total = Counter()

    for aa, true_codon, pred_codon in zip(true_aas, true_codons, pred_codons):
        if pred_codon == true_codon:
            aa_correct[aa] += 1
        aa_total[aa] += 1

    aa_acc = {}
    for aa in sorted(aa_total.keys()):
        aa_acc[aa] = aa_correct[aa] / aa_total[aa]

    return aa_acc


sense_codons = [c for c in all_codons if genetic_code[c] != "*"]


def codon_frequency_vector(codon_list, codon_order):
    counts = Counter(codon_list)
    vec = np.array([counts[c] for c in codon_order], dtype=float)
    if vec.sum() > 0:
        vec = vec / vec.sum()
    return vec


def cosine_similarity_vec(a, b):
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return 0.0
    return float(np.dot(a, b) / denom)


def gc_content_from_codons(codon_list):
    seq = "".join(codon_list)
    if len(seq) == 0:
        return 0.0
    gc = sum(1 for b in seq if b in ["G", "C"])
    return gc / len(seq)


def build_codon_weights_from_training(train_dataset):
    aa_to_counts = defaultdict(Counter)

    for i in range(len(train_dataset)):
        src_tensor, tgt_tensor, _ = train_dataset[i]
        aa_seq = decode_aa_indices(src_tensor)
        codon_seq = decode_codon_indices(tgt_tensor)

        for aa, codon in zip(aa_seq, codon_seq):
            aa_to_counts[aa][codon] += 1

    codon_weights = {}
    for aa, codons in aa_to_synonymous_codons.items():
        max_count = max(aa_to_counts[aa][c] for c in codons)
        for c in codons:
            codon_weights[c] = aa_to_counts[aa][c] / max_count if max_count > 0 else 0.0

    return codon_weights


def cai_of_codon_sequence(codon_list, codon_weights):
    weights = []
    for codon in codon_list:
        if codon in codon_weights and codon_weights[codon] > 0:
            weights.append(codon_weights[codon])

    if len(weights) == 0:
        return 0.0

    log_mean = np.mean(np.log(weights))
    return float(np.exp(log_mean))


def codons_to_amino_acids(codon_list):
    return [genetic_code.get(codon, "?") for codon in codon_list]


def simple_accuracy(true_list, pred_list):
    correct = sum(t == p for t, p in zip(true_list, pred_list))
    total = len(true_list)
    return correct / total if total > 0 else 0.0


model_aa_acc = accuracy_by_amino_acid(all_true_aas, all_true_codons, all_pred_codons)
baseline_aa_acc = accuracy_by_amino_acid(
    baseline_true_aas,
    baseline_true_codons,
    baseline_pred_codons
)

true_codon_vec = codon_frequency_vector(all_true_codons, sense_codons)
model_codon_vec = codon_frequency_vector(all_pred_codons, sense_codons)
baseline_codon_vec = codon_frequency_vector(baseline_pred_codons, sense_codons)

model_codon_cosine = cosine_similarity_vec(true_codon_vec, model_codon_vec)
baseline_codon_cosine = cosine_similarity_vec(true_codon_vec, baseline_codon_vec)

true_gc = gc_content_from_codons(all_true_codons)
model_gc = gc_content_from_codons(all_pred_codons)
baseline_gc = gc_content_from_codons(baseline_pred_codons)

model_gc_diff = abs(model_gc - true_gc)
baseline_gc_diff = abs(baseline_gc - true_gc)

codon_weights = build_codon_weights_from_training(train_dataset)

tai_weight_dict = tai_weight_dict_for_selection

true_cai = cai_of_codon_sequence(all_true_codons, codon_weights)
model_cai = cai_of_codon_sequence(all_pred_codons, codon_weights)
baseline_cai = cai_of_codon_sequence(baseline_pred_codons, codon_weights)

model_cai_diff = abs(model_cai - true_cai)
baseline_cai_diff = abs(baseline_cai - true_cai)

true_tai = tai_of_codon_sequence(all_true_codons, tai_weight_dict)
model_tai = tai_of_codon_sequence(all_pred_codons, tai_weight_dict)
baseline_tai = tai_of_codon_sequence(baseline_pred_codons, tai_weight_dict)
model_tai_diff = abs(model_tai - true_tai)
baseline_tai_diff = abs(baseline_tai - true_tai)

true_mfe_proxy = mfe_proxy_of_codon_sequence(all_true_codons)
model_mfe_proxy = mfe_proxy_of_codon_sequence(all_pred_codons)
baseline_mfe_proxy = mfe_proxy_of_codon_sequence(baseline_pred_codons)
model_mfe_proxy_diff = abs(model_mfe_proxy - true_mfe_proxy)
baseline_mfe_proxy_diff = abs(baseline_mfe_proxy - true_mfe_proxy)

if RUN_EXACT_RNAFOLD_ON_TEST_SET:
    true_rna = "".join(all_true_codons).replace("T", "U")
    model_rna = "".join(all_pred_codons).replace("T", "U")
    baseline_rna = "".join(baseline_pred_codons).replace("T", "U")
    true_rnafold_mfe = run_rnafold_mfe(true_rna)
    model_rnafold_mfe = run_rnafold_mfe(model_rna)
    baseline_rnafold_mfe = run_rnafold_mfe(baseline_rna)
else:
    true_rnafold_mfe = None
    model_rnafold_mfe = None
    baseline_rnafold_mfe = None

true_aa_from_codons = codons_to_amino_acids(all_true_codons)
model_aa_from_codons = codons_to_amino_acids(all_pred_codons)
baseline_aa_from_codons = codons_to_amino_acids(baseline_pred_codons)

model_amino_acid_accuracy = simple_accuracy(true_aa_from_codons, model_aa_from_codons)
baseline_amino_acid_accuracy = simple_accuracy(true_aa_from_codons, baseline_aa_from_codons)


print("\nGC meatures:")
print(f"True GC:                                        {true_gc:.4f}")
print(f"GC-controlled Transformer GC:                   {model_gc:.4f}")
print(f"Baseline GC:                                    {baseline_gc:.4f}")
print(f"Transformer GC difference from target:          {abs(model_gc - TARGET_GC):.4f}")
print(f"Baseline GC difference from target:             {abs(baseline_gc - TARGET_GC):.4f}")

print("\ntAI-like measures from tRNAscan-SE gene-copy weights:")
print(f"True tAI-like score:                                 {true_tai:.4f}")
print(f"GC-controlled Transformer tAI-like score:            {model_tai:.4f}")
print(f"Baseline tAI-like score:                             {baseline_tai:.4f}")
print(f"Transformer tAI-like difference from true:     {model_tai_diff:.4f}")
print(f"Baseline tAI-like difference from true:        {baseline_tai_diff:.4f}")

print("\nMFE-proxy measures:")
print(f"True MFE proxy:                                 {true_mfe_proxy:.4f}")
print(f"GC-controlled Transformer MFE proxy:            {model_mfe_proxy:.4f}")
print(f"Baseline MFE proxy:                             {baseline_mfe_proxy:.4f}")
print(f"Transformer MFE-proxy difference from true:     {model_mfe_proxy_diff:.4f}")
print(f"Baseline MFE-proxy difference from true:        {baseline_mfe_proxy_diff:.4f}")

if true_rnafold_mfe is not None:
    print("\nRNAfold exact MFE measures:")
    print(f"True RNAfold MFE:                               {true_rnafold_mfe:.4f}")
    print(f"GC-controlled Transformer RNAfold MFE:          {model_rnafold_mfe:.4f}")
    print(f"Baseline RNAfold MFE:                           {baseline_rnafold_mfe:.4f}")
else:
    print("\nRNAfold exact MFE measures: skipped for this run; using MFE proxy only.")

print(f"Transformer amino acid fidelity:                {model_amino_acid_accuracy:.4f}")
print(f"Baseline amino acid fidelity:                   {baseline_amino_acid_accuracy:.4f}")

print("\n=== FINAL CHECK SUMMARY ===")
print(f"Test Total Loss:          {avg_test_loss:.4f}")
print(f"Test Soft CE Loss:        {avg_test_ce:.4f}")
print(f"Test Entropy Loss:        {avg_test_entropy:.4f}")
print(f"Test GC Loss:             {avg_test_gc_loss:.4f}")
print(f"Test tAI Loss:            {avg_test_tai_loss:.4f}")
print(f"Test MFE Proxy Loss:      {avg_test_mfe_proxy_loss:.4f}")
print(f"Test Expected GC:         {avg_test_expected_gc:.4f}")
print(f"Test Expected tAI:        {avg_test_expected_tai:.4f}")
print(f"Test Expected MFE Proxy:  {avg_test_expected_mfe_proxy:.4f}")
print(f"Target GC:                {TARGET_GC:.4f}")
print(f"Model Actual GC:          {model_gc:.4f}")
print(f"GC Difference to Target:  {abs(model_gc - TARGET_GC):.4f}")
print(f"AA Fidelity:              {model_amino_acid_accuracy:.4f}")
print(f"Codon Accuracy:           {model_accuracy:.4f}")
print(f"Codon Cosine Similarity:  {model_codon_cosine:.4f}")
print(f"CAI Difference:           {model_cai_diff:.4f}")
print(f"tAI-like Difference:     {model_tai_diff:.4f}")
print(f"MFE Proxy Difference:     {model_mfe_proxy_diff:.4f}")


def generate_ranked_candidates_for_sequence(
    logits_seq,
    valid_positions,
    gc_weights,
    tai_weights,
    tai_weight_dict,
    codon_weights,
    top_k=4,
    num_candidates=12
):
    adjusted_logits = (
        logits_seq
        + GC_DECODING_BONUS * gc_weights.view(1, -1)
        + TAI_DECODING_BONUS * tai_weights.view(1, -1)
    )

    position_choices = []
    for pos in valid_positions.tolist():
        scores = adjusted_logits[pos]
        finite_idx = torch.where(scores > -1e8)[0]
        if finite_idx.numel() == 0:
            position_choices.append([])
            continue
        k = min(top_k, finite_idx.numel())
        top_values, top_indices_local = torch.topk(scores[finite_idx], k=k)
        top_indices = finite_idx[top_indices_local]
        position_choices.append(top_indices.detach().cpu().tolist())

    raw_candidates = []

    greedy = []
    for choices in position_choices:
        if choices:
            greedy.append(idx_to_codon[choices[0]])
    raw_candidates.append(greedy)

    for _ in range(max(0, num_candidates - 1)):
        cand = []
        for choices in position_choices:
            if not choices:
                continue
            chosen_idx = choices[0] if len(choices) == 1 else random.choice(choices)
            cand.append(idx_to_codon[chosen_idx])
        raw_candidates.append(cand)

    unique = []
    seen = set()
    for cand in raw_candidates:
        seq = "".join(cand)
        if seq and seq not in seen:
            seen.add(seq)
            unique.append(cand)

    ranked = []
    for cand in unique:
        dna = "".join(cand)
        gc = gc_content_from_codons(cand)
        tai = tai_of_codon_sequence(cand, tai_weight_dict)
        cai = cai_of_codon_sequence(cand, codon_weights)
        mfe_proxy = mfe_proxy_of_codon_sequence(cand)
        motif_penalty, restriction_hits, polya_hits, restriction_total, polya_total = motif_penalty_for_codons(cand)

        score = (
            4.0 * abs(gc - TARGET_GC)
            + 1.0 * abs(tai - TARGET_TAI)
            - 0.5 * cai
            + motif_penalty
        )

        ranked.append({
            "codons": cand,
            "dna": dna,
            "score": score,
            "gc": gc,
            "tai_like": tai,
            "cai": cai,
            "mfe_proxy": mfe_proxy,
            "restriction_site_total": restriction_total,
            "polyA_signal_total": polya_total,
            "restriction_hits": restriction_hits,
            "polyA_hits": polya_hits,
        })

    ranked.sort(key=lambda x: x["score"])
    return ranked


if GENERATE_RNAFOLD_CANDIDATES:
    print("\nGenerating candidate CDS sequences for later RNAfold...")
    model.eval()
    fasta_records = []
    tsv_rows = []
    selected_sequences = 0

    with torch.no_grad():
        for src_batch, tgt_batch, gc_batch in test_loader:
            src_batch = src_batch.to(device)
            tgt_batch = tgt_batch.to(device)
            gc_batch = gc_batch.to(device)

            aligned_target = build_aligned_targets(src_batch, tgt_batch)
            src_pad_mask = (src_batch == aa_vocab["<PAD>"])

            output = model(
                src=src_batch,
                gc_target=gc_batch,
                src_key_padding_mask=src_pad_mask
            )
            mask_for_batch = valid_codon_mask[src_batch]
            masked_output = output.masked_fill(~mask_for_batch, -1e9)
            valid_mask = aligned_target != codon_vocab["<PAD>"]

            for i in range(src_batch.size(0)):
                if selected_sequences >= NUM_SEQUENCES_FOR_CANDIDATES:
                    break

                valid_positions = torch.where(valid_mask[i])[0]
                true_codons = [idx_to_codon[aligned_target[i, pos].item()] for pos in valid_positions.tolist()]
                true_aas = [idx_to_aa[src_batch[i, pos].item()] for pos in valid_positions.tolist()]

                ranked_candidates = generate_ranked_candidates_for_sequence(
                    logits_seq=masked_output[i],
                    valid_positions=valid_positions,
                    gc_weights=gc_weights,
                    tai_weights=tai_weights,
                    tai_weight_dict=tai_weight_dict,
                    codon_weights=codon_weights,
                    top_k=CANDIDATE_TOP_K,
                    num_candidates=CANDIDATES_PER_SEQUENCE
                )

                for cand_rank, cand in enumerate(ranked_candidates, start=1):
                    candidate_id = f"seq{selected_sequences + 1}_candidate{cand_rank}"
                    pred_aas = codons_to_amino_acids(cand["codons"])
                    aa_fidelity = simple_accuracy(true_aas, pred_aas)
                    codon_acc = simple_accuracy(true_codons, cand["codons"])

                    fasta_records.append((candidate_id, cand["dna"]))
                    tsv_rows.append({
                        "candidate_id": candidate_id,
                        "source_sequence_index": selected_sequences + 1,
                        "rank": cand_rank,
                        "length_nt": len(cand["dna"]),
                        "score_without_rnafold": cand["score"],
                        "gc": cand["gc"],
                        "gc_difference_to_target": abs(cand["gc"] - TARGET_GC),
                        "tai_like": cand["tai_like"],
                        "tai_difference_to_target": abs(cand["tai_like"] - TARGET_TAI),
                        "cai": cand["cai"],
                        "mfe_proxy_report_only": cand["mfe_proxy"],
                        "restriction_site_total": cand["restriction_site_total"],
                        "polyA_signal_total": cand["polyA_signal_total"],
                        "restriction_hits": cand["restriction_hits"],
                        "polyA_hits": cand["polyA_hits"],
                        "aa_fidelity": aa_fidelity,
                        "codon_accuracy_vs_reference": codon_acc,
                    })

                selected_sequences += 1

            if selected_sequences >= NUM_SEQUENCES_FOR_CANDIDATES:
                break

    with open(RNAFOLD_CANDIDATE_FASTA, "w") as fasta_out:
        for candidate_id, dna in fasta_records:
            fasta_out.write(f">{candidate_id}\n")
            for start in range(0, len(dna), 80):
                fasta_out.write(dna[start:start + 80] + "\n")

    with open(RNAFOLD_CANDIDATE_TSV, "w") as tsv_out:
        header = [
            "candidate_id",
            "source_sequence_index",
            "rank",
            "length_nt",
            "score_without_rnafold",
            "gc",
            "gc_difference_to_target",
            "tai_like",
            "tai_difference_to_target",
            "cai",
            "mfe_proxy_report_only",
            "restriction_site_total",
            "polyA_signal_total",
            "restriction_hits",
            "polyA_hits",
            "aa_fidelity",
            "codon_accuracy_vs_reference",
        ]
        tsv_out.write("\t".join(header) + "\n")
        for row in tsv_rows:
            tsv_out.write("\t".join(str(row[h]) for h in header) + "\n")

    print(f"Saved candidate FASTA for later RNAfold: {RNAFOLD_CANDIDATE_FASTA}")
    print(f"Saved candidate summary: {RNAFOLD_CANDIDATE_TSV}")
    print("RNAfold was not run in this script. Run RNAfold later on the saved FASTA candidates.")


epochs = list(range(1, len(train_losses) + 1))

plt.figure(figsize=(8, 5))
plt.plot(epochs, train_losses, label="Train Total Loss")
plt.plot(epochs, val_losses, label="Validation Total Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("GC-Controlled Transformer Total Loss")
plt.legend()
plt.tight_layout()
plt.savefig("gc_controlled_total_loss_NbT2T.png")
plt.close()

plt.figure(figsize=(8, 5))
plt.plot(epochs, train_ce_losses, label="Train CE Loss")
plt.plot(epochs, val_ce_losses, label="Validation CE Loss")
plt.xlabel("Epoch")
plt.ylabel("Cross-Entropy Loss")
plt.title("GC-Controlled Transformer Cross-Entropy Loss")
plt.legend()
plt.tight_layout()
plt.savefig("gc_controlled_ce_loss_NbT2T.png")
plt.close()


torch.save(model.state_dict(), "transformer_model_NbT2T.pt")

with open("predictions_NbT2T.pkl", "wb") as f:
    pickle.dump({
        "true_codons": all_true_codons,
        "pred_codons": all_pred_codons,
        "true_aas": all_true_aas
    }, f)

with open("training_history_NbT2T.pkl", "wb") as f:
    pickle.dump({
        "train_losses": train_losses,
        "val_losses": val_losses,
        "train_ce_losses": train_ce_losses,
        "val_ce_losses": val_ce_losses,
        "train_gc_losses": train_gc_losses,
        "val_gc_losses": val_gc_losses,
        "train_expected_gcs": train_expected_gcs,
        "val_expected_gcs": val_expected_gcs,
        "train_entropy_losses": train_entropy_losses,
        "val_entropy_losses": val_entropy_losses,
        "train_tai_losses": train_tai_losses,
        "val_tai_losses": val_tai_losses,
        "train_mfe_proxy_losses": train_mfe_proxy_losses,
        "val_mfe_proxy_losses": val_mfe_proxy_losses,
        "train_expected_tais": train_expected_tais,
        "val_expected_tais": val_expected_tais,
        "train_expected_mfe_proxies": train_expected_mfe_proxies,
        "val_expected_mfe_proxies": val_expected_mfe_proxies
    }, f)

with open("metrics_NbT2T.txt", "w") as f:
    f.write("=== FINAL RESULTS ===\n")
    f.write(f"GC-controlled Transformer Total Loss: {avg_test_loss:.4f}\n")
    f.write(f"GC-controlled Transformer CE Loss: {avg_test_ce:.4f}\n")
    f.write(f"GC-controlled Transformer GC Loss: {avg_test_gc_loss:.4f}\n")
    f.write(f"GC-controlled Transformer Expected GC: {avg_test_expected_gc:.4f}\n")
    f.write(f"GC-controlled Transformer Expected tAI: {avg_test_expected_tai:.4f}\n")
    f.write(f"GC-controlled Transformer Expected MFE Proxy: {avg_test_expected_mfe_proxy:.4f}\n")
    f.write(f"Target GC: {TARGET_GC:.4f}\n")
    f.write(f"Transformer Accuracy: {model_accuracy:.4f}\n")
    f.write(f"Baseline Accuracy: {baseline_accuracy:.4f}\n")
    f.write(f"Codon Cosine: {model_codon_cosine:.4f}\n")
    f.write(f"GC Difference: {model_gc_diff:.4f}\n")
    f.write(f"CAI Difference: {model_cai_diff:.4f}\n")
    f.write(f"tAI-like Difference: {model_tai_diff:.4f}\n")
    f.write(f"MFE Proxy Difference: {model_mfe_proxy_diff:.4f}\n")
    f.write(f"Amino Acid Fidelity: {model_amino_acid_accuracy:.4f}\n")
    f.write(f"USE_TRNASCAN_TAI: {USE_TRNASCAN_TAI}\n")
    f.write(f"TRNA_SCAN_FILE: {TRNA_SCAN_FILE}\n")
    f.write(f"USE_POST_MFE_SELECTION: {USE_POST_MFE_SELECTION}\n")
    f.write(f"POST_TOP_K: {POST_TOP_K}\n")
    f.write(f"POST_NUM_CANDIDATES: {POST_NUM_CANDIDATES}\n")
    f.write(f"LAMBDA_MFE: {LAMBDA_MFE}\n")
    f.write(f"RUN_EXACT_RNAFOLD_ON_TEST_SET: {RUN_EXACT_RNAFOLD_ON_TEST_SET}\n")
    f.write(f"FORBIDDEN_RESTRICTION_SITES: {FORBIDDEN_RESTRICTION_SITES}\n")
    f.write(f"POLYA_SIGNALS: {POLYA_SIGNALS}\n")
    f.write(f"GENERATE_RNAFOLD_CANDIDATES: {GENERATE_RNAFOLD_CANDIDATES}\n")
    f.write(f"RNAFOLD_CANDIDATE_FASTA: {RNAFOLD_CANDIDATE_FASTA}\n")
    f.write(f"RNAFOLD_CANDIDATE_TSV: {RNAFOLD_CANDIDATE_TSV}\n")
