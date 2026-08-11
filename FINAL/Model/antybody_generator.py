#!/usr/bin/env python3
"""
Generate codon-optimized CDS candidates for Adalimumab heavy/light chains
using the current trained NbT2T encoder-only Transformer model.

DETAILED NOTES:
- Candidate quality is evaluated with GC content, a tRNA-derived tAI-like score,
  an MFE proxy/RNAfold result when available, and selected sequence-motif penalties.
- The Transformer is encoder-only. It receives the amino-acid sequence plus a
  target GC value and predicts a codon distribution at each amino-acid position.
- One greedy candidate is generated from the highest-scoring codon at each
  position. Additional candidates are probability-weighted samples from the
  top-K valid codons, providing diversity without allowing arbitrary codons.


Required files in the same folder:
  - NbT2T_model.pt  OR  codon_model_NbT2T.pt  OR  model_NbT2T.pt
  - clean_aa_codon_pairs_NbT2T.txt
  - nbenthamiana_NbT2T_trna.out

Run:
  python antibody_generator.py

Outputs:
  - adalimumab_generated_candidates_NbT2T.fasta
  - adalimumab_generated_candidates_NbT2T.tsv
"""

import math
import random
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn


# Settings

# the GC content can be changed here for 41-48 per cent depends what you need the weights on select5ion

RANDOM_SEED = 1234
TARGET_GC = 0.44
TARGET_TAI = 0.3765
TARGET_MFE_PROXY = 0.48
GC_FLOOR = 0.43
GC_CEILING = 0.48
GC_CEILING_WEIGHT = 50.0
POST_GC_WEIGHT = 8.0

GC_DECODING_BONUS = 0.0
TAI_DECODING_BONUS = 0.0
TOP_K = 4
NUM_CANDIDATES_PER_CHAIN = 120
SAMPLING_TEMPERATURE = 0.8

TRNA_SCAN_FILE = "nbenthamiana_NbT2T_trna.out"
TRAINING_PAIR_FILE = "clean_aa_codon_pairs_NbT2T.txt"


CHECKPOINT_CANDIDATES = [
    "transformer_model_NbT2T.pt",
    "NbT2T_model.pt",
    "codon_model_NbT2T.pt",
    "model_NbT2T.pt",
    "trained_model_NbT2T.pt",
]

OUTPUT_FASTA = "antibody_candidates.fasta"
OUTPUT_TSV = "antibody_candidates.tsv"

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

random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)


# Target antibody AA sequences
# Stop symbol removed. Signal peptide kept.


TARGET_PROTEINS = {
    "Adalimumab_HC": (
        "MAKTNLFLFLIFSLLLSLSSAEVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGK"
        "GLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSL"
        "DYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHT"
        "FPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLG"
        "GPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTY"
        "RVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSL"
        "TCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMH"
        "EALHNHYTQKSLSLSPGK"
    ),
    "Adalimumab_LC": (
        "MAKTNLFLFLIFSLLLSLSSADIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPK"
        "LLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIKRTVAAP"
        "SVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTL"
        "TLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
    ),
}


# Vocab and genetic code



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
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

aa_to_synonymous_codons = defaultdict(list)
for codon, aa in genetic_code.items():
    if aa != "*":
        aa_to_synonymous_codons[aa].append(codon)

bases = ["T", "C", "A", "G"]
all_codons = [a + b + c for a in bases for b in bases for c in bases]

codon_vocab = {"<PAD>": 0, "<SOS>": 1, "<EOS>": 2}
for i, codon in enumerate(all_codons, start=3):
    codon_vocab[codon] = i
idx_to_codon = {v: k for k, v in codon_vocab.items()}

aa_vocab = {
    "<PAD>": 0, "<SOS>": 1, "<EOS>": 2,
    "A": 3, "R": 4, "N": 5, "D": 6, "C": 7,
    "Q": 8, "E": 9, "G": 10, "H": 11, "I": 12,
    "L": 13, "K": 14, "M": 15, "F": 16, "P": 17,
    "S": 18, "T": 19, "W": 20, "Y": 21, "V": 22,
}
idx_to_aa = {v: k for k, v in aa_vocab.items()}

# Model definition: must match training script
# IMPORTANT: this architecture must match the architecture used to create the checkpoint. Embedding size, attention heads, encoder layers, feed-forward size, and vocabulary dimensions determine the shape/meaning of the saved parameters.

class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=6000):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer("pe", pe.unsqueeze(0))

    def forward(self, x):
        return x + self.pe[:, : x.size(1), :]


class EncoderOnlyCodonTransformer(nn.Module):
    def __init__(self, src_vocab_size, tgt_vocab_size, d_model=128, nhead=4,
                 num_encoder_layers=4, dim_feedforward=128, dropout=0.1):
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
            batch_first=True,
        )
        self.encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers=num_encoder_layers,
            enable_nested_tensor=False,
        )
        self.fc_out = nn.Linear(d_model, tgt_vocab_size)

    def forward(self, src, gc_target, src_key_padding_mask=None):
        x = self.src_embedding(src) * math.sqrt(self.d_model)
        gc_embed = self.gc_projection(gc_target).unsqueeze(1).expand(-1, src.size(1), -1)
        x = x + gc_embed
        x = self.positional_encoding(x)
        x = self.encoder(x, src_key_padding_mask=src_key_padding_mask)
        return self.fc_out(x)


# Helper functions:
# The helper layer separates validation/encoding, tRNA processing, sequence metrics, motif detection, scoring, model loading, decoding, and file output.


def valid_aa(seq):
    """Validate that a protein sequence contains only standard amino-acid symbols"""
    return all(a in aa_vocab and a not in ["<PAD>", "<SOS>", "<EOS>"] for a in seq)


def encode_aa(seq):
    """Convert an amino-acid sequence into SOS + amino acids + EOS token IDs"""
    return [aa_vocab["<SOS>"]] + [aa_vocab[a] for a in seq] + [aa_vocab["<EOS>"]]


def build_valid_codon_mask(device):
    """Create an amino-acid/codon compatibility mask"""
    valid_mask = torch.zeros(len(aa_vocab), len(codon_vocab), dtype=torch.bool, device=device)
    for aa, codons in aa_to_synonymous_codons.items():
        aa_idx = aa_vocab[aa]
        for codon in codons:
            valid_mask[aa_idx, codon_vocab[codon]] = True
    return valid_mask


def build_gc_weight_vector(device):
    """Calculate the fractional GC content of every codon"""
    gc_weights = torch.zeros(len(codon_vocab), device=device)
    for codon, idx in codon_vocab.items():
        if codon not in ["<PAD>", "<SOS>", "<EOS>"]:
            gc_weights[idx] = sum(1 for b in codon if b in "GC") / 3.0
    return gc_weights


def reverse_complement(seq):
    """Return the DNA reverse complement of a nucleotide sequence"""
    return seq.upper().translate(str.maketrans("ATCG", "TAGC"))[::-1]


aa1_to_aa3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}
aa3_to_aa1 = {v: k for k, v in aa1_to_aa3.items()}
aa3_to_aa1["iMet"] = "M"


def parse_trnascan_to_codon_counts(trna_file):
    """Parse tRNAscan output and count valid tRNA codons"""
    codon_counts = Counter()
    standard_types = set(aa3_to_aa1.keys())
    with open(trna_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("Sequence", "Name", "-")):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            trna_type = parts[4]
            anticodon = parts[5].upper()
            if "pseudo" in line.lower():
                continue
            if trna_type not in standard_types:
                continue
            if len(anticodon) != 3 or any(b not in "ATCG" for b in anticodon):
                continue
            aa1 = aa3_to_aa1[trna_type]
            codon = reverse_complement(anticodon)
            if genetic_code.get(codon) != aa1:
                continue
            codon_counts[codon] += 1
    return codon_counts


def build_tai_weights_from_trnascan(trna_file, device):
    """Convert tRNA codon counts into relative synonymous-codon weights"""
    if not Path(trna_file).exists():
        raise FileNotFoundError(f"Missing {trna_file}. Put your tRNAscan output in this folder.")
    codon_counts = parse_trnascan_to_codon_counts(trna_file)
    pseudocount = 0.5
    weights = torch.zeros(len(codon_vocab), device=device)
    for aa, codons in aa_to_synonymous_codons.items():
        max_count = max(codon_counts[c] for c in codons)
        if max_count == 0:
            for codon in codons:
                weights[codon_vocab[codon]] = 1.0
        else:
            for codon in codons:
                weights[codon_vocab[codon]] = (codon_counts[codon] + pseudocount) / (max_count + pseudocount)
    return weights


def tai_of_codon_sequence(codon_list, tai_weight_dict):
    """Calculate the geometric-mean tAI-like score of a codon sequence"""
    weights = [tai_weight_dict.get(c, 0.0) for c in codon_list]
    weights = [w for w in weights if w > 0]
    if not weights:
        return 0.0
    return float(np.exp(np.mean(np.log(weights))))


def mfe_proxy_of_codon_sequence(codon_list):
    """Calculate a simple codon-composition proxy for RNA folding tendency"""
    if not codon_list:
        return 0.0
    scores = []
    for codon in codon_list:
        rna = codon.replace("T", "U")
        gc = sum(1 for b in rna if b in "GC") / 3.0
        au = sum(1 for b in rna if b in "AU") / 3.0
        scores.append(0.7 * gc + 0.3 * au)
    return float(np.mean(scores))


def run_rnafold_mfe(dna, max_len=5000):
    """Run ViennaRNA RNAfold when available and return its parsed MFE"""
    if shutil.which("RNAfold") is None or len(dna) > max_len:
        return None
    rna = dna.replace("T", "U")
    result = subprocess.run(["RNAfold", "--noPS"], input=rna + "\n", text=True, capture_output=True, check=False)
    if result.returncode != 0:
        return None
    match = re.search(r"\(([-+]?[0-9]*\.?[0-9]+)\)", result.stdout)
    return float(match.group(1)) if match else None


def count_motif_hits(dna_sequence, motif_dict):
    """Count occurrences of each motif and return both detail and total counts"""
    hits = {}
    total = 0
    for name, motif in motif_dict.items():
        count = dna_sequence.count(motif)
        hits[name] = count
        total += count
    return hits, total


def sequence_metrics(codon_list, tai_weight_dict):
    """Compute the complete set of post-decoding candidate metrics"""
    dna = "".join(codon_list)
    gc = sum(1 for b in dna if b in "GC") / len(dna) if dna else 0.0
    tai = tai_of_codon_sequence(codon_list, tai_weight_dict)
    mfe_proxy = mfe_proxy_of_codon_sequence(codon_list)
    rnafold_mfe = run_rnafold_mfe(dna)
    r_hits, r_total = count_motif_hits(dna, FORBIDDEN_RESTRICTION_SITES)
    p_hits, p_total = count_motif_hits(dna, POLYA_SIGNALS)
    return {
        "dna": dna,
        "length_nt": len(dna),
        "gc": gc,
        "tai": tai,
        "mfe_proxy": mfe_proxy,
        "rnafold_mfe": rnafold_mfe,
        "rnafold_mfe_per_nt": (rnafold_mfe / len(dna)) if rnafold_mfe is not None and dna else None,
        "restriction_site_total": r_total,
        "polyA_signal_total": p_total,
        "restriction_hits": r_hits,
        "polyA_hits": p_hits,
    }


def selection_score(metrics):
    """Combine candidate-quality penalties into one lower-is-better score"""

    mfe_term = abs(metrics["mfe_proxy"] - TARGET_MFE_PROXY)
    gc_ceiling_penalty = max(0.0, metrics["gc"] - GC_CEILING)
    gc_floor_penalty = max(0.0, GC_FLOOR - metrics["gc"])
    return (
        POST_GC_WEIGHT * abs(metrics["gc"] - TARGET_GC)
        + GC_CEILING_WEIGHT * gc_ceiling_penalty
        + GC_CEILING_WEIGHT * gc_floor_penalty
        + 1.0 * abs(metrics["tai"] - TARGET_TAI)
        + 0.3 * mfe_term
        + 10.0 * metrics["restriction_site_total"]
        + 10.0 * metrics["polyA_signal_total"]
    )


def find_checkpoint():
    """Return the first existing trained-model checkpoint from the candidate list"""
    for name in CHECKPOINT_CANDIDATES:
        if Path(name).exists():
            return name
    raise FileNotFoundError(
        "No model checkpoint found. Save your trained model first, for example:\n"
        "torch.save(model.state_dict(), 'NbT2T_model.pt')\n"
        f"Then place it in this folder. Tried: {CHECKPOINT_CANDIDATES}"
    )


def decode_candidates(model, protein_name, aa_seq, device, valid_codon_mask, gc_weights, tai_weights, tai_weight_dict):
    """Generate synonymous DNA candidates for one protein sequence."""
    if aa_seq.endswith("*") or aa_seq.endswith("."):
        aa_seq = aa_seq[:-1]
    aa_seq = aa_seq.replace("\n", "").replace(" ", "").upper()
    if not valid_aa(aa_seq):
        bad = sorted(set(a for a in aa_seq if a not in aa_vocab))
        raise ValueError(f"{protein_name} contains invalid amino acids: {bad}")

    src_ids = encode_aa(aa_seq)
    src = torch.tensor([src_ids], dtype=torch.long, device=device)
    gc_tensor = torch.tensor([[TARGET_GC]], dtype=torch.float32, device=device)
    src_pad_mask = src == aa_vocab["<PAD>"]


    with torch.no_grad():
        logits = model(src, gc_tensor, src_key_padding_mask=src_pad_mask)[0]
        mask = valid_codon_mask[src[0]]
        masked_logits = logits.masked_fill(~mask, -1e9)
        adjusted_logits = masked_logits + GC_DECODING_BONUS * gc_weights.view(1, -1) + TAI_DECODING_BONUS * tai_weights.view(1, -1)

    valid_positions = list(range(1, len(src_ids) - 1))


# Greedy candidate
    candidates = []
    greedy_codons = []
    for pos in valid_positions:
        idx = int(torch.argmax(adjusted_logits[pos]).item())
        greedy_codons.append(idx_to_codon[idx])
    candidates.append(("greedy", greedy_codons))

# Top-k random candidates: retain only the best TOP_K valid codons at each position, then sample according to their model probabilities. This adds diversity while avoiding uniform random selection among all synonymous codons.

    position_topk = []
    for pos in valid_positions:
        scores = adjusted_logits[pos]
        finite_idx = torch.where(scores > -1e8)[0]
        k = min(TOP_K, finite_idx.numel())
        _, local_idx = torch.topk(scores[finite_idx], k=k)
        top_indices = finite_idx[local_idx].detach().cpu().tolist()
        position_topk.append(top_indices)

    for cand_i in range(1, NUM_CANDIDATES_PER_CHAIN):
        codons = []
        for pos, choices in zip(valid_positions, position_topk):
            if len(choices) == 1:
                chosen_idx = choices[0]
            else:

                choice_tensor = torch.tensor(choices, dtype=torch.long, device=device)
                choice_scores = adjusted_logits[pos, choice_tensor] / SAMPLING_TEMPERATURE
                probs = torch.softmax(choice_scores, dim=0)
                sampled_local = torch.multinomial(probs, num_samples=1).item()
                chosen_idx = choices[sampled_local]
            codons.append(idx_to_codon[chosen_idx])
        candidates.append((f"topk_{cand_i}", codons))

    rows = []
    for cand_label, codons in candidates:
        metrics = sequence_metrics(codons, tai_weight_dict)
        score = selection_score(metrics)
        aa_from_codons = "".join(genetic_code[c] for c in codons)
        aa_identity = sum(a == b for a, b in zip(aa_seq, aa_from_codons)) / len(aa_seq)
        rows.append({
            "candidate_id": f"{protein_name}_{cand_label}",
            "protein_name": protein_name,
            "aa_length": len(aa_seq),
            "dna": metrics["dna"],
            "length_nt": metrics["length_nt"],
            "gc_content": metrics["gc"],
            "tai_like_score": metrics["tai"],
            "mfe_proxy": metrics["mfe_proxy"],
            "rnafold_mfe": metrics["rnafold_mfe"],
            "rnafold_mfe_per_nt": metrics["rnafold_mfe_per_nt"],
            "restriction_site_total": metrics["restriction_site_total"],
            "polyA_signal_total": metrics["polyA_signal_total"],
            "aa_identity": aa_identity,
            "selection_score": score,
            "restriction_hits": metrics["restriction_hits"],
            "polyA_hits": metrics["polyA_hits"],
        })

    rows.sort(key=lambda x: x["selection_score"])
    return rows


def write_outputs(all_rows):
    """Write ranked candidates to FASTA and TSV files"""
    with open(OUTPUT_FASTA, "w") as fasta:
        for row in all_rows:
            fasta.write(f">{row['candidate_id']} gc={row['gc_content']:.4f} tai={row['tai_like_score']:.4f} score={row['selection_score']:.4f}\n")
            seq = row["dna"]
            for i in range(0, len(seq), 80):
                fasta.write(seq[i:i+80] + "\n")

    cols = [
        "candidate_id", "protein_name", "aa_length", "length_nt", "gc_content",
        "tai_like_score", "mfe_proxy", "rnafold_mfe", "rnafold_mfe_per_nt",
        "restriction_site_total", "polyA_signal_total", "aa_identity",
        "selection_score", "restriction_hits", "polyA_hits", "dna",
    ]
    with open(OUTPUT_TSV, "w") as tsv:
        tsv.write("\t".join(cols) + "\n")
        for row in all_rows:
            vals = []
            for c in cols:
                v = row[c]
                if isinstance(v, float):
                    vals.append(f"{v:.6f}")
                else:
                    vals.append(str(v))
            tsv.write("\t".join(vals) + "\n")


def main():
    """Run the complete model-loading, decoding, scoring, and output workflow"""
# Windows CUDA first. Falls back to Apple MPS or CPU if CUDA is unavailable.
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print("Using device:", device)
        print("CUDA version:", torch.version.cuda)
        print("GPU:", torch.cuda.get_device_name(0))
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        device = torch.device("mps")
        print("Using device:", device)
    else:
        device = torch.device("cpu")
        print("Using device:", device)
    print("TARGET_GC:", TARGET_GC)
    print("GC_DECODING_BONUS:", GC_DECODING_BONUS)
    print("GC_FLOOR:", GC_FLOOR)
    print("GC_CEILING:", GC_CEILING)
    print("NUM_CANDIDATES_PER_CHAIN:", NUM_CANDIDATES_PER_CHAIN)
    print("SAMPLING_TEMPERATURE:", SAMPLING_TEMPERATURE)

    checkpoint = find_checkpoint()
    print("Loading checkpoint:", checkpoint)

    model = EncoderOnlyCodonTransformer(len(aa_vocab), len(codon_vocab)).to(device)
    state = torch.load(checkpoint, map_location=device)
    model.load_state_dict(state)
    model.eval()

    valid_codon_mask = build_valid_codon_mask(device)
    gc_weights = build_gc_weight_vector(device)
    tai_weights = build_tai_weights_from_trnascan(TRNA_SCAN_FILE, device)
    tai_weight_dict = {
        codon: float(tai_weights[idx].detach().cpu().item())
        for codon, idx in codon_vocab.items()
        if codon not in ["<PAD>", "<SOS>", "<EOS>"]
    }

    all_rows = []
    for name, aa_seq in TARGET_PROTEINS.items():
        print(f"Generating {name}: {len(aa_seq)} aa")
        rows = decode_candidates(model, name, aa_seq, device, valid_codon_mask, gc_weights, tai_weights, tai_weight_dict)
        all_rows.extend(rows)
        best = rows[0]
        print(
            f"  Best: {best['candidate_id']} | GC={best['gc_content']:.4f} | "
            f"tAI={best['tai_like_score']:.4f} | motifs={best['restriction_site_total'] + best['polyA_signal_total']} | "
            f"AA identity={best['aa_identity']:.4f} | score={best['selection_score']:.4f}"
        )

    write_outputs(all_rows)
    print(f"\nSaved FASTA: {OUTPUT_FASTA}")
    print(f"Saved TSV:   {OUTPUT_TSV}")


if __name__ == "__main__":
    main()
