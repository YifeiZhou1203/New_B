#!/usr/bin/env python3
"""
#########
Screen selected Adalimumab top candidates for plant miRNA seed matches.

DETAILED NOTES:
- This is a post-selection screen. It does not train the Transformer and does not
  generate new codon sequences.
- Each mature miRNA contributes a short seed pattern based on the first
  MIRNA_SEED_LEN nucleotides. That seed is reversely complemented into a DNA pattern
  and searched for exact matches in each candidate CDS.
- An exact 7-nt match is only a sequence flag. It is not sufficient evidence for
  real miRNA regulation in a plant cell; the original script intentionally treats
  this as weak evidence.
- The miRNA penalty is added to the previous selection score rather than replacing
  it. This keeps the earlier GC/tAI/MFE/motif criteria as the primary ranking layer.
- Candidates are still grouped as HC and LC for ranking/reporting.
- HARD_REMOVE_MIRNA_HITS is False by default so a seed hit is reported and mildly
  penalized rather than automatically discarded.


This script does NOT generate new candidates and does NOT require the Transformer model.
It starts from:
  - adalimumab_top5_candidates_NbT2T.tsv
  - adalimumab_top5_candidates_NbT2T.fasta
  - arabidopsis_mature_miRNA_major_families.fa  (or another miRNA FASTA)

Outputs:
  - adalimumab_top5_candidates_miRNA_screened_NbT2T.tsv
  - adalimumab_top5_candidates_miRNA_screened_NbT2T.fasta
  - adalimumab_top5_candidates_miRNA_hits_NbT2T.tsv
"""

from pathlib import Path
import ast
import pandas as pd


# Input / output files

INPUT_TSV = "antibody_top5_candidates_0527.tsv"
INPUT_FASTA = "antibody_top5_candidates_0527.fasta"
MIRNA_FILE = "arabidopsis_mature_miRNA_major_families.fa"




OUTPUT_TSV = "adalimumab_top5_candidates_miRNA_screened_NbT2T.tsv"
OUTPUT_FASTA = "adalimumab_top5_candidates_miRNA_screened_NbT2T.fasta"
OUTPUT_HITS_TSV = "adalimumab_top5_candidates_miRNA_hits_NbT2T.tsv"


# miRNA settings
# MIRNA_SEED_LEN controls the exact seed length searched. The penalty weight is small because short exact matches can occur by chance.

MIRNA_SEED_LEN = 7
MIRNA_SEED_PENALTY_WEIGHT = 1.0  # keep weak; 7-nt exact seed hits can be false positives

# if True, remove candidates with any seed hit from the FASTA output if alternatives exist.
HARD_REMOVE_MIRNA_HITS = False

# =========================
# Helpers
# =========================
# These helpers handle FASTA I/O, DNA normalization, reverse complementation,
# miRNA seed construction, hit counting, flexible ID handling, and numeric conversion.

def read_fasta(path):
    """Read FASTA records as (identifier, sequence) pairs"""
    records = []
    current_id = None
    seq_parts = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(seq_parts)))
                current_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    if current_id is not None:
        records.append((current_id, "".join(seq_parts)))

    return records


def write_fasta(records, path, width=80):
    """Write screened DNA records to FASTA with normalized sequences"""
    with open(path, "w") as out:
        for record_id, seq in records:
            out.write(f">{record_id}\n")
            seq = clean_dna(seq)
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def clean_dna(seq):
    """Normalize a nucleotide string to uppercase DNA using ATCG only"""
    seq = str(seq).upper().replace("U", "T")
    return "".join(b for b in seq if b in "ATCG")


def reverse_complement_dna(seq):
    """Return the reverse complement of a cleaned DNA sequence"""
    seq = clean_dna(seq)
    comp = str.maketrans("ATCG", "TAGC")
    return seq.translate(comp)[::-1]


def mirna_seed_to_target_dna(mirna_seq, seed_len=7):
    mirna_seq = clean_dna(mirna_seq)
    if len(mirna_seq) < seed_len:
        return None
    return reverse_complement_dna(mirna_seq[:seed_len])


def read_mirna_seeds(path, seed_len=7):
    """Read mature-miRNA FASTA records and create ID-to-seed mappings"""
    if not Path(path).exists():
        raise FileNotFoundError(
            f"Missing {path}. Put your miRNA FASTA in this folder or change MIRNA_FILE."
        )

    seeds = {}
    current_id = None
    seq_parts = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None and seq_parts:
                    seed = mirna_seed_to_target_dna("".join(seq_parts), seed_len)
                    if seed:
                        seeds[current_id] = seed
                current_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    if current_id is not None and seq_parts:
        seed = mirna_seed_to_target_dna("".join(seq_parts), seed_len)
        if seed:
            seeds[current_id] = seed

    return seeds


def count_seed_hits(dna, mirna_seeds):
    """Count exact occurrences of every miRNA seed in one candidate DNA sequence"""
    dna = clean_dna(dna)
    hits = {}
    total = 0
    for mirna_id, seed in mirna_seeds.items():
        count = dna.count(seed)
        if count > 0:
            hits[mirna_id] = count
            total += count
    return hits, total


def get_candidate_id_column(df):
    """Find an accepted candidate-ID column in the input DataFrame"""
    for col in ["candidate_id", "id", "name", "sequence_id"]:
        if col in df.columns:
            return col
    raise ValueError("Could not find candidate ID column. Expected 'candidate_id'.")


def get_dna_for_row(row, fasta_dict):
    """Retrieve a candidate's DNA from FASTA, with TSV DNA as fallback"""
    cid = str(row["candidate_id"])
    if cid in fasta_dict:
        return clean_dna(fasta_dict[cid])
    if "dna" in row and pd.notna(row["dna"]):
        return clean_dna(row["dna"])
    return ""


def safe_float_column(df, col, default=0.0):
    """Return a numeric DataFrame column with safe defaults for missing values"""
    if col not in df.columns:
        return pd.Series([default] * len(df), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


#######################
# Main fucntion


def main():
    """Run the complete miRNA seed-screening and ranking workflow.Candidates are loaded, miRNA seeds are generated, exact matches are counted,a weak penalty is added to the previous score, HC/LC rankings are produced,and screened TSV/FASTA plus a detailed hit table are written."""
    if not Path(INPUT_TSV).exists():
        raise FileNotFoundError(f"Missing {INPUT_TSV}")
    if not Path(INPUT_FASTA).exists():
        print(f"WARNING: {INPUT_FASTA} not found. Will try to use dna column from TSV only.")
        fasta_records = []
    else:
        fasta_records = read_fasta(INPUT_FASTA)

    fasta_dict = {rid: clean_dna(seq) for rid, seq in fasta_records}
    df = pd.read_csv(INPUT_TSV, sep="\t")

    id_col = get_candidate_id_column(df)
    if id_col != "candidate_id":
        df = df.rename(columns={id_col: "candidate_id"})

    mirna_seeds = read_mirna_seeds(MIRNA_FILE, MIRNA_SEED_LEN)
    print(f"Loaded candidates: {len(df)}")
    print(f"Loaded miRNA seed patterns: {len(mirna_seeds)}")
    print(f"MIRNA_SEED_LEN: {MIRNA_SEED_LEN}")
    print(f"MIRNA_SEED_PENALTY_WEIGHT: {MIRNA_SEED_PENALTY_WEIGHT}")

    hit_rows = []
    screened_dna = []
    total_hits = []
    unique_mirna_hits = []
    hit_detail_strings = []

# add one miRNA seed may occur more than once in the same CDS.
    for _, row in df.iterrows():
        cid = str(row["candidate_id"])
        dna = get_dna_for_row(row, fasta_dict)
        screened_dna.append(dna)

        hits, total = count_seed_hits(dna, mirna_seeds)
        total_hits.append(total)
        unique_mirna_hits.append(len(hits))
        hit_detail_strings.append(str(hits))

        for mirna_id, count in hits.items():
            hit_rows.append({
                "candidate_id": cid,
                "mirna_id": mirna_id,
                "target_seed_dna": mirna_seeds[mirna_id],
                "hit_count": count,
            })

    df["dna_for_screening"] = screened_dna
    df["mirna_seed_total"] = total_hits
    df["mirna_seed_unique_count"] = unique_mirna_hits
    df["mirna_seed_hits"] = hit_detail_strings
    df["mirna_seed_penalty"] = df["mirna_seed_total"] * MIRNA_SEED_PENALTY_WEIGHT


# Preserve the existing final score if available; otherwise use selection_score if present.
    if "final_selection_score" in df.columns:
        base_score_col = "final_selection_score"
    elif "selection_score" in df.columns:
        base_score_col = "selection_score"
    else:
        base_score_col = None

    if base_score_col:
        df["base_score_before_miRNA"] = safe_float_column(df, base_score_col)
        df["miRNA_adjusted_score"] = df["base_score_before_miRNA"] + df["mirna_seed_penalty"]
    else:
        print("WARNING: No final_selection_score/selection_score found. Ranking by miRNA hits only.")
        df["base_score_before_miRNA"] = 0.0
        df["miRNA_adjusted_score"] = df["mirna_seed_penalty"]

# Rank HC and LC separately if candidate IDs contain HC/LC.
    df["chain_type"] = "Other"
    df.loc[df["candidate_id"].str.contains("HC", case=False, na=False), "chain_type"] = "HC"
    df.loc[df["candidate_id"].str.contains("LC", case=False, na=False), "chain_type"] = "LC"

    df_ranked = (
        df.sort_values(["chain_type", "miRNA_adjusted_score", "mirna_seed_total"], ascending=[True, True, True])
          .reset_index(drop=True)
    )


    if HARD_REMOVE_MIRNA_HITS:
        keep_parts = []
        for chain, group in df_ranked.groupby("chain_type", sort=False):
            clean_group = group[group["mirna_seed_total"] == 0]
            keep_parts.append(clean_group if len(clean_group) > 0 else group)
        df_output = pd.concat(keep_parts, ignore_index=True)
    else:
        df_output = df_ranked

    df_output.to_csv(OUTPUT_TSV, sep="\t", index=False)

# Write FASTA in ranked order.
    out_records = []
    for _, row in df_output.iterrows():
        cid = str(row["candidate_id"])
        dna = row["dna_for_screening"]
        if dna:
            out_records.append((cid, dna))
    write_fasta(out_records, OUTPUT_FASTA)

    hits_df = pd.DataFrame(hit_rows)
    if len(hits_df) == 0:
        hits_df = pd.DataFrame(columns=["candidate_id", "mirna_id", "target_seed_dna", "hit_count"])
    hits_df.to_csv(OUTPUT_HITS_TSV, sep="\t", index=False)

    print("\nSummary by chain:")
    summary_cols = ["candidate_id", "mirna_seed_total", "mirna_seed_unique_count", "base_score_before_miRNA", "miRNA_adjusted_score"]
    for chain in ["HC", "LC", "Other"]:
        sub = df_output[df_output["chain_type"] == chain]
        if len(sub) > 0:
            print(f"\n{chain} ranked candidates:")
            print(sub[summary_cols].to_string(index=False))

    print(f"\nSaved screened table: {OUTPUT_TSV}")
    print(f"Saved screened FASTA: {OUTPUT_FASTA}")
    print(f"Saved miRNA hit details: {OUTPUT_HITS_TSV}")
    print("\nNote: exact seed hits are weak evidence only; validate important candidates experimentally or with a stronger plant miRNA-target tool.")


if __name__ == "__main__":
    main()
