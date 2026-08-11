#!/usr/bin/env python3
"""
Run RNAfold on generated Adalimumab candidate CDS sequences and select the best 5. Windows-compatible if RNAfold is available on PATH.

Expected input from antibody_generator.py:
    adalimumab_generated_candidates_NbT2T.fasta
    adalimumab_generated_candidates_NbT2T.tsv

Outputs:
    adalimumab_rnafold_results_NbT2T.tsv
    adalimumab_top5_candidates_NbT2T.tsv
    adalimumab_top5_candidates_NbT2T.fasta

Install RNAfold if needed:
    conda install -c conda-forge -c bioconda viennarna
"""

from pathlib import Path
import re
import shutil
import subprocess
import pandas as pd

#
# Input / output files

INPUT_FASTA = "antibody_candidates.fasta"
INPUT_SUMMARY = "antibody_candidates.tsv"

RNAFOLD_OUTPUT = "antibody_rnafold_results.tsv"
TOP5_TSV = "antibody_top5_candidates.tsv"
TOP5_FASTA = "antibody_top5_candidates.fasta"


# Selection targets


TARGET_GC = 0.43
TARGET_TAI = 0.3765
TARGET_MFE_PER_NT = -0.30
GC_CEILING = 0.48
GC_FLOOR = 0.41
HARD_GC_FILTER = True

# Scoring weights. Lower final_selection_score is better.
GC_WEIGHT = 4.0
GC_CEILING_WEIGHT = 20.0  # strong penalty only when GC exceeds 48%
TAI_WEIGHT = 1.0
MFE_WEIGHT = 1.0
RESTRICTION_PENALTY_WEIGHT = 10.0
POLYA_PENALTY_WEIGHT = 10.0

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


# HELPER
# These functions isolate file parsing, nucleotide normalization, motif counting,RNAfold execution, and flexible column-name handling from the ranking workflow.

def read_fasta(path):
    """Read FASTA records as (identifier, sequence) tuples.

    Multi-line sequence records are concatenated. Only the first whitespace-
    separated token of each FASTA header is used as the identifier.
    """
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
    """Write sequence records as wrapped FASTA entries"""
    with open(path, "w") as out:
        for record_id, seq in records:
            out.write(f">{record_id}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def clean_dna(seq):
    """Normalize an input sequence to uppercase DNA containing ATCG only"""
    seq = seq.upper().replace("U", "T")
    return "".join(b for b in seq if b in "ATCG")


def gc_content(dna):
    """Return the fraction of nucleotides that are G or C.Empty input is assigned 0.0 to avoid division by zero"""
    if not dna:
        return 0.0
    return sum(1 for b in dna if b in "GC") / len(dna)


def count_motifs(dna, motif_dict):
    """Count non-overlapping occurrences of each motif in a DNA sequence"""
    total = 0
    details = {}
    for name, motif in motif_dict.items():
        count = dna.count(motif)
        details[name] = count
        total += count
    return total, details


def run_rnafold(rna):
    """Run ViennaRNA RNAfold and parse structure text and MFE"""
    result = subprocess.run(
        ["RNAfold", "--noPS"],
        input=rna + "\n",
        text=True,
        capture_output=True,
        check=False,
    )

    if result.returncode != 0:
        return None, None, result.stderr.strip()

    lines = result.stdout.strip().splitlines()
    if len(lines) < 2:
        return None, None, "RNAfold output missing structure line"

    structure_line = lines[1]
    match = re.search(r"\(([-+]?[0-9]*\.?[0-9]+)\)", structure_line)
    if match is None:
        return structure_line, None, "Could not parse MFE"

    mfe = float(match.group(1))
    return structure_line, mfe, None


def find_column(df, candidates):
    """Return the first matching column name from a list of accepted aliases"""
    for c in candidates:
        if c in df.columns:
            return c
    return None


# Main ###################

def main():
    """Run RNAfold analysis, scoring, GC filtering, and HC/LC selection."""
    if shutil.which("RNAfold") is None:
        raise RuntimeError(
            "RNAfold not found on PATH. Install ViennaRNA/RNAfold first, then reopen terminal. "
            "Conda option: conda install -c conda-forge -c bioconda viennarna"
        )

    fasta_path = Path(INPUT_FASTA)
    if not fasta_path.exists():
        raise FileNotFoundError(f"Missing {INPUT_FASTA}. Run generate_adalimumab_current_model.py first.")

    records = read_fasta(INPUT_FASTA)
    if not records:
        raise ValueError(f"No FASTA records found in {INPUT_FASTA}")

    fasta_dict = {}
    rnafold_rows = []

    print(f"Running RNAfold on {len(records)} candidates...")
    print("TARGET_GC:", TARGET_GC)
    print("GC_CEILING:", GC_CEILING)
    print("GC_FLOOR:", GC_FLOOR)
    print("HARD_GC_FILTER:", HARD_GC_FILTER)

# Fold candidates independently.
    for candidate_id, seq in records:
        dna = clean_dna(seq)
        fasta_dict[candidate_id] = dna

        if not dna:
            rnafold_rows.append({
                "candidate_id": candidate_id,
                "length_nt": 0,
                "gc_content_rnafold": None,
                "mfe": None,
                "mfe_per_nt": None,
                "structure_or_error": "invalid_sequence",
            })
            continue

        rna = dna.replace("T", "U")
        structure, mfe, error = run_rnafold(rna)

        rnafold_rows.append({
            "candidate_id": candidate_id,
            "length_nt": len(dna),
            "gc_content_rnafold": gc_content(dna),
            "mfe": mfe,
            "mfe_per_nt": (mfe / len(dna)) if mfe is not None else None,
            "structure_or_error": error if error else structure,
        })

    rnafold_df = pd.DataFrame(rnafold_rows)
    rnafold_df.to_csv(RNAFOLD_OUTPUT, sep="\t", index=False)
    print(f"Saved RNAfold results: {RNAFOLD_OUTPUT}")

# Merge optional summary file if present.
    if Path(INPUT_SUMMARY).exists():
        summary = pd.read_csv(INPUT_SUMMARY, sep="\t")
        df = summary.merge(rnafold_df, on="candidate_id", how="right")
    else:
        print(f"WARNING: {INPUT_SUMMARY} not found. Ranking will use FASTA + RNAfold only.")
        df = rnafold_df.copy()

# Choose usable GC column.
    gc_col = find_column(df, ["gc_content", "generated_gc", "gc", "gc_content_rnafold"])
    if gc_col is None:
        df["gc_content_for_ranking"] = df["candidate_id"].map(lambda cid: gc_content(fasta_dict.get(cid, "")))
        gc_col = "gc_content_for_ranking"

# Choose usable tAI column if present.
    tai_col = find_column(df, ["tai_like_score", "generated_tai", "tai", "tAI", "tai_score"])
    if tai_col is None:
        df["tai_diff"] = 0.0
        print("WARNING: No tAI column found. tAI will not affect ranking.")
    else:
        df["tai_diff"] = (df[tai_col] - TARGET_TAI).abs()

# Count motifs directly from FASTA so ranking is robust.
    restriction_totals = []
    polya_totals = []
    for cid in df["candidate_id"]:
        dna = fasta_dict.get(cid, "")
        restriction_total, _ = count_motifs(dna, FORBIDDEN_RESTRICTION_SITES)
        polya_total, _ = count_motifs(dna, POLYA_SIGNALS)
        restriction_totals.append(restriction_total)
        polya_totals.append(polya_total)

    df["restriction_site_total"] = restriction_totals
    df["polyA_signal_total"] = polya_totals

    df["gc_diff"] = (df[gc_col] - TARGET_GC).abs()
# Asymmetric ceiling penalty: candidates above 46% GC are pushed down the ranking.
# Candidates below 46% receive no ceiling penalty and are ranked by closeness to TARGET_GC.
    df["gc_ceiling_penalty"] = (df[gc_col] - GC_CEILING).clip(lower=0)
    df["gc_floor_penalty"] = (GC_FLOOR - df[gc_col]).clip(lower=0)
    df["mfe_per_nt_diff"] = (df["mfe_per_nt"] - TARGET_MFE_PER_NT).abs()

# Invalid RNAfold candidates should rank last........
    df["mfe_per_nt_diff"] = df["mfe_per_nt_diff"].fillna(999.0)

    df["final_selection_score"] = (
        GC_WEIGHT * df["gc_diff"]
        + GC_CEILING_WEIGHT * df["gc_ceiling_penalty"]
        + GC_CEILING_WEIGHT * df["gc_floor_penalty"]
        + TAI_WEIGHT * df["tai_diff"]
        + MFE_WEIGHT * df["mfe_per_nt_diff"]
        + RESTRICTION_PENALTY_WEIGHT * df["restriction_site_total"]
        + POLYA_PENALTY_WEIGHT * df["polyA_signal_total"]
    )
#############################
# Apply the hard GC eligibility rule before final ranking when requested.
# This prevents a candidate outside the allowed GC interval from winning because it has a favorable RNAfold/MFE value.
# Hard biological GC filter first, then ranking.

    df_all_ranked = df.sort_values("final_selection_score", ascending=True).reset_index(drop=True)

    if HARD_GC_FILTER:
        df_valid_gc = df[(df[gc_col] >= GC_FLOOR) & (df[gc_col] <= GC_CEILING)].copy()
        print(f"GC hard filter: kept {len(df_valid_gc)} / {len(df)} candidates within {GC_FLOOR:.2f}-{GC_CEILING:.2f}")
        if len(df_valid_gc) == 0:
            print("WARNING: No candidates passed hard GC filter. Falling back to soft-ranked candidates.")
            df_for_ranking = df_all_ranked.copy()
        else:
            df_for_ranking = df_valid_gc.sort_values("final_selection_score", ascending=True).reset_index(drop=True)
    else:
        df_for_ranking = df_all_ranked.copy()

# Rank HC and LC separately.
    hc_ranked = df_for_ranking[df_for_ranking["candidate_id"].str.contains("HC", case=False, na=False)].copy()
    lc_ranked = df_for_ranking[df_for_ranking["candidate_id"].str.contains("LC", case=False, na=False)].copy()

    if len(hc_ranked) == 0:
        print("WARNING: No HC candidates passed hard GC filter. Falling back to best soft-ranked HC candidates.")
        hc_ranked = df_all_ranked[df_all_ranked["candidate_id"].str.contains("HC", case=False, na=False)].copy()
    if len(lc_ranked) == 0:
        print("WARNING: No LC candidates passed hard GC filter. Falling back to best soft-ranked LC candidates.")
        lc_ranked = df_all_ranked[df_all_ranked["candidate_id"].str.contains("LC", case=False, na=False)].copy()

    hc_top5 = hc_ranked.head(5).copy()
    lc_top5 = lc_ranked.head(5).copy()
    top5 = pd.concat([hc_top5, lc_top5], ignore_index=True)
    top5.to_csv(TOP5_TSV, sep="\t", index=False)

    top5_records = [(cid, fasta_dict[cid]) for cid in top5["candidate_id"] if cid in fasta_dict]
    write_fasta(top5_records, TOP5_FASTA)

    display_cols = [
        "candidate_id",
        gc_col,
        "gc_diff",
        "gc_ceiling_penalty",
        "gc_floor_penalty",
        tai_col if tai_col else None,
        "tai_diff",
        "mfe",
        "mfe_per_nt",
        "mfe_per_nt_diff",
        "restriction_site_total",
        "polyA_signal_total",
        "final_selection_score",
    ]
    display_cols = [c for c in display_cols if c is not None and c in top5.columns]

    print("\nTop HC candidates:")
    print(hc_top5[display_cols].to_string(index=False))
    print("\nTop LC candidates:")
    print(lc_top5[display_cols].to_string(index=False))
    print(f"\nSaved selected candidates table: {TOP5_TSV}")
    print(f"Saved selected candidates FASTA: {TOP5_FASTA}")

if __name__ == "__main__":
    main()
