from pathlib import Path
import math
from itertools import permutations

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path("/mnt/d/hanjia/2.10/dnds")
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "output"

REFERENCE_FASTA = DATA_DIR / "reference.fasta"
ANNOTATION_GFF = DATA_DIR / "genome_annotation.gff3"
GLOBAL_FASTA = DATA_DIR / "global.fasta"
SHENZHEN_FASTA = DATA_DIR / "197_shenzhen_rename.fasta"

MIN_GENOME_LEN_FRACTION = 0.95
MAX_GENOME_AMBIGUOUS_FRACTION = 0.05
MIN_VALID_CODON_FRACTION = 0.8
PSEUDOCOUNT = 0.5


CODON_TABLE = {
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

VALID_BASES = set("ACGT")
RC_TABLE = str.maketrans("ACGTN-", "TGCAN-")


def parse_fasta(path):
    name = None
    seq = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if name is not None:
            yield name, "".join(seq).upper()


def reverse_complement(seq):
    return seq.translate(RC_TABLE)[::-1]


def parse_attributes(attr_text):
    attributes = {}
    for item in attr_text.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key] = value
    return attributes


def parse_cds_features(gff_path):
    cds_features = []
    seen_names = set()
    with open(gff_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parse_attributes(parts[8])
            gene_name = attributes.get("Name") or attributes.get("gene")
            if not gene_name or gene_name in seen_names:
                continue
            cds_features.append(
                {
                    "gene": gene_name,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "product": attributes.get("product", ""),
                }
            )
            seen_names.add(gene_name)
    cds_features.sort(key=lambda item: item["start"])
    return cds_features


def translate_codon(codon):
    return CODON_TABLE.get(codon)


def synonymous_sites_for_codon(codon):
    aa = translate_codon(codon)
    if aa is None or aa == "*":
        return (0.0, 0.0)

    synonymous = 0.0
    for idx, base in enumerate(codon):
        per_position_syn = 0
        for alt in "ACGT":
            if alt == base:
                continue
            mutated = codon[:idx] + alt + codon[idx + 1 :]
            mutated_aa = translate_codon(mutated)
            if mutated_aa == aa:
                per_position_syn += 1
        synonymous += per_position_syn / 3.0
    nonsynonymous = 3.0 - synonymous
    return synonymous, nonsynonymous


def classify_codon_change(ref_codon, alt_codon):
    diffs = [idx for idx, (r, a) in enumerate(zip(ref_codon, alt_codon)) if r != a]
    if not diffs:
        return 0.0, 0.0

    if any(base not in VALID_BASES for base in ref_codon + alt_codon):
        return None

    if len(diffs) == 1:
        ref_aa = translate_codon(ref_codon)
        alt_aa = translate_codon(alt_codon)
        if ref_aa is None or alt_aa is None or "*" in {ref_aa, alt_aa}:
            return None
        return (1.0, 0.0) if ref_aa == alt_aa else (0.0, 1.0)

    syn_total = 0.0
    nonsyn_total = 0.0
    valid_paths = 0
    for order in permutations(diffs):
        current = ref_codon
        syn_count = 0.0
        nonsyn_count = 0.0
        valid = True
        for idx in order:
            next_codon = current[:idx] + alt_codon[idx] + current[idx + 1 :]
            current_aa = translate_codon(current)
            next_aa = translate_codon(next_codon)
            if current_aa is None or next_aa is None or "*" in {current_aa, next_aa}:
                valid = False
                break
            if current_aa == next_aa:
                syn_count += 1.0
            else:
                nonsyn_count += 1.0
            current = next_codon
        if valid:
            syn_total += syn_count
            nonsyn_total += nonsyn_count
            valid_paths += 1

    if valid_paths == 0:
        return None
    return syn_total / valid_paths, nonsyn_total / valid_paths


def extract_oriented_gene(seq, feature):
    start = feature["start"] - 1
    end = feature["end"]
    if end > len(seq):
        return None
    gene_seq = seq[start:end]
    if feature["strand"] == "-":
        gene_seq = reverse_complement(gene_seq)
    return gene_seq


def gene_metrics_from_pair(ref_gene, sample_gene):
    codon_count = len(ref_gene) // 3
    syn_sites = 0.0
    nonsyn_sites = 0.0
    syn_changes = 0.0
    nonsyn_changes = 0.0
    valid_codons = 0

    for idx in range(0, codon_count * 3, 3):
        ref_codon = ref_gene[idx : idx + 3]
        alt_codon = sample_gene[idx : idx + 3]
        if any(base not in VALID_BASES for base in ref_codon + alt_codon):
            continue
        if "*" in {translate_codon(ref_codon), translate_codon(alt_codon)}:
            continue
        valid_codons += 1
        s_sites, n_sites = synonymous_sites_for_codon(ref_codon)
        syn_sites += s_sites
        nonsyn_sites += n_sites
        change = classify_codon_change(ref_codon, alt_codon)
        if change is None:
            continue
        syn_changes += change[0]
        nonsyn_changes += change[1]

    return {
        "valid_codons": valid_codons,
        "total_codons": codon_count,
        "syn_sites": syn_sites,
        "nonsyn_sites": nonsyn_sites,
        "syn_changes": syn_changes,
        "nonsyn_changes": nonsyn_changes,
    }


def genome_quality_ok(seq, reference_len):
    seq = seq.upper()
    good_len = len(seq) >= reference_len * MIN_GENOME_LEN_FRACTION
    ambiguous = sum(base not in VALID_BASES for base in seq)
    ambiguous_fraction = ambiguous / len(seq) if seq else 1.0
    return good_len and ambiguous_fraction <= MAX_GENOME_AMBIGUOUS_FRACTION


def summarize_group(group_name, fasta_path, reference_seq, features):
    aggregates = {
        feature["gene"]: {
            "group": group_name,
            "gene": feature["gene"],
            "start": feature["start"],
            "end": feature["end"],
            "strand": feature["strand"],
            "product": feature["product"],
            "gene_length_nt": feature["end"] - feature["start"] + 1,
            "gene_length_aa": (feature["end"] - feature["start"] + 1) // 3,
            "samples_used": 0,
            "valid_codons_total": 0,
            "syn_sites": 0.0,
            "nonsyn_sites": 0.0,
            "syn_changes": 0.0,
            "nonsyn_changes": 0.0,
        }
        for feature in features
    }

    reference_genes = {
        feature["gene"]: extract_oriented_gene(reference_seq, feature)
        for feature in features
    }

    total_records = 0
    passed_records = 0
    reference_len = len(reference_seq)

    for _, sample_seq in parse_fasta(fasta_path):
        total_records += 1
        if not genome_quality_ok(sample_seq, reference_len):
            continue
        passed_records += 1
        for feature in features:
            ref_gene = reference_genes[feature["gene"]]
            sample_gene = extract_oriented_gene(sample_seq, feature)
            if sample_gene is None or len(sample_gene) != len(ref_gene):
                continue
            metrics = gene_metrics_from_pair(ref_gene, sample_gene)
            if metrics["total_codons"] == 0:
                continue
            valid_fraction = metrics["valid_codons"] / metrics["total_codons"]
            if valid_fraction < MIN_VALID_CODON_FRACTION:
                continue

            gene_stats = aggregates[feature["gene"]]
            gene_stats["samples_used"] += 1
            gene_stats["valid_codons_total"] += metrics["valid_codons"]
            gene_stats["syn_sites"] += metrics["syn_sites"]
            gene_stats["nonsyn_sites"] += metrics["nonsyn_sites"]
            gene_stats["syn_changes"] += metrics["syn_changes"]
            gene_stats["nonsyn_changes"] += metrics["nonsyn_changes"]

    rows = []
    for feature in features:
        row = aggregates[feature["gene"]]
        syn_sites = row["syn_sites"]
        nonsyn_sites = row["nonsyn_sites"]
        syn_changes = row["syn_changes"]
        nonsyn_changes = row["nonsyn_changes"]

        dn = nonsyn_changes / nonsyn_sites if nonsyn_sites > 0 else math.nan
        ds = syn_changes / syn_sites if syn_sites > 0 else math.nan
        raw_dnds = dn / ds if ds and not math.isnan(ds) else math.nan

        dn_adj = (nonsyn_changes + PSEUDOCOUNT) / (nonsyn_sites + PSEUDOCOUNT)
        ds_adj = (syn_changes + PSEUDOCOUNT) / (syn_sites + PSEUDOCOUNT)
        dnds_adj = dn_adj / ds_adj if ds_adj > 0 else math.nan

        row["dn"] = dn
        row["ds"] = ds
        row["dnds_raw"] = raw_dnds
        row["dn_adj"] = dn_adj
        row["ds_adj"] = ds_adj
        row["dnds_adj"] = dnds_adj
        row["records_total"] = total_records
        row["records_passed_qc"] = passed_records
        rows.append(row)

    summary = {
        "group": group_name,
        "records_total": total_records,
        "records_passed_qc": passed_records,
        "qc_pass_rate": passed_records / total_records if total_records else 0.0,
    }
    return rows, summary


def write_results(rows, summaries, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    results_path = output_dir / "dnds_by_gene.tsv"
    pd.DataFrame(rows).to_csv(results_path, sep="\t", index=False)

    summary_path = output_dir / "run_summary.tsv"
    pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False)
    return results_path, summary_path


def make_plot(results_df, output_dir):
    ordered_genes = (
        results_df[["gene", "start"]]
        .drop_duplicates()
        .sort_values("start")
        .reset_index(drop=True)
    )
    gene_positions = {gene: idx for idx, gene in enumerate(ordered_genes["gene"])}
    results_df = results_df.copy()
    results_df["x"] = results_df["gene"].map(gene_positions)

    fig, axes = plt.subplots(2, 1, figsize=(34, 12), sharex=True, constrained_layout=True)
    panel_specs = [
        ("Global", "#3b82b8", axes[0], "A"),
        ("Shenzhen", "#d96b45", axes[1], "B"),
    ]

    for group_name, color, ax, panel_label in panel_specs:
        subset = results_df[results_df["group"] == group_name].sort_values("start")
        ax.scatter(
            subset["x"],
            subset["dnds_adj"],
            s=38,
            color=color,
            edgecolor="black",
            linewidth=0.3,
            alpha=0.9,
        )
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
        ax.set_ylabel("dN/dS", fontsize=13)
        ax.set_title(f"{panel_label}  {group_name}", loc="left", fontsize=15, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle=":", alpha=0.35)

        ymax = subset["dnds_adj"].replace([math.inf, -math.inf], math.nan).max()
        ymax = 1.2 if pd.isna(ymax) else max(1.2, min(ymax * 1.15, 8.0))
        ax.set_ylim(0, ymax)

    axes[1].set_xlabel("MPXV genes / protein regions", fontsize=13)
    tick_step = 5
    tick_positions = [idx for idx in range(len(ordered_genes)) if idx % tick_step == 0]
    tick_labels = [ordered_genes.loc[idx, "gene"] for idx in tick_positions]
    axes[1].set_xticks(tick_positions)
    axes[1].set_xticklabels(tick_labels, rotation=90, fontsize=7)

    fig.suptitle(
        "Approximate gene-level dN/dS comparison between global and Shenzhen MPXV genomes",
        fontsize=16,
        y=1.02,
    )

    svg_path = output_dir / "dnds_global_vs_shenzhen.svg"
    pdf_path = output_dir / "dnds_global_vs_shenzhen.pdf"
    plt.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return svg_path, pdf_path


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    reference_records = list(parse_fasta(REFERENCE_FASTA))
    if len(reference_records) != 1:
        raise ValueError("reference.fasta should contain exactly one sequence.")
    reference_name, reference_seq = reference_records[0]

    features = parse_cds_features(ANNOTATION_GFF)
    global_rows, global_summary = summarize_group("Global", GLOBAL_FASTA, reference_seq, features)
    shenzhen_rows, shenzhen_summary = summarize_group("Shenzhen", SHENZHEN_FASTA, reference_seq, features)

    results_df = pd.DataFrame(global_rows + shenzhen_rows)
    results_path, summary_path = write_results(results_df.to_dict("records"), [global_summary, shenzhen_summary], OUTPUT_DIR)
    svg_path, pdf_path = make_plot(results_df, OUTPUT_DIR)

    print("Reference:", reference_name)
    print("CDS features:", len(features))
    print("Results table:", results_path)
    print("Run summary:", summary_path)
    print("Figure SVG:", svg_path)
    print("Figure PDF:", pdf_path)
    print()
    print("QC summary:")
    for summary in [global_summary, shenzhen_summary]:
        print(
            f"  {summary['group']}: {summary['records_passed_qc']}/{summary['records_total']} passed "
            f"({summary['qc_pass_rate']:.1%})"
        )


if __name__ == "__main__":
    main()
