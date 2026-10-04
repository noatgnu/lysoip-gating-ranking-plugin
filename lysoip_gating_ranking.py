#!/usr/bin/env python3
"""Gates proteins on hard statistics, then z-scored weighted composite-ranks each gated pool."""

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

SIGNED_METRICS = {"fold_change", "organelle_specificity"}


def read_metric_column(path, metric_name, value_column):
    """Returns {protein: value} and {protein: gene} from a scorer output file."""
    values = {}
    genes = {}
    if not path:
        return values, genes
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            values[row["protein"]] = float(row[value_column])
            genes[row["protein"]] = row.get("gene", "")
    return values, genes


def main():
    parser = argparse.ArgumentParser(description="Gate proteins on hard statistics, then rank within each pool")
    parser.add_argument("--differential_expression_file", required=True)
    parser.add_argument("--reproducibility_file", default=None)
    parser.add_argument("--organelle_specificity_file", default=None)
    parser.add_argument("--go_evidence_file", default=None)
    parser.add_argument("--compass_file", default=None)
    parser.add_argument("--crapome_file", default=None)
    parser.add_argument("--fold_change_min_abs_log2", type=float, default=0.585)
    parser.add_argument("--significance_alpha", type=float, default=0.05)
    parser.add_argument("--weight_fold_change", type=float, default=1.0)
    parser.add_argument("--weight_reproducibility", type=float, default=1.0)
    parser.add_argument("--weight_organelle_specificity", type=float, default=1.0)
    parser.add_argument("--weight_go_evidence", type=float, default=1.0)
    parser.add_argument("--weight_compass", type=float, default=1.0)
    parser.add_argument("--weight_crapome", type=float, default=1.0)
    parser.add_argument("--output_folder", required=True)
    args = parser.parse_args()

    output_folder = Path(args.output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # @step: Loading differential expression
    fold_change, gene_lookup = read_metric_column(args.differential_expression_file, "fold_change", "log2_fold_change")
    significance, _ = read_metric_column(args.differential_expression_file, "significance", "p_adjusted")

    # @step: Loading optional scorer outputs
    reproducibility, repro_genes = read_metric_column(args.reproducibility_file, "reproducibility", "reproducibility")
    organelle_specificity, organelle_genes = read_metric_column(args.organelle_specificity_file, "organelle_specificity", "organelle_specificity")
    go_evidence, go_genes = read_metric_column(args.go_evidence_file, "go_evidence", "go_evidence")
    compass, compass_genes = read_metric_column(args.compass_file, "compass", "compass")
    crapome, crapome_genes = read_metric_column(args.crapome_file, "crapome", "crapome")
    gene_lookup.update(repro_genes)
    gene_lookup.update(organelle_genes)
    gene_lookup.update(go_genes)
    gene_lookup.update(compass_genes)
    gene_lookup.update(crapome_genes)

    metric_columns = {
        "fold_change": fold_change,
        "reproducibility": reproducibility,
        "organelle_specificity": organelle_specificity,
        "go_evidence": go_evidence,
        "compass": compass,
        "crapome": crapome,
    }
    weights = {
        "fold_change": args.weight_fold_change,
        "reproducibility": args.weight_reproducibility,
        "organelle_specificity": args.weight_organelle_specificity,
        "compass": args.weight_compass,
        "crapome": args.weight_crapome,
        "go_evidence": args.weight_go_evidence,
    }

    matrix = pd.DataFrame(metric_columns)
    matrix = matrix.dropna(how="all")

    # @step: Gating proteins on significance and fold-change magnitude
    positive_ids, negative_ids = set(), set()
    for protein in fold_change:
        fc = fold_change.get(protein)
        sig = significance.get(protein)
        if fc is None or sig is None:
            continue
        if sig > args.significance_alpha:
            continue
        if abs(fc) < args.fold_change_min_abs_log2:
            continue
        (positive_ids if fc >= 0 else negative_ids).add(protein)

    results = {}
    for protein in matrix.index:
        if protein not in positive_ids and protein not in negative_ids:
            results[protein] = ("Nonsignificant", None)

    # @step: Ranking within each gated pool
    def score_pool(pool_ids, classification):
        if not pool_ids:
            return
        pool = matrix.loc[matrix.index.isin(pool_ids)].copy()
        for signed_column in SIGNED_METRICS:
            if signed_column in pool.columns:
                pool[signed_column] = pool[signed_column].abs()

        present_weights = pd.Series({metric: weight for metric, weight in weights.items() if metric in pool.columns})
        features = pool[present_weights.index]
        normalized = (features - features.mean()) / features.std(ddof=0).replace(0, 1)
        normalized = normalized.fillna(0.0)

        composite = (normalized * present_weights).sum(axis=1) / present_weights.sum()
        for protein, value in composite.items():
            results[protein] = (classification, float(value))

    score_pool(positive_ids, "Positive")
    score_pool(negative_ids, "Negative")

    # @step: Writing verdicts
    rows = [
        {"protein": protein, "gene": gene_lookup.get(protein, ""), "classification": classification,
         "composite_value": "" if composite_value is None else composite_value}
        for protein, (classification, composite_value) in results.items()
    ]
    with open(output_folder / "verdict.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["protein", "gene", "classification", "composite_value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} verdicts ({len(positive_ids)} positive, {len(negative_ids)} negative)", file=sys.stderr)
    print("Gating and ranking complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
