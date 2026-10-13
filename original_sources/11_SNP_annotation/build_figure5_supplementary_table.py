import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
WORKDIR = Path(__file__).resolve().parent

SUMMARY = WORKDIR / "v4_full_workflow" / "figure5a_v4_snp_summary.csv"
CANDIDATES = WORKDIR / "v4_full_workflow" / "figure5a_v4_candidate_genes_by_window.csv"
GFF = WORKDIR / "Data_for_Menghui" / "LD_analysis" / "Morex_gff_file_all_chromosomes.txt"
OUT = WORKDIR / "Supplementary_Table_Figure5A_10_SNP_annotations.csv"


def parse_attrs(attr_text: str) -> dict:
    attrs = {}
    for item in attr_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def choose_primary_candidate(summary_row: dict) -> tuple[str, str]:
    if int(summary_row["overlap_gene_n"]) > 0:
        return summary_row["nearest_gene_id"], "overlapping_gene"
    if summary_row["best_50kb_gene"]:
        return summary_row["best_50kb_gene"], "best_gene_within_50kb"
    if summary_row["best_100kb_gene"]:
        return summary_row["best_100kb_gene"], "best_gene_within_100kb"
    return summary_row["best_250kb_gene"], "best_gene_within_250kb"


def main() -> None:
    summary_rows = [
        row
        for row in csv.DictReader(SUMMARY.open(newline="", encoding="utf-8"))
        if row["selection_reason"] == "Figure5A_positive_network"
    ]

    candidate_rows = list(csv.DictReader(CANDIDATES.open(newline="", encoding="utf-8")))
    overlap_gene_info = {}
    window_counts = {}

    for row in candidate_rows:
        snp_id = row["snp_id"]
        window_counts.setdefault(snp_id, {})[row["window_bp"]] = window_counts.setdefault(snp_id, {}).get(
            row["window_bp"], 0
        ) + 1
        if row["distance_to_snp_bp"] == "0" and snp_id not in overlap_gene_info:
            overlap_gene_info[snp_id] = row

    gene_meta = {}
    with GFF.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            attrs = parse_attrs(parts[8])
            if parts[2] == "gene":
                gene_id = attrs.get("ID", "")
                if gene_id:
                    gene_meta.setdefault(gene_id, {}).update(
                        {
                            "confidence": attrs.get("confidence", ""),
                            "biotype": attrs.get("biotype", ""),
                            "gene_name": attrs.get("Name", gene_id),
                        }
                    )
            elif parts[2] == "mRNA":
                parent = attrs.get("Parent", "")
                if parent:
                    gene_meta.setdefault(parent, {}).update(
                        {
                            "description": attrs.get("description", ""),
                            "ontology": attrs.get("ontology", ""),
                            "pfam": attrs.get("pfam", ""),
                        }
                    )

    out_rows = []
    for row in summary_rows:
        snp_id = row["snp_id"]
        primary_gene, primary_basis = choose_primary_candidate(row)
        primary_meta = gene_meta.get(primary_gene, {})
        overlap_row = overlap_gene_info.get(snp_id, {})
        overlap_gene = overlap_row.get("gene_id", "")
        overlap_meta = gene_meta.get(overlap_gene, {})

        out_rows.append(
            {
                "figure_context": "Figure 5A positive network",
                "snp_id": snp_id,
                "trial": row["trial"],
                "associated_biomarker_n": row["target_n"],
                "associated_biomarkers": row["targets"],
                "target_organisms": row["target_organisms"],
                "chromosome": row["chromosome"],
                "position_bp_v4": row["position_bp"],
                "alleles": row["alleles"],
                "overlap_gene_n": row["overlap_gene_n"],
                "overlapping_gene_id": overlap_gene,
                "overlapping_gene_confidence": overlap_meta.get("confidence", ""),
                "overlapping_gene_biotype": overlap_meta.get("biotype", ""),
                "overlapping_gene_description": overlap_meta.get("description", ""),
                "nearest_gene_id": row["nearest_gene_id"],
                "nearest_gene_distance_bp": row["nearest_gene_distance_bp"],
                "genes_within_50kb_n": row["window_50kb_gene_n"],
                "genes_within_100kb_n": row["window_100kb_gene_n"],
                "genes_within_250kb_n": row["window_250kb_gene_n"],
                "primary_candidate_gene_id": primary_gene,
                "primary_candidate_basis": primary_basis,
                "primary_candidate_confidence": primary_meta.get("confidence", ""),
                "primary_candidate_biotype": primary_meta.get("biotype", ""),
                "primary_candidate_description": primary_meta.get("description", ""),
                "primary_candidate_GO": primary_meta.get("ontology", ""),
                "primary_candidate_Pfam": primary_meta.get("pfam", ""),
                "annotation_note": (
                    "Direct overlap with annotated gene."
                    if int(row["overlap_gene_n"]) > 0
                    else "No overlapping gene; candidate assigned from nearest/flanking annotated gene."
                ),
            }
        )

    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
