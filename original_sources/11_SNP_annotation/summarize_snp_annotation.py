import csv
import re
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent / "Data_for_Menghui" / "Menghui_paper_LD"
OUT_DIR = Path(__file__).resolve().parent


def parse_details(details: str) -> dict:
    fields = {}
    for item in details.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            fields[key] = value
    return fields


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    snp_rows = read_csv(BASE_DIR / "SNP4_annotations_with_Block_info.csv")

    no_gene_by_block = {}
    gene_by_block = {}

    for chr_dir in sorted(BASE_DIR.glob("Chr*")):
        chr_name = chr_dir.name

        gene_file = chr_dir / f"{chr_name}_Significant_SNPs_with_genes.csv"
        if gene_file.exists():
            for row in read_csv(gene_file):
                details = parse_details(row["Details"])
                gene_by_block.setdefault((chr_name, row["Block"]), []).append(
                    {
                        "block_start": row["Block_start"],
                        "block_end": row["Block_end"],
                        "marker_start": row["Marker_start"],
                        "marker_end": row["Marker_end"],
                        "gene_id": details.get("ID", ""),
                        "gene_name": details.get("Name", ""),
                        "confidence": details.get("confidence", ""),
                        "biotype": details.get("biotype", ""),
                    }
                )

        no_gene_file = chr_dir / f"{chr_name}_Significant_SNPs_with_no_genes.csv"
        if no_gene_file.exists():
            for row in read_csv(no_gene_file):
                no_gene_by_block[(chr_name, row["Block"])] = {
                    "block_start": row["Start"],
                    "block_end": row["End"],
                    "block_size": row["Block_Size"],
                    "snp_count_in_block": row["Count"],
                }

    candidate_rows = []
    no_gene_rows = []

    for snp in snp_rows:
        key = (snp["Chr"], snp["Block"])
        gene_hits = gene_by_block.get(key, [])

        if gene_hits:
            block_span = ""
            try:
                block_span = str(int(gene_hits[0]["block_end"]) - int(gene_hits[0]["block_start"]))
            except ValueError:
                pass

            for gene in gene_hits:
                candidate_rows.append(
                    {
                        "snp_id": snp["rs"],
                        "chromosome": snp["Chr"],
                        "block": snp["Block"],
                        "block_start": gene["block_start"],
                        "block_end": gene["block_end"],
                        "block_span_bp": block_span,
                        "gene_id": gene["gene_id"],
                        "gene_name": gene["gene_name"],
                        "gene_confidence": gene["confidence"],
                        "gene_biotype": gene["biotype"],
                        "gene_start": gene["marker_start"],
                        "gene_end": gene["marker_end"],
                    }
                )
        else:
            no_gene_info = no_gene_by_block.get(key, {})
            no_gene_rows.append(
                {
                    "snp_id": snp["rs"],
                    "chromosome": snp["Chr"],
                    "block": snp["Block"],
                    "block_start": no_gene_info.get("block_start", ""),
                    "block_end": no_gene_info.get("block_end", ""),
                    "block_span_bp": no_gene_info.get("block_size", ""),
                    "snp_count_in_block": no_gene_info.get("snp_count_in_block", ""),
                }
            )

    candidate_rows.sort(key=lambda row: (row["chromosome"], int(row["block"]), row["snp_id"], row["gene_id"]))
    no_gene_rows.sort(key=lambda row: (row["chromosome"], int(row["block"]), row["snp_id"]))

    candidate_path = OUT_DIR / "candidate_gene_table.csv"
    with candidate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(candidate_rows[0].keys()))
        writer.writeheader()
        writer.writerows(candidate_rows)

    no_gene_path = OUT_DIR / "snps_without_gene_overlap.csv"
    with no_gene_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(no_gene_rows[0].keys()))
        writer.writeheader()
        writer.writerows(no_gene_rows)

    print(f"Wrote {candidate_path}")
    print(f"Wrote {no_gene_path}")
    print(f"Candidate rows: {len(candidate_rows)}")
    print(f"SNPs without gene overlap: {len(no_gene_rows)}")


if __name__ == "__main__":
    main()
