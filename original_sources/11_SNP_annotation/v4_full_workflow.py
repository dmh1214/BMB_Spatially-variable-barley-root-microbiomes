import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional


WORKDIR = Path(__file__).resolve().parent
ROOT = WORKDIR.parent

V4_HAPMAP = Path(
    "/Users/gaoqi/Library/CloudStorage/OneDrive-Aarhusuniversitet/BarleyMicroBreed/data/phenotype/BMB_Geno_Hapmap v4.txt"
)
GFF_PATH = WORKDIR / "Data_for_Menghui" / "LD_analysis" / "Morex_gff_file_all_chromosomes.txt"
EDGE_PATH = ROOT / "Network" / "mgwas_sum_edge.csv"

PANEL_HAPMAPS = {
    "SpringPan": ROOT / "Bac_Aut" / "mGWAS" / "bac_spring_snp_hapmap.txt",
    "WinterPan_Bac": ROOT / "Bac_Aut" / "mGWAS" / "bac_winter_snp_hapmap.txt",
    "WinterPan_Fun": ROOT / "Fun_Aut" / "mGWAS" / "fun_winter_snp_hapmap.txt",
}

CHROM_MAP = {
    "1": "chr1H",
    "2": "chr2H",
    "3": "chr3H",
    "4": "chr4H",
    "5": "chr5H",
    "6": "chr6H",
    "7": "chr7H",
}

FOCUS_SNPS = {
    "JHI-Hv50k-2016-201999",
    "JHI-Hv50k-2016-27968",
    "JHI-Hv50k-2016-47293",
    "SCRI_RS_149566",
}

WINDOWS = [50_000, 100_000, 250_000]
LD_RADIUS = 1_000_000
LD_R2_THRESHOLD = 0.2


def normalize_sample(name: str) -> str:
    return name.strip().lower()


def parse_gene_attributes(attr: str) -> dict:
    values = {}
    for field in attr.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
    return values


def read_panel_samples(path: Path) -> list[str]:
    with path.open() as handle:
        header = next(handle).rstrip("\n").split("\t")
    return [normalize_sample(x) for x in header[11:]]


def determine_target_snps() -> list[dict]:
    selected = {}
    with EDGE_PATH.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    source_col = "\ufeffSource" if "\ufeffSource" in rows[0] else "Source"
    for row in rows:
        snp = row[source_col]
        target = row["Target"]
        trial = row["Trial"]
        positive = row["shap_sign"] == "Positive"
        organism = "Bacteria" if target.startswith("ASV") else "Fungi"

        if positive:
            record = selected.setdefault(
                snp,
                {
                    "snp_id": snp,
                    "selection_reason": "Figure5A_positive_network",
                    "trial": trial,
                    "targets": [],
                    "target_organisms": set(),
                },
            )
            record["targets"].append(target)
            record["target_organisms"].add(organism)

        if snp in FOCUS_SNPS and snp not in selected:
            record = selected.setdefault(
                snp,
                {
                    "snp_id": snp,
                    "selection_reason": "manuscript_focus",
                    "trial": trial,
                    "targets": [],
                    "target_organisms": set(),
                },
            )
            record["targets"].append(target)
            record["target_organisms"].add(organism)

    out = []
    for record in selected.values():
        organisms = sorted(record["target_organisms"])
        if record["trial"] == "SpringPan":
            panel_key = "SpringPan"
        elif organisms == ["Bacteria"]:
            panel_key = "WinterPan_Bac"
        else:
            panel_key = "WinterPan_Fun"
        out.append(
            {
                "snp_id": record["snp_id"],
                "selection_reason": record["selection_reason"],
                "trial": record["trial"],
                "panel_key": panel_key,
                "target_n": len(record["targets"]),
                "targets": ";".join(sorted(record["targets"])),
                "target_organisms": ";".join(organisms),
            }
        )
    return sorted(out, key=lambda row: (row["trial"], row["snp_id"]))


def read_gff_genes() -> dict[str, list[dict]]:
    genes_by_chr = defaultdict(list)
    with GFF_PATH.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "gene":
                continue
            attrs = parse_gene_attributes(parts[8])
            genes_by_chr[parts[0]].append(
                {
                    "chrom_gff": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "gene_id": attrs.get("ID", ""),
                    "gene_name": attrs.get("Name", ""),
                    "confidence": attrs.get("confidence", ""),
                    "biotype": attrs.get("biotype", ""),
                }
            )
    for genes in genes_by_chr.values():
        genes.sort(key=lambda row: row["start"])
    return genes_by_chr


def load_v4_rows(target_ids: set[str], sample_lookup: dict[str, list[int]]) -> tuple[dict, list[dict]]:
    target_rows = {}
    nearby_rows = []
    with V4_HAPMAP.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_names = header[11:]
        index_by_sample = {normalize_sample(name): idx for idx, name in enumerate(sample_names)}
        for sample in sample_lookup:
            if sample not in index_by_sample:
                raise KeyError(f"Sample {sample} not present in v4 hapmap")
        for row in reader:
            snp_id = row[0]
            chrom = row[2]
            pos = int(row[3])
            alleles = row[1]
            full_row = {
                "snp_id": snp_id,
                "chrom": chrom,
                "pos": pos,
                "alleles": alleles,
                "genotypes": row[11:],
            }
            if snp_id in target_ids:
                target_rows[snp_id] = full_row
            if chrom in {"1", "3"}:
                nearby_rows.append(full_row)
    return target_rows, nearby_rows


def genotype_to_dosage(gt: str, alt: str) -> Optional[int]:
    gt = gt.upper()
    if gt in {"N", "NN", "NA", ""}:
        return None
    if len(gt) == 1:
        return 2 if gt == alt else 0
    if len(gt) == 2 and gt[0] == gt[1]:
        return 2 if gt[0] == alt else 0
    return None


def dosage_vector(row: dict, sample_positions: list[int]) -> list[Optional[int]]:
    parts = row["alleles"].split("/")
    if len(parts) != 2:
        return [None] * len(sample_positions)
    alt = parts[1].upper()
    return [genotype_to_dosage(row["genotypes"][idx], alt) for idx in sample_positions]


def compute_r2(x: list[Optional[int]], y: list[Optional[int]]) -> Optional[float]:
    paired = [(a, b) for a, b in zip(x, y) if a is not None and b is not None]
    if len(paired) < 10:
        return None
    xs = [a for a, _ in paired]
    ys = [b for _, b in paired]
    if len(set(xs)) < 2 or len(set(ys)) < 2:
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    cov = sum((a - mean_x) * (b - mean_y) for a, b in paired)
    var_x = sum((a - mean_x) ** 2 for a in xs)
    var_y = sum((b - mean_y) ** 2 for b in ys)
    if var_x == 0 or var_y == 0:
        return None
    r = cov / math.sqrt(var_x * var_y)
    return r * r


def find_flanking_genes(genes: list[dict], pos: int) -> tuple[list[dict], Optional[dict], Optional[dict]]:
    overlaps = [gene for gene in genes if gene["start"] <= pos <= gene["end"]]
    left = None
    right = None
    for gene in genes:
        if gene["end"] < pos:
            left = gene
        elif gene["start"] > pos:
            right = gene
            break
    return overlaps, left, right


def genes_within_window(genes: list[dict], pos: int, window: int) -> list[dict]:
    lower = pos - window
    upper = pos + window
    return [gene for gene in genes if gene["end"] >= lower and gene["start"] <= upper]


def best_gene(genes: list[dict]) -> Optional[dict]:
    if not genes:
        return None
    ranked = sorted(
        genes,
        key=lambda gene: (
            gene["confidence"] != "high",
            gene["biotype"] != "protein_coding",
            gene["start"],
        ),
    )
    return ranked[0]


def main() -> None:
    outdir = WORKDIR / "v4_full_workflow"
    outdir.mkdir(exist_ok=True)

    target_snp_rows = determine_target_snps()
    panel_samples = {key: read_panel_samples(path) for key, path in PANEL_HAPMAPS.items()}

    sample_positions_lookup = {}
    # Build from the master v4 header once.
    with V4_HAPMAP.open() as handle:
        header = next(handle).rstrip("\n").split("\t")
    v4_samples = header[11:]
    v4_index = {normalize_sample(name): idx for idx, name in enumerate(v4_samples)}
    for panel_key, samples in panel_samples.items():
        sample_positions_lookup[panel_key] = [v4_index[sample] for sample in samples if sample in v4_index]

    genes_by_chr = read_gff_genes()
    target_ids = {row["snp_id"] for row in target_snp_rows}
    target_rows, nearby_rows = load_v4_rows(target_ids, {})

    metadata_rows = []
    candidate_rows = []
    ld_rows = []

    nearby_by_chr = defaultdict(list)
    for row in nearby_rows:
        nearby_by_chr[row["chrom"]].append(row)
    for rows in nearby_by_chr.values():
        rows.sort(key=lambda row: row["pos"])

    for snp_info in target_snp_rows:
        snp_id = snp_info["snp_id"]
        row = target_rows[snp_id]
        chrom_gff = CHROM_MAP[row["chrom"]]
        genes = genes_by_chr[chrom_gff]

        overlaps, left_gene, right_gene = find_flanking_genes(genes, row["pos"])
        nearest_gene = None
        nearest_distance = ""
        if overlaps:
            nearest_gene = overlaps[0]
            nearest_distance = 0
        else:
            candidates = []
            if left_gene:
                candidates.append((row["pos"] - left_gene["end"], left_gene))
            if right_gene:
                candidates.append((right_gene["start"] - row["pos"], right_gene))
            if candidates:
                nearest_distance, nearest_gene = min(candidates, key=lambda item: item[0])

        window_hits = {}
        for window in WINDOWS:
            genes_in_window = genes_within_window(genes, row["pos"], window)
            window_hits[window] = genes_in_window
            for gene in genes_in_window:
                candidate_rows.append(
                    {
                        "snp_id": snp_id,
                        "selection_reason": snp_info["selection_reason"],
                        "trial": snp_info["trial"],
                        "window_bp": window,
                        "chromosome": row["chrom"],
                        "snp_pos": row["pos"],
                        "gene_id": gene["gene_id"],
                        "gene_name": gene["gene_name"],
                        "gene_start": gene["start"],
                        "gene_end": gene["end"],
                        "distance_to_snp_bp": 0
                        if gene["start"] <= row["pos"] <= gene["end"]
                        else min(abs(row["pos"] - gene["start"]), abs(row["pos"] - gene["end"])),
                        "confidence": gene["confidence"],
                        "biotype": gene["biotype"],
                        "strand": gene["strand"],
                    }
                )

        panel_key = snp_info["panel_key"]
        sample_positions = sample_positions_lookup[panel_key]
        target_dosage = dosage_vector(row, sample_positions)

        local_r2 = []
        for other in nearby_by_chr[row["chrom"]]:
            if abs(other["pos"] - row["pos"]) > LD_RADIUS:
                continue
            other_dosage = dosage_vector(other, sample_positions)
            r2 = compute_r2(target_dosage, other_dosage)
            if r2 is None:
                continue
            local_r2.append((other["snp_id"], other["pos"], r2))

        local_r2.sort(key=lambda item: item[1])
        support = [item for item in local_r2 if item[2] >= LD_R2_THRESHOLD]
        support_start = min((item[1] for item in support), default="")
        support_end = max((item[1] for item in support), default="")
        peak_partner = max(local_r2, key=lambda item: item[2]) if local_r2 else None

        for partner_id, partner_pos, r2 in sorted(local_r2, key=lambda item: item[2], reverse=True)[:20]:
            ld_rows.append(
                {
                    "snp_id": snp_id,
                    "trial": snp_info["trial"],
                    "panel_key": panel_key,
                    "partner_snp": partner_id,
                    "partner_pos": partner_pos,
                    "distance_bp": partner_pos - row["pos"],
                    "r2": round(r2, 6),
                }
            )

        best_50 = best_gene(window_hits[50_000])
        best_100 = best_gene(window_hits[100_000])
        best_250 = best_gene(window_hits[250_000])

        metadata_rows.append(
            {
                "snp_id": snp_id,
                "selection_reason": snp_info["selection_reason"],
                "trial": snp_info["trial"],
                "panel_key": panel_key,
                "target_n": snp_info["target_n"],
                "targets": snp_info["targets"],
                "target_organisms": snp_info["target_organisms"],
                "chromosome": row["chrom"],
                "position_bp": row["pos"],
                "alleles": row["alleles"],
                "overlap_gene_n": len(overlaps),
                "left_gene_id": left_gene["gene_id"] if left_gene else "",
                "left_gene_distance_bp": row["pos"] - left_gene["end"] if left_gene else "",
                "right_gene_id": right_gene["gene_id"] if right_gene else "",
                "right_gene_distance_bp": right_gene["start"] - row["pos"] if right_gene else "",
                "nearest_gene_id": nearest_gene["gene_id"] if nearest_gene else "",
                "nearest_gene_distance_bp": nearest_distance,
                "window_50kb_gene_n": len(window_hits[50_000]),
                "window_100kb_gene_n": len(window_hits[100_000]),
                "window_250kb_gene_n": len(window_hits[250_000]),
                "best_50kb_gene": best_50["gene_id"] if best_50 else "",
                "best_100kb_gene": best_100["gene_id"] if best_100 else "",
                "best_250kb_gene": best_250["gene_id"] if best_250 else "",
                "ld_support_start_bp_r2_ge_0_2": support_start,
                "ld_support_end_bp_r2_ge_0_2": support_end,
                "ld_support_span_bp_r2_ge_0_2": support_end - support_start if support_start != "" else "",
                "ld_partner_count_within_1mb": len(local_r2),
                "ld_partner_count_r2_ge_0_2": len(support),
                "peak_ld_partner": peak_partner[0] if peak_partner else "",
                "peak_ld_partner_pos": peak_partner[1] if peak_partner else "",
                "peak_ld_r2": round(peak_partner[2], 6) if peak_partner else "",
            }
        )

    def write_csv(path: Path, rows: list[dict]) -> None:
        if not rows:
            return
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    write_csv(outdir / "figure5a_v4_snp_summary.csv", metadata_rows)
    write_csv(outdir / "figure5a_v4_candidate_genes_by_window.csv", candidate_rows)
    write_csv(outdir / "figure5a_v4_top_ld_partners.csv", ld_rows)

    print(f"Wrote {outdir / 'figure5a_v4_snp_summary.csv'}")
    print(f"Wrote {outdir / 'figure5a_v4_candidate_genes_by_window.csv'}")
    print(f"Wrote {outdir / 'figure5a_v4_top_ld_partners.csv'}")


if __name__ == "__main__":
    main()
