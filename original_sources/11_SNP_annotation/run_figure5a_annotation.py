import subprocess
import sys
from pathlib import Path


WORKDIR = Path(__file__).resolve().parent

SCRIPTS = [
    WORKDIR / "v4_full_workflow.py",
    WORKDIR / "build_figure5_supplementary_table.py",
    WORKDIR / "reformat_figure5_supplementary_table.py",
]


def run_script(script: Path) -> None:
    print(f"[run] {script.name}")
    subprocess.run([sys.executable, str(script)], check=True)


def main() -> None:
    for script in SCRIPTS:
        run_script(script)

    print("\n[done] Core Figure 5A SNP annotation workflow completed.")
    print("[output] 11_SNP_annotation/v4_full_workflow/figure5a_v4_snp_summary.csv")
    print("[output] 11_SNP_annotation/v4_full_workflow/figure5a_v4_candidate_genes_by_window.csv")
    print("[output] 11_SNP_annotation/v4_full_workflow/figure5a_v4_top_ld_partners.csv")
    print("[output] 11_SNP_annotation/Supplementary_Table_Figure5A_10_SNP_annotations.csv")
    print("[output] 11_SNP_annotation/Supplementary_Table_Figure5A_10_SNP_annotations_clean.csv")


if __name__ == "__main__":
    main()
