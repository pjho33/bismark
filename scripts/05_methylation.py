# scripts/05_methylation.py
from pathlib import Path
from utils import load_paths_and_params, run_cmd


def run_methylation_extractor():
    paths, params = load_paths_and_params()

    output_root = Path(paths.get("output_dir", "results"))
    bam_dir = output_root / "bam"
    if not bam_dir.exists():
        print(f"[METH] BAM directory does not exist: {bam_dir}")
        return

    meth_dir = output_root / "methylation"
    meth_dir.mkdir(parents=True, exist_ok=True)

    threads = params.get("methylation", {}).get("threads", 8)

    bam_files = sorted(bam_dir.glob("*.deduplicated.bam"))
    if not bam_files:
        print(f"[METH] No deduplicated BAMs found in {bam_dir}")
        return

    print(f"[METH] Found {len(bam_files)} deduplicated BAMs")

    for bam in bam_files:
        cmd = [
            "bismark_methylation_extractor",
            "--bedGraph",
            "--gzip",
            "--multicore", str(threads),
            "--comprehensive",
            "--cytosine_report",
            "--genome_folder", str(Path(paths["reference_genome"]).parent),
            "-o", str(meth_dir),
            str(bam),
        ]
        run_cmd(cmd, log_name=f"05_methylation_extractor_{bam.stem}.log")

    print(f"[METH] Methylation extraction finished. Results in: {meth_dir}")


if __name__ == "__main__":
    run_methylation_extractor()
