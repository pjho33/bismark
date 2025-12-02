# scripts/04_dedup.py
from pathlib import Path
from utils import load_paths_and_params, run_cmd


def run_deduplicate():
    paths, _ = load_paths_and_params()

    output_root = Path(paths.get("output_dir", "results"))
    bam_dir = output_root / "bam"
    if not bam_dir.exists():
        print(f"[DEDUP] BAM directory does not exist: {bam_dir}")
        return

    bam_files = sorted(list(bam_dir.glob("*.bam")))

    if not bam_files:
        print(f"[DEDUP] No BAM files found in {bam_dir}")
        return

    print(f"[DEDUP] Found {len(bam_files)} BAMs")

    for bam in bam_files:
        cmd = [
            "deduplicate_bismark",
            "--bam",
            "--output_dir", str(bam_dir),
            str(bam)
        ]
        run_cmd(cmd, log_name=f"04_dedup_{bam.stem}.log")

    print("[DEDUP] Deduplication finished. Dedup BAMs (*.deduplicated.bam) created.")


if __name__ == "__main__":
    run_deduplicate()
