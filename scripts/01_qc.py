# scripts/01_qc.py
from pathlib import Path
from utils import load_paths_and_params, run_cmd


def run_fastqc():
    paths, params = load_paths_and_params()

    fastq_dir = Path(paths["fastq_dir"])
    if not fastq_dir.exists():
        print(f"[QC] FASTQ directory does not exist: {fastq_dir}")
        return

    out_dir = Path(paths.get("output_dir", "results")) / "qc"
    out_dir.mkdir(parents=True, exist_ok=True)

    threads = params.get("qc", {}).get("threads", 4)

    fastq_files = sorted(
        list(fastq_dir.glob("*.fastq")) +
        list(fastq_dir.glob("*.fastq.gz")) +
        list(fastq_dir.glob("*.fq")) +
        list(fastq_dir.glob("*.fq.gz"))
    )

    if not fastq_files:
        print(f"[QC] No FASTQ files found in {fastq_dir}")
        return

    print(f"[QC] Found {len(fastq_files)} FASTQ files.")

    for fq in fastq_files:
        cmd = [
            "fastqc",
            "-t", str(threads),
            "-o", str(out_dir),
            str(fq)
        ]
        run_cmd(cmd, log_name=f"01_qc_fastqc_{fq.stem}.log")

    print(f"[QC] FastQC finished. Results in: {out_dir}")


if __name__ == "__main__":
    run_fastqc()
