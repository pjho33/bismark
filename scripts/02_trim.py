# scripts/02_trim.py
from pathlib import Path
from utils import load_paths_and_params, run_cmd


def pair_fastq(fastq_dir: Path):
    """
    R1/R2 짝을 맞춰서 (sample, R1, R2) 리스트를 반환.
    규칙: *_R1.fastq(.gz), *_R2.fastq(.gz)
    """
    r1_files = sorted(list(fastq_dir.glob("*_R1.fastq")) +
                      list(fastq_dir.glob("*_R1.fastq.gz")))
    pairs = []

    for r1 in r1_files:
        base = r1.name.replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        candidates_r2 = [
            fastq_dir / f"{base}_R2.fastq",
            fastq_dir / f"{base}_R2.fastq.gz",
        ]
        r2 = None
        for c in candidates_r2:
            if c.exists():
                r2 = c
                break
        if r2 is None:
            print(f"[TRIM] WARNING: R2 not found for {r1.name}, skipping")
            continue
        pairs.append((base, r1, r2))

    return pairs


def run_trim_galore():
    paths, params = load_paths_and_params()

    fastq_dir = Path(paths["fastq_dir"])
    out_dir = Path(paths.get("output_dir", "results")) / "trimmed"
    out_dir.mkdir(parents=True, exist_ok=True)

    threads = params.get("trim", {}).get("threads", 8)
    adapter = params.get("trim", {}).get("adapter", None)

    pairs = pair_fastq(fastq_dir)
    if not pairs:
        print(f"[TRIM] No paired FASTQ found in {fastq_dir}")
        return

    print(f"[TRIM] Found {len(pairs)} paired samples")

    for sample, r1, r2 in pairs:
        cmd = [
            "trim_galore",
            "--paired",
            f"--cores={threads}",
            "-o", str(out_dir),
        ]
        if adapter:
            cmd += ["--adapter", adapter]

        cmd += [str(r1), str(r2)]

        run_cmd(cmd, log_name="02_trim_galore.log")

    print(f"[TRIM] Trim Galore finished. Results in: {out_dir}")


if __name__ == "__main__":
    run_trim_galore()
