# scripts/03_align.py
from pathlib import Path
from utils import load_paths_and_params, run_cmd


def pair_trimmed(trim_dir: Path):
    """
    Trim Galore 결과의 val_1 / val_2를 매칭.
    예: sample_R1_val_1.fq.gz, sample_R2_val_2.fq.gz
    -> base: sample
    """
    r1_files = sorted(list(trim_dir.glob("*val_1.fq")) +
                      list(trim_dir.glob("*val_1.fq.gz")))
    pairs = []

    for r1 in r1_files:
        name = r1.name
        base = name.replace("_R1_val_1.fq.gz", "").replace("_R1_val_1.fq", "")
        # R2 후보
        candidates_r2 = [
            trim_dir / f"{base}_R2_val_2.fq",
            trim_dir / f"{base}_R2_val_2.fq.gz",
        ]
        r2 = None
        for c in candidates_r2:
            if c.exists():
                r2 = c
                break
        if r2 is None:
            print(f"[ALIGN] WARNING: R2 not found for {r1.name}, skipping")
            continue
        pairs.append((base, r1, r2))

    return pairs


def run_bismark_align():
    paths, params = load_paths_and_params()

    output_root = Path(paths.get("output_dir", "results"))
    trim_dir = output_root / "trimmed"
    bam_dir = output_root / "bam"
    bam_dir.mkdir(parents=True, exist_ok=True)

    threads = params.get("align", {}).get("threads", 8)

    ref_fa = Path(paths["reference_genome"])
    genome_folder = ref_fa.parent  # 예: data/reference

    pairs = pair_trimmed(trim_dir)
    if not pairs:
        print(f"[ALIGN] No trimmed paired files found in {trim_dir}")
        return

    print(f"[ALIGN] Genome folder: {genome_folder}")
    print(f"[ALIGN] Found {len(pairs)} samples")

    for sample, r1, r2 in pairs:
        cmd = [
            "bismark",
            str(genome_folder),
            "-1", str(r1),
            "-2", str(r2),
            "-p", str(threads),
            "--output_dir", str(bam_dir),
        ]
        # 필요시 옵션 추가 (e.g. --gzip, --non_directional 등)
        run_cmd(cmd, log_name="03_bismark_align.log")

    print(f"[ALIGN] Bismark alignment finished. BAMs in: {bam_dir}")


if __name__ == "__main__":
    run_bismark_align()
