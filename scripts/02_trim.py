# scripts/02_trim.py
from pathlib import Path
from utils import load_paths_and_params, run_cmd


def pair_fastq(fastq_dir: Path):
    """
    R1/R2 짝을 맞춰서 (sample, R1, R2) 리스트를 반환.
    지원 패턴:
      - *_R1.fastq(.gz), *_R2.fastq(.gz)
      - *_R1_*.fastq(.gz), *_R2_*.fastq(.gz)  (예: sample_R1_001.fastq.gz)
    """
    import re
    
    # R1 파일 찾기 (R1 뒤에 추가 문자가 있을 수 있음)
    all_files = list(fastq_dir.iterdir())
    r1_pattern = re.compile(r'(.+)_R1([_.].*)?\.f(ast)?q(\.gz)?$')
    
    r1_files = sorted([f for f in all_files if r1_pattern.match(f.name)])
    pairs = []

    for r1 in r1_files:
        # R1 -> R2 변환
        r2_name = re.sub(r'_R1([_.])', r'_R2\1', r1.name)
        r2 = fastq_dir / r2_name
        
        if not r2.exists():
            print(f"[TRIM] WARNING: R2 not found for {r1.name}, skipping")
            continue
        
        # 샘플명 추출 (R1 앞부분)
        match = r1_pattern.match(r1.name)
        base = match.group(1) if match else r1.stem
        pairs.append((base, r1, r2))

    return pairs


def run_trim_galore():
    paths, params = load_paths_and_params()

    fastq_dir = Path(paths["fastq_dir"])
    if not fastq_dir.exists():
        print(f"[TRIM] FASTQ directory does not exist: {fastq_dir}")
        return

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

        run_cmd(cmd, log_name=f"02_trim_galore_{sample}.log")

    print(f"[TRIM] Trim Galore finished. Results in: {out_dir}")


if __name__ == "__main__":
    run_trim_galore()
