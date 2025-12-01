# scripts/utils.py
import subprocess
from pathlib import Path
import yaml
from typing import Tuple, Dict, Any, List, Optional


CONFIG_DIR = Path("config")
PATHS_FILE = CONFIG_DIR / "paths.yaml"
PARAMS_FILE = CONFIG_DIR / "params.yaml"
LOG_DIR = Path("logs")
LOG_DIR.mkdir(parents=True, exist_ok=True)


def load_yaml(path: Path) -> Dict[str, Any]:
    with open(path) as f:
        return yaml.safe_load(f)


def load_paths_and_params() -> Tuple[Dict[str, Any], Dict[str, Any]]:
    paths = load_yaml(PATHS_FILE)
    params = load_yaml(PARAMS_FILE)
    return paths, params


def run_cmd(cmd: List[str],
            log_name: Optional[str] = None,
            check: bool = True) -> None:
    """
    공통적으로 shell 명령을 실행하는 함수.
    cmd: ["fastqc", "-t", "8", ...] 처럼 리스트로 넘김.
    log_name: logs 폴더에 저장할 로그 파일 이름 (예: '01_qc.log')
    """
    cmd_str = " ".join(cmd)
    print(f"[RUN] {cmd_str}")

    if log_name:
        log_path = LOG_DIR / log_name
        with open(log_path, "a") as log:
            log.write(f"\n\n===== CMD =====\n{cmd_str}\n")
            log.flush()
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True
            )
    else:
        result = subprocess.run(cmd)

    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed with code {result.returncode}: {cmd_str}")
