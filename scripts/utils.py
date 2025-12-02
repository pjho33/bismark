# scripts/utils.py
import subprocess
import sys
from pathlib import Path

import yaml


def get_project_root() -> Path:
    """프로젝트 루트 디렉토리 반환 (scripts의 상위 폴더)"""
    return Path(__file__).resolve().parent.parent


def load_paths_and_params() -> tuple[dict, dict]:
    """config/paths.yaml과 config/params.yaml을 로드하여 반환"""
    root = get_project_root()
    
    paths_file = root / "config" / "paths.yaml"
    params_file = root / "config" / "params.yaml"
    
    if not paths_file.exists():
        raise FileNotFoundError(f"paths.yaml not found: {paths_file}")
    if not params_file.exists():
        raise FileNotFoundError(f"params.yaml not found: {params_file}")
    
    with open(paths_file, "r") as f:
        paths = yaml.safe_load(f) or {}
    
    with open(params_file, "r") as f:
        params = yaml.safe_load(f) or {}
    
    # 상대 경로를 프로젝트 루트 기준 절대 경로로 변환
    for key in ["fastq_dir", "reference_genome", "output_dir"]:
        if key in paths and not Path(paths[key]).is_absolute():
            paths[key] = str(root / paths[key])
    
    return paths, params


def run_cmd(cmd: list, log_name: str = None, cwd: str = None) -> int:
    """
    명령어 실행 및 로그 기록
    
    Args:
        cmd: 실행할 명령어 리스트
        log_name: 로그 파일 이름 (logs/ 디렉토리에 저장)
        cwd: 작업 디렉토리
    
    Returns:
        프로세스 반환 코드
    """
    root = get_project_root()
    log_dir = root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    cmd_str = " ".join(cmd)
    print(f"[CMD] {cmd_str}")
    
    log_path = log_dir / log_name if log_name else None
    
    try:
        if log_path:
            with open(log_path, "a") as log_file:
                result = subprocess.run(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    cwd=cwd,
                    text=True
                )
        else:
            result = subprocess.run(
                cmd,
                cwd=cwd,
                text=True
            )
        return result.returncode
    except FileNotFoundError:
        print(f"[ERROR] Command not found: {cmd[0]}")
        return 1
    except Exception as e:
        print(f"[ERROR] Failed to run command: {e}")
        return 1
