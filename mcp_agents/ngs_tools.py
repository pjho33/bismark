from mcp.server.fastmcp import FastMCP
import subprocess
import pysam
import os

# 서버 이름 정의
mcp = FastMCP("NGS-Bio-Toolbox")

@mcp.tool()
def get_bam_flagstat(bam_path: str) -> str:
    """
    BAM 파일 경로를 받아 samtools flagstat을 실행하고 결과를 반환합니다.
    매핑 통계, QC 등의 정보를 확인할 때 사용합니다.
    """
    if not os.path.exists(bam_path):
        return f"Error: 파일을 찾을 수 없습니다: {bam_path}"

    try:
        # samtools가 시스템 PATH에 있어야 합니다.
        # 바이오마커 발굴 시 매핑 퀄리티 확인용
        result = subprocess.run(
            ["samtools", "flagstat", bam_path],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"Error running samtools: {e.stderr}"
    except FileNotFoundError:
        return "Error: 'samtools' 명령어를 찾을 수 없습니다. 설치되어 있는지 확인해주세요."

@mcp.tool()
def query_vcf_region(vcf_path: str, chrom: str, start: int, end: int) -> str:
    """
    VCF 파일에서 특정 영역(Locus)의 변이 정보를 조회합니다.
    pysam을 사용하여 빠르게 해당 영역만 fetch합니다.
    결과는 사람이 읽기 쉬운 요약 텍스트로 반환됩니다.
    
    Args:
        vcf_path: VCF 파일 경로 (반드시 인덱싱(.tbi) 되어 있어야 함)
        chrom: 염색체 번호 (예: 'chr1', '1')
        start: 시작 위치 (1-based coordinate 권장)
        end: 끝 위치
    """
    if not os.path.exists(vcf_path):
        return f"Error: 파일을 찾을 수 없습니다: {vcf_path}"

    output_lines = []
    
    try:
        # pysam을 이용한 VCF 읽기 (인덱싱 필수)
        vcf = pysam.VariantFile(vcf_path)
        
        # 해당 영역 fetch
        # 주의: pysam은 0-based, half-open interval을 사용하지만, 
        # 사용 편의를 위해 입력받은 좌표를 그대로 넘기거나 필요시 보정합니다.
        # 여기서는 안전하게 pysam의 fetch 기능을 그대로 씁니다.
        
        count = 0
        MAX_RECORDS = 50 # Context Window 보호를 위해 최대 조회 개수 제한
        
        output_lines.append(f"Querying {vcf_path} at {chrom}:{start}-{end}\n")
        
        for record in vcf.fetch(chrom, start, end):
            # 필요한 정보만 포맷팅 (Chrom, Pos, Ref, Alt, Filter, Quality)
            alts = ",".join(record.alts) if record.alts else "."
            filt = ",".join(record.filter.keys()) if record.filter.keys() else "PASS"
            
            # Info 필드에서 중요한 것만 뽑을 수도 있습니다 (예: AF, DP)
            dp = record.info.get('DP', 'NA')
            af = record.info.get('AF', 'NA')
            
            line = (f"[Pos: {record.pos}] {record.ref} -> {alts} | "
                    f"Qual: {record.qual} | Filter: {filt} | DP: {dp} | AF: {af}")
            output_lines.append(line)
            
            count += 1
            if count >= MAX_RECORDS:
                output_lines.append(f"\n... (Too many variants. Stopping at {MAX_RECORDS} records)")
                break
        
        vcf.close()
        
        if count == 0:
            return "No variants found in this region."
            
        return "\n".join(output_lines)

    except ValueError as e:
        return f"Error parsing VCF (인덱스 파일 .tbi가 있는지 확인하세요): {str(e)}"
    except Exception as e:
        return f"Unexpected error: {str(e)}"

if __name__ == "__main__":
    # Windsurf가 stdio를 통해 이 스크립트와 통신합니다.
    mcp.run(transport='stdio')