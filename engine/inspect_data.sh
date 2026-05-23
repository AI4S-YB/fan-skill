#!/bin/bash
# inspect_data.sh — 探查输入数据特征，输出 data_profile.json
# 用法: bash inspect_data.sh <input_path> [format_hint]
# 输出: JSON 到 stdout

INPUT="$1"
FORMAT_HINT="${2:-auto}"

# 检测格式
detect_format() {
    local path="$1"
    if [ -f "${path}.bed" ] || [ -f "${path}.bim" ] || [ -f "${path}.fam" ]; then
        echo "bed"
    elif [ -f "$path" ] && [[ "$path" == *.vcf.gz || "$path" == *.vcf ]]; then
        echo "vcf"
    elif [ -f "$path" ] && [[ "$path" == *.hmp.txt || "$path" == *.hapmap ]]; then
        echo "hapmap"
    else
        echo "unknown"
    fi
}

# 从染色体名推断物种
infer_species() {
    local bim_file="$1"
    if [ ! -f "$bim_file" ]; then
        echo "unknown"
        return
    fi
    local chr_names
    chr_names=$(cut -f1 "$bim_file" 2>/dev/null | sort -u | head -20)

    if echo "$chr_names" | grep -qE '^chr[0-9]+$'; then
        local max_chr=$(echo "$chr_names" | sed 's/chr//' | sort -n | tail -1)
        case "$max_chr" in
            10) echo "zea_mays" ;;
            12) echo "oryza_sativa" ;;
            5)  echo "arabidopsis_thaliana" ;;
            *)  echo "unknown" ;;
        esac
    elif echo "$chr_names" | grep -qE '^Chr0[0-9]|^Chr1[0-2]$'; then
        echo "oryza_sativa"
    elif echo "$chr_names" | grep -qE '^[0-9]+[ABD]$'; then
        echo "triticum_aestivum"
    elif echo "$chr_names" | grep -qE '^[0-9]+$'; then
        local max_chr=$(echo "$chr_names" | sort -n | tail -1)
        case "$max_chr" in
            10) echo "zea_mays" ;;
            12) echo "oryza_sativa" ;;
            20) echo "glycine_max" ;;
            5)  echo "arabidopsis_thaliana" ;;
            *)  echo "unknown" ;;
        esac
    else
        echo "unknown"
    fi
}

# 主逻辑
FORMAT="$FORMAT_HINT"
if [ "$FORMAT" = "auto" ]; then
    FORMAT=$(detect_format "$INPUT")
fi

# 获取统计量
SNP_COUNT=0
SAMPLE_COUNT=0
SPECIES="unknown"

if [ "$FORMAT" = "bed" ]; then
    [ -f "${INPUT}.bim" ] && SNP_COUNT=$(wc -l < "${INPUT}.bim" | tr -d ' ')
    [ -f "${INPUT}.fam" ] && SAMPLE_COUNT=$(wc -l < "${INPUT}.fam" | tr -d ' ')
    SPECIES=$(infer_species "${INPUT}.bim")
elif [ "$FORMAT" = "vcf" ]; then
    if [[ "$INPUT" == *.vcf.gz ]]; then
        SNP_COUNT=$(zgrep -v '^#' "$INPUT" 2>/dev/null | wc -l | tr -d ' ')
    else
        SNP_COUNT=$(grep -v '^#' "$INPUT" 2>/dev/null | wc -l | tr -d ' ')
    fi
    # 从 VCF 头推断物种（简化）
    SPECIES="unknown"
fi

# 评估标记密度
DENSITY="unknown"
if [ "$SNP_COUNT" -gt 0 ] 2>/dev/null; then
    if [ "$SNP_COUNT" -lt 1000 ]; then
        DENSITY="very_low"
    elif [ "$SNP_COUNT" -lt 10000 ]; then
        DENSITY="low"
    elif [ "$SNP_COUNT" -lt 100000 ]; then
        DENSITY="medium"
    else
        DENSITY="high"
    fi
fi

# 输出 JSON
python3 -c "
import json
data = {
    'input_path': '${INPUT}',
    'format': '${FORMAT}',
    'snp_count': ${SNP_COUNT:-0},
    'sample_count': ${SAMPLE_COUNT:-0},
    'species_detected': '${SPECIES}',
    'species_source': 'inferred_from_chromosome_names' if '${SPECIES}' != 'unknown' else 'none',
    'ploidy': 'unknown',
    'breeding_system': 'unknown',
    'chromosome_count': 0,
    'subgenome_info': None,
    'snp_density_quality': '${DENSITY}'
}
print(json.dumps(data, indent=2))
"
