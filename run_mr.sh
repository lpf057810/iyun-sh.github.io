#!/bin/bash

# =============================================================================
# MR Analysis Unified Entry Script
# =============================================================================
# Supports: eQTL, pQTL, Traditional MR (selfharm + violence)
# Optional: Colocalization analysis
# Standardized: logs/, result/, report/ directories
# =============================================================================

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"
R_SCRIPT_DIR="$SCRIPT_DIR/scripts"
LOG_DIR="$SCRIPT_DIR/logs"
RESULT_DIR="$SCRIPT_DIR/results"
REPORT_DIR="$SCRIPT_DIR/report"

# Generate timestamp (YYYYMMDD_HHMMSS)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Default values
MR_TYPE=""
INPUT_FILE=""
OUTPUT_DIR=""
COLOC_ENABLED=0
PVAL_THRESHOLD="5e-8"
KB_DISTANCE="10000"
R2_THRESHOLD="0.001"
FSTAT_THRESHOLD="10"
EAF_THRESHOLD="0.01"
MIN_SNPS="3"
OUTCOME_ID="ebi-a-GCST90091033"
PQTLDIR="/media/desk16/iyunlpf/NAF/pqtl/data/cis_clump"
HELP=0

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print colored message
print_msg() {
    local color=$1
    local msg=$2
    echo -e "${color}${msg}${NC}"
}

# Log message to both console and log file
log_message() {
    local msg="$1"
    echo "[$(date '+%H:%M:%S')] $msg" | tee -a "$LOG_FILE"
}

# Show help
show_help() {
    cat << EOF
MR Analysis Unified Entry Script

USAGE: $(basename "$0") [OPTIONS]

OPTIONS:
    -t, --type TYPE       MR analysis type: eqtl, pqtl, traditional
    -i, --input FILE      Input file (gene list CSV)
    -o, --output DIR      Output directory (default: result/MR_TYPE)
    -c, --coloc NUM       Run colocalization analysis (0 or 1)
    -p, --pval VALUE      P-value threshold (default: 5e-8)
    -k, --kb VALUE        Distance in kb for clumping (default: 10000)
    -r, --r2 VALUE        LD r2 threshold (default: 0.001)
    -f, --fstat VALUE     F-statistic threshold (default: 10)
    -e, --eaf VALUE       EAF threshold (default: 0.01)
    -m, --minsnps VALUE   Minimum SNPs (default: 3)
    -d, --outcome VALUE   Outcome GWAS ID (default: ebi-a-GCST90091033)
    -q, --pqtl-dir VALUE pQTL data directory (default: /media/desk16/iyunlpf/NAF/pqtl/data/cis_clump)
    -h, --help            Show this help message

EXAMPLES:
    # Interactive mode - select MR type from menu
    $(basename "$0")

    # Run eQTL analysis
    $(basename "$0") -t eqtl -i /path/to/genes.csv -o result/MR_eqtl

    # Run pQTL analysis with colocalization
    $(basename "$0") -t pqtl -i /path/to/genes.csv -o result/MR_pqtl -c 1

    # Run traditional MR
    $(basename "$0") -t traditional -o result/MR_traditional

EOF
}

# Interactive menu
show_menu() {
    clear
    print_msg "$BLUE" "=============================================="
    print_msg "$BLUE" "       MR Analysis - Type Selection"
    print_msg "$BLUE" "=============================================="
    echo ""
    echo "Select MR Analysis Type:"
    echo ""
    echo "  1) eQTL Analysis"
    echo "  2) pQTL Analysis"
    echo "  3) eQTL + pQTL Combined Analysis"
    echo "  4) Traditional MR (selfharm + violence)"
    echo "  5) Exit"
    echo ""
    read -p "Enter your choice [1-5]: " choice
    echo ""

    case $choice in
        1) MR_TYPE="eqtl" ;;
        2) MR_TYPE="pqtl" ;;
        3) MR_TYPE="eqtl_pqtl" ;;
        4) MR_TYPE="traditional" ;;
        5) exit 0 ;;
        *)
            print_msg "$RED" "Invalid choice. Please try again."
            sleep 2
            show_menu
            ;;
    esac
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -t|--type)
                MR_TYPE="$2"
                shift 2
                ;;
            -i|--input)
                INPUT_FILE="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -c|--coloc)
                COLOC_ENABLED="$2"
                shift 2
                ;;
            -p|--pval)
                PVAL_THRESHOLD="$2"
                shift 2
                ;;
            -k|--kb)
                KB_DISTANCE="$2"
                shift 2
                ;;
            -r|--r2)
                R2_THRESHOLD="$2"
                shift 2
                ;;
            -f|--fstat)
                FSTAT_THRESHOLD="$2"
                shift 2
                ;;
            -e|--eaf)
                EAF_THRESHOLD="$2"
                shift 2
                ;;
            -m|--minsnps)
                MIN_SNPS="$2"
                shift 2
                ;;
            -d|--outcome)
                OUTCOME_ID="$2"
                shift 2
                ;;
            -q|--pqtl-dir)
                PQTLDIR="$2"
                shift 2
                ;;
            -h|--help)
                HELP=1
                shift
                ;;
            *)
                print_msg "$RED" "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

# Validate inputs
validate_inputs() {
    local valid_types=("eqtl" "pqtl" "eqtl_pqtl" "traditional")

    # Check MR type
    if [[ -z "$MR_TYPE" ]]; then
        show_menu
    fi

    if [[ ! " ${valid_types[*]} " =~ " ${MR_TYPE} " ]]; then
        print_msg "$RED" "Error: Invalid MR type '$MR_TYPE'"
        print_msg "$RED" "Valid types: eqtl, pqtl, traditional"
        exit 1
    fi

    # Check input file for eQTL/pQTL
    if [[ "$MR_TYPE" != "traditional" && -z "$INPUT_FILE" ]]; then
        print_msg "$RED" "Error: Input file is required for eQTL/pQTL analysis"
        exit 1
    fi

    # Set default output directory if not specified
    if [[ -z "$OUTPUT_DIR" ]]; then
        OUTPUT_DIR="$RESULT_DIR/MR_${MR_TYPE}"
    fi

    # Make output path absolute if relative
    if [[ "$OUTPUT_DIR" != /* ]]; then
        OUTPUT_DIR="$SCRIPT_DIR/$OUTPUT_DIR"
    fi
}

# Create log file
init_log() {
    LOG_FILE="$LOG_DIR/MR_${MR_TYPE}.${TIMESTAMP}.log"

    # Create logs directory if not exists
    mkdir -p "$LOG_DIR"

    # Write log header
    cat > "$LOG_FILE" << EOF
================================================================================
MR分析日志 - Mendelian Randomization Analysis
================================================================================
时间戳: $(date '+%Y-%m-%d %H:%M:%S')
运行ID: $TIMESTAMP
脚本: run_mr.sh

----------------------------------------
输入参数
----------------------------------------
MR类型: $MR_TYPE
输入文件: ${INPUT_FILE:-N/A}
输出目录: $OUTPUT_DIR
共定位分析: $COLOC_ENABLED
P值阈值: $PVAL_THRESHOLD
KB距离: $KB_DISTANCE

----------------------------------------
环境信息
----------------------------------------
工作目录: $SCRIPT_DIR
R脚本目录: $R_SCRIPT_DIR

EOF

    log_message "日志文件已创建: $LOG_FILE"
}

# Run MR analysis
run_mr() {
    local type=$1
    local input=$2
    local output=$3
    local coloc=$4
    local pval=$5
    local kb=$6
    local r2=$7
    local fstat=$8
    local eaf=$9
    local minsnps=${10}
    local outcome=${11}
    local pqtl_dir=${12}

    log_message "开始MR分析..."
    log_message "  类型: $type"
    log_message "  输入: ${input:-N/A}"
    log_message "  输出: $output"
    log_message "  共定位: $coloc"
    log_message "  P值阈值: $pval"
    log_message "  KB距离: $kb"
    log_message "  R2阈值: $r2"
    log_message "  F统计量阈值: $fstat"
    log_message "  EAF阈值: $eaf"
    log_message "  最少SNPs: $minsnps"
    log_message "  结局ID: $outcome"
    if [ "$type" = "pqtl" ]; then
        log_message "  pQTL目录: $pqtl_dir"
    fi

    # Create output directory structure
    mkdir -p "$output/tables"
    mkdir -p "$output/figures"

    # Create hidden timestamp file for tracking
    echo "$TIMESTAMP" > "$output/.$TIMESTAMP"

    log_message "输出目录已创建: $output"

    # Run R script based on type
    case $type in
        eqtl)
            log_message "开始eQTL MR分析..."
            Rscript "$R_SCRIPT_DIR/MR_eqtl.R" \
                --input "$input" \
                --output "$output" \
                --pval "$pval" \
                --kb "$kb" \
                --r2 "$r2" \
                --fstat "$fstat" \
                --eaf "$eaf" \
                --outcome "$outcome" \
                --timestamp "$TIMESTAMP" \
                2>&1 | tee -a "$LOG_FILE"
            ;;
        pqtl)
            log_message "开始pQTL MR分析..."
            Rscript "$R_SCRIPT_DIR/MR_pqtl.R" \
                --input "$input" \
                --output "$output" \
                --pval "$pval" \
                --kb "$kb" \
                --r2 "$r2" \
                --fstat "$fstat" \
                --eaf "$eaf" \
                --outcome "$outcome" \
                --pqtl-dir "$pqtl_dir" \
                --timestamp "$TIMESTAMP" \
                2>&1 | tee -a "$LOG_FILE"
            ;;
        eqtl_pqtl)
            # Create subdirectories for eQTL and pQTL
            local eqtl_output="$output/eqtl"
            local pqtl_output="$output/pqtl"
            mkdir -p "$eqtl_output/tables" "$eqtl_output/figures"
            mkdir -p "$pqtl_output/tables" "$pqtl_output/figures"

            # Run eQTL analysis
            log_message "开始eQTL MR分析..."
            Rscript "$R_SCRIPT_DIR/MR_eqtl.R" \
                --input "$input" \
                --output "$eqtl_output" \
                --pval "$pval" \
                --kb "$kb" \
                --r2 "$r2" \
                --fstat "$fstat" \
                --eaf "$eaf" \
                --outcome "$outcome" \
                --timestamp "$TIMESTAMP" \
                2>&1 | tee -a "$LOG_FILE"

            # Run pQTL analysis
            log_message "开始pQTL MR分析..."
            Rscript "$R_SCRIPT_DIR/MR_pqtl.R" \
                --input "$input" \
                --output "$pqtl_output" \
                --pval "$pval" \
                --kb "$kb" \
                --r2 "$r2" \
                --fstat "$fstat" \
                --eaf "$eaf" \
                --outcome "$outcome" \
                --pqtl-dir "$pqtl_dir" \
                --timestamp "$TIMESTAMP" \
                2>&1 | tee -a "$LOG_FILE"

            # Find common genes
            log_message "寻找eQTL和pQTL共有基因..."
            Rscript "$R_SCRIPT_DIR/find_common_genes.R" \
                --eqtl-dir "$eqtl_output" \
                --pqtl-dir "$pqtl_output" \
                --output "$output" \
                --timestamp "$TIMESTAMP" \
                2>&1 | tee -a "$LOG_FILE"
            ;;
        traditional)
            log_message "开始Traditional MR分析..."
            Rscript "$R_SCRIPT_DIR/MR_traditional.R" \
                --output "$output" \
                --pval "$pval" \
                --kb "$kb" \
                --r2 "$r2" \
                --fstat "$fstat" \
                --minsnps "$minsnps" \
                --timestamp "$TIMESTAMP" \
                2>&1 | tee -a "$LOG_FILE"
            ;;
    esac

    local mr_status=$?

    if [[ $mr_status -eq 0 ]]; then
        log_message "MR分析完成!"
    else
        log_message "MR分析失败 (退出码: $mr_status)"
    fi

    # Run colocalization if enabled
    if [[ "$coloc" == "1" ]]; then
        log_message "开始共定位分析..."
        local coloc_output="$output/../coloc"
        mkdir -p "$coloc_output"

        Rscript "$R_SCRIPT_DIR/coloc_from_mr.R" \
            --input "$output" \
            --output "$coloc_output" \
            --timestamp "$TIMESTAMP" \
            2>&1 | tee -a "$LOG_FILE"

        if [[ $? -eq 0 ]]; then
            log_message "共定位分析完成!"
        else
            log_message "警告: 共定位分析失败"
        fi
    fi

    return $mr_status
}

# Generate config file
generate_config() {
    local output=$1
    local type=$2

    mkdir -p "$output"
    local config_file="$output/${type}.config.ini"

    cat > "$config_file" << EOF
[Analysis]
timestamp = $TIMESTAMP
type = $type

[Input]
gene_list = ${INPUT_FILE:-N/A}

[Output]
base_dir = $output

[Parameters]
pval_threshold = $PVAL_THRESHOLD
kb_distance = $KB_DISTANCE
r2_threshold = $R2_THRESHOLD
fstat_threshold = $FSTAT_THRESHOLD
eaf_threshold = $EAF_THRESHOLD
min_snps = $MIN_SNPS

[Coloc]
enabled = $COLOC_ENABLED
EOF

    log_message "配置文件已生成: $config_file"
}

# Generate R-based Word report
generate_report() {
    local output=$1
    local type=$2

    mkdir -p "$REPORT_DIR"

    local report_script="$R_SCRIPT_DIR/generate_mr_report.R"
    local report_file="$REPORT_DIR/MR_${type}.${TIMESTAMP}.docx"

    # Check if R script exists
    if [[ ! -f "$report_script" ]]; then
        log_message "警告: 报告生成脚本不存在，跳过报告生成"
        return 1
    fi

    # Handle eqtl_pqtl combined analysis
    if [[ "$type" == "eqtl_pqtl" ]]; then
        log_message "开始生成eQTL+pQTL联合报告..."

        # Generate eQTL report
        local eqtl_output="$output/eqtl"
        local eqtl_report="$REPORT_DIR/MR_eqtl.${TIMESTAMP}.docx"
        if [[ -d "$eqtl_output/tables" ]]; then
            log_message "生成eQTL分析报告..."
            Rscript "$report_script" \
                --input "$eqtl_output" \
                --output "$eqtl_report" \
                --type "eqtl" \
                --timestamp "$TIMESTAMP" \
                2>&1 | tee -a "$LOG_FILE"
        fi

        # Generate pQTL report
        local pqtl_output="$output/pqtl"
        local pqtl_report="$REPORT_DIR/MR_pqtl.${TIMESTAMP}.docx"
        if [[ -d "$pqtl_output/tables" ]]; then
            log_message "生成pQTL分析报告..."
            Rscript "$report_script" \
                --input "$pqtl_output" \
                --output "$pqtl_report" \
                --type "pqtl" \
                --timestamp "$TIMESTAMP" \
                2>&1 | tee -a "$LOG_FILE"
        fi

        # Generate combined report
        log_message "生成联合报告..."
        Rscript "$R_SCRIPT_DIR/generate_combined_report.R" \
            --eqtl-dir "$eqtl_output" \
            --pqtl-dir "$pqtl_output" \
            --output "$report_file" \
            --timestamp "$TIMESTAMP" \
            2>&1 | tee -a "$LOG_FILE"

        if [[ $? -eq 0 && -f "$report_file" ]]; then
            log_message "联合报告已生成: $report_file"
        else
            log_message "联合报告生成失败"
        fi
        return 0
    fi

    # Check if result directory exists
    if [[ ! -d "$output/tables" ]]; then
        log_message "警告: 结果目录不存在，跳过报告生成"
        return 1
    fi

    log_message "开始生成Word报告..."

    # Run R script to generate report
    Rscript "$report_script" \
        --input "$output" \
        --output "$report_file" \
        --type "$type" \
        --timestamp "$TIMESTAMP" \
        2>&1 | tee -a "$LOG_FILE"

    if [[ $? -eq 0 && -f "$report_file" ]]; then
        log_message "Word报告已生成: $report_file"
    else
        log_message "报告生成失败，使用简单文本报告"
        generate_simple_report "$output" "$type"
    fi
}

# Generate simple text report (fallback)
generate_simple_report() {
    local output=$1
    local type=$2

    local report_file="$REPORT_DIR/MR_${type}.${TIMESTAMP}.txt"

    cat > "$report_file" << EOF
================================================================================
MR Analysis Report
================================================================================

Analysis Type: $type
Timestamp: $TIMESTAMP
Input File: ${INPUT_FILE:-N/A}
Output Directory: $output

Parameters:
  - P-value threshold: $PVAL_THRESHOLD
  - KB distance: $KB_DISTANCE
  - R2 threshold: $R2_THRESHOLD
  - F-statistic threshold: $FSTAT_THRESHOLD
  - EAF threshold: $EAF_THRESHOLD
  - Colocalization: $COLOC_ENABLED

Result Files:
  - tables/: Analysis results
  - figures/: Visualization plots

Log File: logs/MR_${type}.${TIMESTAMP}.log
Config: ${type}.config.ini

================================================================================
EOF

    # 读取并汇总结果
    if [[ -f "$output/00.Complete_MR_Results.csv" ]]; then
        echo "" >> "$report_file"
        echo "Summary:" >> "$report_file"
        echo "----------" >> "$report_file"

        total_genes=$(wc -l < "$output/00.Complete_MR_Results.csv")
        total_genes=$((total_genes - 1))

        echo "Total genes analyzed: $total_genes" >> "$report_file"

        # 统计显著基因数（ivw_significant == "YES"）
        sig_count=$(awk -F',' 'NR>1 && /YES/' "$output/06.Three_Filter_Results.csv" | wc -l)
        echo "Significant genes (IVW p<0.05): $sig_count" >> "$report_file"
    fi

    log_message "简单报告已生成: $report_file"
}

# Main function
main() {
    # Parse arguments
    parse_args "$@"

    # Show help if requested
    if [[ $HELP -eq 1 ]]; then
        show_help
        exit 0
    fi

    # Validate inputs
    validate_inputs

    # Initialize log file
    init_log

    # Print header
    echo ""
    print_msg "$BLUE" "=============================================="
    print_msg "$BLUE" "       MR Analysis Pipeline"
    print_msg "$BLUE" "  Timestamp: $TIMESTAMP"
    print_msg "$BLUE" "=============================================="
    echo ""

    # Generate config file
    generate_config "$OUTPUT_DIR" "$MR_TYPE"

    # Run MR analysis
    run_mr "$MR_TYPE" "$INPUT_FILE" "$OUTPUT_DIR" "$COLOC_ENABLED" "$PVAL_THRESHOLD" "$KB_DISTANCE" "$R2_THRESHOLD" "$FSTAT_THRESHOLD" "$EAF_THRESHOLD" "$MIN_SNPS" "$OUTCOME_ID" "$PQTLDIR"
    local status=$?

    # Generate report
    generate_report "$OUTPUT_DIR" "$MR_TYPE"

    # Print summary
    echo ""
    if [[ $status -eq 0 ]]; then
        print_msg "$GREEN" "=============================================="
        print_msg "$GREEN" "       分析完成!"
        print_msg "$GREEN" "  运行ID: $TIMESTAMP"
        print_msg "$GREEN" "  结果: $OUTPUT_DIR"
        print_msg "$GREEN" "  日志: $LOG_FILE"
        print_msg "$GREEN" "=============================================="
    else
        print_msg "$RED" "=============================================="
        print_msg "$RED" "       分析失败 (退出码: $status)"
        print_msg "$RED" "  日志: $LOG_FILE"
        print_msg "$RED" "=============================================="
    fi
}

# Run main function
main "$@"
