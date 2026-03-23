#!/bin/bash
# ==============================================================================
# Prognosis/run_prognosis.sh
# 预后分析一键化运行脚本
# 功能：Cox回归 + Lasso特征选择 + 风险模型 + KM生存分析
# ==============================================================================

set -e

TARGET_USER=$(stat -c "%U" "$0")
SCRIPT_DIR="/media/desk16/share/secure/Prognosis"
FULL_ARGS=("$@")

WORK_DIR=""
EXPR_FILE=""
GROUP_FILE=""
SURVIVAL_FILE=""
GENE_FILE=""
INPUT_FILE=""
OUTPUT_DIR="results"
OUTPUT_PATH=""
DO_RM=false

CONFIG_FILE="$SCRIPT_DIR/scripts/config/Prognosis.config.ini"
if [ -f "$CONFIG_FILE" ]; then
    source <(grep -E '^[a-zA-Z_]+=' "$CONFIG_FILE" | sed 's/=/="/' | sed 's/$/"/')
    echo "已加载配置文件: $CONFIG_FILE"
fi

show_usage() {
    cat << 'EOF'
用法: run_prognosis.sh [参数]

【合并格式（推荐）】 仅需一个文件
  -i, --input FILE       合并格式文件（含OS.time/OS/基因表达，df.csv）
  -o, --output-dir DIR   输出目录

【分离格式（兼容）】 需要四个文件
  -w, --work-dir DIR         工作目录
  -e, --expr-file FILE       表达矩阵文件
  -g, --group-file FILE      分组文件
  -s, --survival-file FILE   生存数据文件
  -G, --gene-file FILE       基因列表文件
  -o, --output-dir DIR       输出目录

【其他】
  --rm                     删除输出目录（仅创建者可执行）
  -h, --help               显示帮助信息

【示例】
  # 合并格式（一行搞定）
  ./run_prognosis.sh -i df.csv -o results

  # 分离格式
  ./run_prognosis.sh -w /path/to/project \
    -e expression.csv -g group.csv -s survival.csv -G genes.csv -o results
EOF
}

if [ "$(whoami)" != "$TARGET_USER" ]; then
    exec sudo -u "$TARGET_USER" "$SCRIPT_DIR/run_prognosis.sh" "${FULL_ARGS[@]}"
fi

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -w|--work-dir)       WORK_DIR="$2"; shift 2 ;;
            -e|--expr-file)      EXPR_FILE="$2"; shift 2 ;;
            -g|--group-file)     GROUP_FILE="$2"; shift 2 ;;
            -s|--survival-file)  SURVIVAL_FILE="$2"; shift 2 ;;
            -G|--gene-file)      GENE_FILE="$2"; shift 2 ;;
            -i|--input)          INPUT_FILE="$2"; shift 2 ;;
            -o|--output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
            --rm)                DO_RM=true; shift ;;
            -h|--help)           show_usage; exit 0 ;;
            *)  echo "错误: 未知参数 '$1'"; show_usage; exit 1 ;;
        esac
    done
}

validate_args() {
    # 合并格式：只需 --input 和 --output-dir
    if [[ -n "$INPUT_FILE" ]]; then
        if [[ ! -f "$INPUT_FILE" ]]; then
            echo "错误: 文件不存在: $INPUT_FILE"
            exit 1
        fi
        # 合并格式不需要四个文件
        return 0
    fi

    # 分离格式：需要工作目录和四个文件
    if [[ -z "$WORK_DIR" ]]; then
        echo "错误: 缺少必需参数 -w/--work-dir（或使用 -i 指定合并格式文件）"
        show_usage; exit 1
    fi
    if [[ -z "$EXPR_FILE" ]]; then
        echo "错误: 缺少必需参数 -e/--expr-file（或使用 -i 指定合并格式文件）"
        show_usage; exit 1
    fi
    if [[ -z "$GROUP_FILE" ]]; then
        echo "错误: 缺少必需参数 -g/--group-file（或使用 -i 指定合并格式文件）"
        show_usage; exit 1
    fi
    if [[ -z "$SURVIVAL_FILE" ]]; then
        echo "错误: 缺少必需参数 -s/--survival-file（或使用 -i 指定合并格式文件）"
        show_usage; exit 1
    fi
    if [[ -z "$GENE_FILE" ]]; then
        echo "错误: 缺少必需参数 -G/--gene-file（或使用 -i 指定合并格式文件）"
        show_usage; exit 1
    fi
}

# 处理 --rm 模式
handle_rm() {
    local delete_dir="$1"
    if [ ! -d "$delete_dir" ]; then
        echo "错误: 输出目录不存在: $delete_dir"
        exit 1
    fi
    if [ -f "$delete_dir/.owner" ]; then
        OWNER=$(cat "$delete_dir/.owner")
        CALLER=${SUDO_USER:-$USER}
        if [ "$OWNER" != "$CALLER" ]; then
            echo "错误: 只能删除自己创建的目录（创建者: $OWNER，当前用户: $CALLER）"
            exit 1
        fi
    fi
    rm -rf "$delete_dir"
    echo "已删除: $delete_dir"
    exit 0
}

# 确定输出路径
resolve_output_path() {
    if [[ -n "$INPUT_FILE" ]]; then
        # 合并格式：输出到 --output-dir 指定路径
        OUTPUT_PATH="$OUTPUT_DIR"
    else
        # 分离格式：默认到 result/，否则到 OUTPUT_DIR
        if [[ -z "$OUTPUT_DIR" ]] || [[ "$OUTPUT_DIR" == "results" ]]; then
            OUTPUT_PATH="$SCRIPT_DIR/result"
        else
            OUTPUT_PATH="$OUTPUT_DIR"
        fi
    fi
}

main() {
    parse_args "$@"
    validate_args

    resolve_output_path

    # --rm 模式
    if [ "$DO_RM" = true ]; then
        handle_rm "$OUTPUT_PATH"
    fi

    # 创建输出目录
    mkdir -p "$OUTPUT_PATH/Cox"
    mkdir -p "$OUTPUT_PATH/RiskModel"

    echo "============================================================"
    echo "        预后分析流程 (Cox + Risk Model)"
    echo "============================================================"
    echo ""

    # ---------- 合并格式 ----------
    if [[ -n "$INPUT_FILE" ]]; then
        echo "输入模式: 合并格式（自动检测）"
        echo "输入文件: $INPUT_FILE"
        echo "输出目录: $OUTPUT_PATH"
        echo ""

        echo "[1/2] 运行 Cox 回归分析..."
        Rscript "$SCRIPT_DIR/scripts/Cox.R" \
            --expr "$INPUT_FILE" \
            --output "$OUTPUT_PATH/Cox"

        echo "[2/2] 运行风险模型分析..."
        Rscript "$SCRIPT_DIR/scripts/RiskModel.R" \
            --expr "$INPUT_FILE" \
            --lasso-coef "$OUTPUT_PATH/Cox/03.Lasso_Coefficients.csv" \
            --output "$OUTPUT_PATH/RiskModel"

        echo "[3/3] 生成 Word 报告..."
        Rscript "$SCRIPT_DIR/scripts/generate_prognosis_report.R" \
            --cox-dir "$OUTPUT_PATH/Cox" \
            --risk-dir "$OUTPUT_PATH/RiskModel" \
            --output "$SCRIPT_DIR/report"

    # ---------- 分离格式 ----------
    else
        cd "$WORK_DIR"
        echo "输入模式: 分离格式"
        echo "工作目录: $(pwd)"
        echo "表达矩阵: $EXPR_FILE"
        echo "分组文件: $GROUP_FILE"
        echo "生存数据: $SURVIVAL_FILE"
        echo "基因列表: $GENE_FILE"
        echo "输出目录: $OUTPUT_PATH"
        echo ""

        echo "[1/2] 运行 Cox 回归分析..."
        Rscript "$SCRIPT_DIR/scripts/Cox.R" \
            --expr "$EXPR_FILE" \
            --group "$GROUP_FILE" \
            --survival "$SURVIVAL_FILE" \
            --genes "$GENE_FILE" \
            --output "$OUTPUT_PATH/Cox"

        echo "[2/2] 运行风险模型分析..."
        Rscript "$SCRIPT_DIR/scripts/RiskModel.R" \
            --expr "$EXPR_FILE" \
            --group "$GROUP_FILE" \
            --survival "$SURVIVAL_FILE" \
            --lasso-coef "$OUTPUT_PATH/Cox/03.Lasso_Coefficients.csv" \
            --output "$OUTPUT_PATH/RiskModel" \
            --report

        echo "[3/3] 生成 Word 报告..."
        Rscript "$SCRIPT_DIR/scripts/generate_prognosis_report.R" \
            --cox-dir "$OUTPUT_PATH/Cox" \
            --risk-dir "$OUTPUT_PATH/RiskModel" \
            --output "$SCRIPT_DIR/report"
    fi

    # 设置权限
    if [ -n "$OUTPUT_PATH" ] && [ -d "$OUTPUT_PATH" ]; then
        CALLER=${SUDO_USER:-$USER}
        echo "$CALLER" > "$OUTPUT_PATH/.owner"
        chmod -R 777 "$OUTPUT_PATH"
        chmod +t "$OUTPUT_PATH"
        chmod 444 "$OUTPUT_PATH/.owner"
    fi

    echo ""
    echo "============================================================"
    echo "                      分析完成"
    echo "============================================================"
    echo "输出目录: $OUTPUT_PATH"
    echo ""
}

main "$@"
