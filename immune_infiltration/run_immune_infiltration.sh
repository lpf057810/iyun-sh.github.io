#!/bin/bash
# run_immune_infiltration.sh - 免疫浸润分析运行脚本
# 通过 sudo 以目标用户身份运行 Python 脚本

TARGET_USER=$(stat -c "%U" "$0")
SCRIPT_DIR="/media/desk16/share/secure/immune_infiltration"
SCRIPT_DIR_ABS="/media/desk16/share/secure/immune_infiltration"

# 生成时间戳
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# 显示帮助信息
show_help() {
    cat << 'EOF'
用法: run_immune_infiltration.sh [选项]

必需选项（二选一）:
  --config CONFIG    配置文件路径(INI格式)
  -i, --interactive  交互式选择要运行的方法（会提示输入文件）

可选参数:
  -o, --output DIR   输出目录(覆盖配置文件中的base_dir)
  --methods METHODS  要执行的方法编号(如: 1,7 或 1-5)
  --workers N        最大并行任务数
  --force            强制重新计算,忽略缓存
  --dry-run          仅显示将要执行的命令,不实际运行
  --verbose          详细输出模式
  --rm               删除输出目录
  -h, --help         显示此帮助信息

可用方法:
 编号  方法        说明
   1    cibersort   CIBERSORT免疫细胞反卷积
   2    epic        EPIC免疫浸润分析
   3    estimate    ESTIMATE免疫评分
   4    ips         IPS免疫表性评分
   5    mcpcounter  MCPcounter免疫细胞计数
   6    ssgsea      ssGSEA基因集富集分析(需要signature文件)
   7    timer       TIMER免疫浸润(需要tissue参数)
   8    xcell       xCell细胞富集评分

示例:
  # 交互式运行（推荐）
  run_immune_infiltration.sh --interactive

  # 使用配置文件运行
  run_immune_infiltration.sh --config config/config.ini -o /path/to/output

  # 仅运行特定方法
  run_immune_infiltration.sh --config config/config.ini -o /path/to/output --methods 1,7

  # 预览将要执行的命令
  run_immune_infiltration.sh --config config/config.ini -o /path/to/output --dry-run

配置文件示例: config/config.example.ini
EOF
}

# 检查是否有参数
if [[ $# -eq 0 ]]; then
    show_help
    exit 0
fi

# 处理 -h/--help
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# 保存原始参数
ORIGINAL_ARGS=("$@")

# 解析参数
OUTPUT_DIR=""
DO_RM=false
INTERACTIVE=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --rm)
            DO_RM=true
            shift
            ;;
        -i|--interactive)
            INTERACTIVE=true
            shift
            ;;
        *)
            shift
            ;;
    esac
done

# 如果不是目标用户，通过 sudo 重新运行
if [ "$(whoami)" != "$TARGET_USER" ]; then
    exec sudo -u "$TARGET_USER" "$SCRIPT_DIR/run_immune_infiltration.sh" "${ORIGINAL_ARGS[@]}"
fi

# 以下以目标用户身份执行
cd "$SCRIPT_DIR"

# 创建日志目录
mkdir -p "$SCRIPT_DIR/logs"

# 创建报告目录
mkdir -p "$SCRIPT_DIR/report"

# 运行日志文件
LOG_FILE="$SCRIPT_DIR/logs/immune_infiltration.${TIMESTAMP}.log"

# 写入日志头部
cat > "$LOG_FILE" << EOF
================================================================================
免疫浸润分析日志 - Immune Infiltration Analysis
================================================================================
时间戳: $(date '+%Y-%m-%d %H:%M:%S')
运行ID: $TIMESTAMP
脚本: run_immune_infiltration.sh

----------------------------------------
环境信息
----------------------------------------
工作目录: $SCRIPT_DIR
R脚本目录: $SCRIPT_DIR/scripts

EOF

echo "开始免疫浸润分析..." | tee -a "$LOG_FILE"

# 运行分析
python3 scripts/run_immune_infiltration.py "${ORIGINAL_ARGS[@]}" 2>&1 | tee -a "$LOG_FILE"

# 生成图表
python3 scripts/generate_plots.py 2>&1 | tee -a "$LOG_FILE"

# 生成报告 (使用R脚本)
echo "开始生成Word报告..." | tee -a "$LOG_FILE"

# 获取输出目录
OUTPUT_DIR=$(grep -E "^base_dir\s*=" "$SCRIPT_DIR/config/config.ini" | cut -d'=' -f2 | tr -d ' ')
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="$SCRIPT_DIR/result"
fi

REPORT_FILE="$SCRIPT_DIR/report/immune_infiltration.${TIMESTAMP}.docx"

# 调用R脚本生成报告
Rscript "$SCRIPT_DIR/scripts/generate_immune_report.R" \
    --input "$OUTPUT_DIR" \
    --output "$REPORT_FILE" \
    --config "$SCRIPT_DIR/config/config.ini" \
    --timestamp "$TIMESTAMP" \
    2>&1 | tee -a "$LOG_FILE"

if [ $? -eq 0 ]; then
    echo "Word报告已生成: $REPORT_FILE" | tee -a "$LOG_FILE"
else
    echo "R报告生成失败，使用Python脚本..." | tee -a "$LOG_FILE"
    python3 scripts/generate_docx.py 2>&1 | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "==============================================" | tee -a "$LOG_FILE"
echo "分析完成!" | tee -a "$LOG_FILE"
echo "运行ID: $TIMESTAMP" | tee -a "$LOG_FILE"
echo "结果: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "日志: $LOG_FILE" | tee -a "$LOG_FILE"
echo "报告: $REPORT_FILE" | tee -a "$LOG_FILE"
echo "==============================================" | tee -a "$LOG_FILE"
