#!/bin/bash
# 启动免疫浸润分析 Dashboard

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# 自动找一个可用端口
for port in 3838 5000 6000 7000 8000 8765 9000; do
    if ! netstat -tuln 2>/dev/null | grep -q ":$port " && \
       ! ss -tuln 2>/dev/null | grep -q ":$port "; then
        echo "🚀 启动 Dashboard，端口: $port"
        Rscript -e "shiny::runApp('.', port = $port, launch.browser = TRUE)"
        exit 0
    fi
done

echo "⚠️ 没有找到可用端口，请手动指定"
