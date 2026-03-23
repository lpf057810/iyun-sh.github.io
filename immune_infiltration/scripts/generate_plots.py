#!/usr/bin/env python3
"""
免疫浸润分析 - 图表生成脚本
从 R 生成的 CSV 结果生成可视化图表
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300

def load_data(expr_file, group_file):
    """加载表达矩阵和分组数据"""
    expr = pd.read_csv(expr_file, index_col=0)
    group = pd.read_csv(group_file)
    return expr, group

def plot_stacked_bar(data, group, output_dir, method_name="CIBERSORT"):
    """绘制堆叠柱状图"""
    # 过滤非数值列
    exclude_cols = ['ID', 'P-value', 'Correlation', 'RMSE', 'sample', 'TumorPurity']
    cell_cols = []
    for c in data.columns:
        if c not in exclude_cols:
            # 检查是否为数值列
            try:
                pd.to_numeric(data[c], errors='coerce')
                cell_cols.append(c)
            except:
                pass

    if not cell_cols:
        print(f"  [警告] {method_name}: 无有效数值列")
        return

    # 合并分组信息
    data_plot = data.copy()
    # 确保合并时使用正确的列名
    group_cols = [c for c in group.columns if c.lower() in ['sample', 'id', 'sample_id']]
    if group_cols:
        merge_col = group_cols[0]
    else:
        merge_col = group.columns[0]
    data_plot = data_plot.merge(group, left_on='ID', right_on=merge_col, how='left')

    # 移除没有分组的行
    data_plot = data_plot.dropna(subset=['group'])
    if len(data_plot) == 0:
        print(f"  [警告] {method_name}: 无有效分组数据")
        return

    # 按分组排序
    data_plot = data_plot.sort_values('group')

    # 准备数据
    plot_data = data_plot[cell_cols].fillna(0)
    samples = data_plot['ID'].values
    groups = data_plot['group'].values

    # 获取唯一分组
    unique_groups = list(set(groups))
    fig, ax = plt.subplots(figsize=(14, 6))

    # 堆叠柱状图
    bottom = np.zeros(len(samples))
    colors = sns.color_palette("husl", n_colors=len(cell_cols))

    for i, col in enumerate(cell_cols):
        ax.bar(range(len(samples)), plot_data[col], bottom=bottom, label=col, color=colors[i], width=0.8)
        bottom += plot_data[col].values

    # 设置标签
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_title(f'Immune Cell Composition - {method_name}', fontsize=14)

    # 添加分组标记
    unique_groups = list(set(groups))
    group_positions = {}
    for g in unique_groups:
        positions = [j for j, gx in enumerate(groups) if gx == g]
        group_positions[g] = (min(positions), max(positions))

    for g, (start, end) in group_positions.items():
        mid = (start + end) / 2
        ax.axvline(x=end + 0.5, color='gray', linestyle='--', alpha=0.5)
        ax.text(mid, ax.get_ylim()[1] * 1.02, g, ha='center', fontsize=10, fontweight='bold')

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, ncol=1)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f'01.{method_name.lower()}_stacked_bar.png')
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  [保存] {output_path}")

def plot_boxplot(data, group, output_dir, method_name="CIBERSORT"):
    """绘制箱线图"""
    # 过滤非数值列
    exclude_cols = ['ID', 'P-value', 'Correlation', 'RMSE', 'sample']
    cell_cols = []
    for c in data.columns:
        if c not in exclude_cols:
            # 检查是否为数值列
            try:
                pd.to_numeric(data[c], errors='coerce')
                cell_cols.append(c)
            except:
                pass

    if not cell_cols:
        print(f"  [警告] {method_name}: 无有效数值列")
        return

    # 转换数据格式
    data_plot = data.copy()
    # 确保合并时使用正确的列名
    group_cols = [c for c in group.columns if c.lower() in ['sample', 'id', 'sample_id']]
    if group_cols:
        merge_col = group_cols[0]
    else:
        merge_col = group.columns[0]
    data_plot = data_plot.merge(group, left_on='ID', right_on=merge_col, how='left')

    # 转换为长格式
    try:
        plot_data = data_plot.melt(id_vars=['ID', 'group'], value_vars=cell_cols,
                                    var_name='Cell_Type', value_name='Score')
    except Exception as e:
        print(f"  [警告] {method_name}: 数据转换失败 - {e}")
        return

    # 创建图形
    n_cells = len(cell_cols)
    n_cols = 4
    n_rows = (n_cells + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, n_rows * 3))

    # 处理 axes 的不同形状
    if n_rows == 1:
        axes = np.array([axes]) if n_cells == 1 else np.array(axes).flatten()
    else:
        axes = axes.flatten()

    for i, cell in enumerate(cell_cols):
        ax = axes[i]
        cell_data = plot_data[plot_data['Cell_Type'] == cell]

        # 绘制箱线图
        sns.boxplot(data=cell_data, x='group', y='Score', ax=ax, palette=['#4ECDC4', '#FF6B6B'])
        sns.stripplot(data=cell_data, x='group', y='Score', ax=ax, color='black', alpha=0.3, size=3)

        # 统计检验
        groups = cell_data['group'].unique()
        if len(groups) == 2:
            g1 = cell_data[cell_data['group'] == groups[0]]['Score']
            g2 = cell_data[cell_data['group'] == groups[1]]['Score']
            if len(g1) > 0 and len(g2) > 0:
                stat, pval = stats.mannwhitneyu(g1, g2, alternative='two-sided')
                sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
                ax.text(0.5, 0.95, f'p={pval:.3f} {sig}', transform=ax.transAxes,
                       ha='center', fontsize=8)

        ax.set_title(cell, fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('Score' if i % n_cols == 0 else '')

    # 隐藏空白子图
    for i in range(len(cell_cols), len(axes)):
        axes[i].axis('off')

    plt.suptitle(f'Immune Cell Comparison - {method_name}', fontsize=14, y=1.02)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f'02.{method_name.lower()}_boxplot.png')
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  [保存] {output_path}")

def main():
    # 路径
    base_dir = "/media/desk16/share/secure/immune_infiltration/immune_infiltration_results"

    # 方法列表
    methods = [
        ('cibersort_output', 'CIBERSORT'),
        ('epic_output', 'EPIC'),
        ('estimate_output', 'ESTIMATE')
    ]

    # 分组文件
    group_file = "/media/desk16/share/secure/immune_infiltration/GSE126124.group.csv"

    print("=" * 50)
    print("生成免疫浸润图表")
    print("=" * 50)

    for method_dir, method_name in methods:
        print(f"\n处理 {method_name}...")

        # 检查数据文件
        data_file = os.path.join(base_dir, method_dir, '01.ciber_res.csv')
        if not os.path.exists(data_file):
            # 尝试其他文件名模式
            possible_files = [
                os.path.join(base_dir, method_dir, '01.ciber_res.csv'),
                os.path.join(base_dir, method_dir, '01.estimate_res.csv'),
                os.path.join(base_dir, method_dir, '01.epic_res.csv'),
            ]
            data_file = None
            for f in possible_files:
                if os.path.exists(f):
                    data_file = f
                    break

        if not data_file or not os.path.exists(data_file):
            print(f"  [跳过] {method_name}: 无结果文件")
            continue

        output_dir = os.path.join(base_dir, method_dir)
        os.makedirs(output_dir, exist_ok=True)

        try:
            # 加载数据
            data = pd.read_csv(data_file)

            # 加载分组
            group = pd.read_csv(group_file)

            # 生成图表
            plot_stacked_bar(data, group, output_dir, method_name)
            plot_boxplot(data, group, output_dir, method_name)

        except Exception as e:
            print(f"  [错误] {method_name}: {e}")

    print("\n完成!")

if __name__ == '__main__':
    main()
