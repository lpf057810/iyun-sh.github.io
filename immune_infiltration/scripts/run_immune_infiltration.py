#!/usr/bin/env python3
import os
import sys
import argparse
import configparser
import subprocess
import shutil
import logging
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple

# 方法编号映射
METHOD_MAP = {
    1: {"name": "cibersort", "script": "cibersort.R", "has_genes": True, "has_signature": False},
    2: {"name": "epic", "script": "epic.R", "has_genes": True, "has_signature": False},
    3: {"name": "estimate", "script": "estimate.R", "has_genes": False, "has_signature": False},
    4: {"name": "ips", "script": "ips.R", "has_genes": False, "has_signature": False},
    5: {"name": "mcpcounter", "script": "mcpcounter.R", "has_genes": True, "has_signature": False},
    6: {"name": "ssgsea", "script": "ssgsea.R", "has_genes": True, "has_signature": True},
    7: {"name": "timer", "script": "timer.R", "has_genes": True, "has_signature": False},
    8: {"name": "xcell", "script": "xcell.R", "has_genes": True, "has_signature": False}
}

# 支持genes参数的方法
GENES_SUPPORT_METHODS = {
    1: "cibersort",
    2: "epic",
    5: "mcpcounter",
    6: "ssgsea",
    7: "timer",
    8: "xcell"
}


def get_timestamp():
    """获取当前时间戳，格式：YYYYMMDD_HHMMSS"""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def setup_logging(log_file: str, verbose: bool = False):
    """设置日志"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )


def write_log_header(config: configparser.ConfigParser, methods: List[int], timestamp: str):
    """写入标准日志头部"""
    import platform
    import subprocess

    # 获取 Python 和 R 版本
    py_version = platform.python_version()
    try:
        r_version = subprocess.run(['Rscript', '--version'], capture_output=True, text=True)
        r_ver = r_version.stderr.split()[2] if r_version.returncode == 0 else "Unknown"
    except:
        r_ver = "Unknown"

    # 获取方法名称
    method_names = [METHOD_MAP[m]['name'] for m in methods if m in METHOD_MAP]

    header = f"""
================================================================================
免疫浸润分析日志 - Immune Infiltration Analysis
================================================================================
时间戳: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
运行ID: {timestamp}
脚本: run_immune_infiltration.py v1.0
运行环境: Python {py_version}, R {r_ver}

----------------------------------------
输入参数
----------------------------------------
方法: {', '.join(method_names)}
CIBERSORT 置换次数: {config.get('Method1.CIBERSORT', 'perm', fallback='100')}
TIMER 组织类型: {config.get('Method7.TIMER', 'tissue', fallback='LUAD')}

----------------------------------------
输入文件
----------------------------------------
表达矩阵: {config['Input']['expression']}
分组文件: {config['Input']['group']}
基因文件: {config['Input'].get('genes', 'Not specified')}
Signature文件: {config['Input'].get('signature', 'Not specified')}

----------------------------------------
输出目录
----------------------------------------
结果目录: {config['Output']['base_dir']}
qs2 缓存: <结果目录>/qs2/ （由程序自动设置，与配置文件中 [Cache] cache_dir 无关）

================================================================================
"""
    logging.info(header)


def parse_methods(methods_str: str) -> List[int]:
    """解析方法字符串,支持逗号分隔和范围"""
    methods = []
    for part in methods_str.split(','):
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            methods.extend(range(start, end + 1))
        elif part:
            methods.append(int(part))
    return methods


def load_config(config_path: str) -> configparser.ConfigParser:
    """加载配置文件"""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"配置文件不存在: {config_path}")
    
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


def validate_config(config: configparser.ConfigParser) -> Tuple[List[str], List[str]]:
    """验证配置文件中的所有必需参数"""
    errors = []
    warnings = []
    
    # 获取脚本所在目录作为基础路径
    script_dir = Path(__file__).parent
    
    # 1. 基础文件验证
    expression = config['Input']['expression']
    group = config['Input']['group']
    
    if not os.path.exists(expression):
        # 尝试相对路径
        rel_path = script_dir / expression
        if rel_path.exists():
            expression = str(rel_path)
        else:
            errors.append(f"表达矩阵文件不存在: {expression}")
    
    if not os.path.exists(group):
        # 尝试相对路径
        rel_path = script_dir / group
        if rel_path.exists():
            group = str(rel_path)
        else:
            errors.append(f"分组信息文件不存在: {group}")
    
    # 2. 可选文件验证
    genes = config['Input']['genes']
    if genes and genes.strip():
        if not os.path.exists(genes):
            # 尝试相对路径
            rel_path = script_dir / genes
            if rel_path.exists():
                config['Input']['genes'] = str(rel_path)
            else:
                errors.append(f"基因列表文件不存在: {genes}")
    
    signature = config['Input']['signature']
    if signature and signature.strip():
        if not os.path.exists(signature):
            # 尝试相对路径
            rel_path = script_dir / signature
            if rel_path.exists():
                config['Input']['signature'] = str(rel_path)
            else:
                errors.append(f"基因签名文件不存在: {signature}")
    
    # 3. 方法特定验证
    try:
        methods = parse_methods(config['Methods']['methods'])
    except ValueError as e:
        errors.append(f"方法列表格式错误: {e}")
        return errors, warnings
    
    # 检查方法编号是否有效
    for m in methods:
        if m not in METHOD_MAP:
            warnings.append(f"未知的方法编号: {m} (有效范围: 1-8)")
    
    # ssgsea需要signature文件
    if 6 in methods:
        if not signature or not signature.strip():
            errors.append(
                "方法6(ssgsea)需要signature文件，但配置文件中未指定！\n"
                "请在[Input]节中配置signature文件，或在方法列表中移除ssgsea(6)"
            )
    
    return errors, warnings


def show_genes_tips(config: configparser.ConfigParser):
    """显示关于genes参数的提示"""
    genes = config['Input']['genes']
    
    if not genes or not genes.strip():
        try:
            methods = parse_methods(config['Methods']['methods'])
        except ValueError:
            return
        
        supported = []
        for m in methods:
            if m in GENES_SUPPORT_METHODS:
                supported.append(GENES_SUPPORT_METHODS[m])
        
        if supported:
            print("\n⚠️  提示: 基因相关性分析")
            print(f"  以下方法支持基因相关性分析: {', '.join(supported)}")
            print(f"  如需分析基因与免疫细胞的相关性，请在配置文件[Input]中指定genes文件")


def validate_r_scripts(script_dir: str, methods: List[int]) -> Tuple[List[str], List[str]]:
    """验证R脚本存在性"""
    errors = []

    for m in methods:
        if m in METHOD_MAP:
            script_name = METHOD_MAP[m]['script']
            script_path = os.path.join(script_dir, script_name)
            if not os.path.exists(script_path):
                errors.append(f"R脚本不存在: {script_path}")

    return errors, []


def validate_r_script_params(script_path: str, method_info: Dict) -> List[str]:
    """验证R脚本支持的参数与所需参数是否匹配

    Args:
        script_path: R脚本路径
        method_info: 方法信息字典

    Returns:
        错误列表（空表示验证通过）
    """
    errors = []

    try:
        with open(script_path, 'r') as f:
            content = f.read()

        # 提取所有 make_option 定义的参数
        import re
        pattern = r'make_option\s*\(\s*c\s*\(\s*"([^"]+)"'
        found_params = re.findall(pattern, content)

        # 扁平化参数列表（处理长选项如 "--force"）
        flat_params = []
        for p in found_params:
            flat_params.extend(p.split(','))

        flat_params = [p.strip().strip('"').strip("'") for p in flat_params]

        # 必需的标准参数
        required_params = ['-e', '-g', '-o']

        for req_param in required_params:
            if req_param not in flat_params:
                errors.append(
                    f"{method_info['name']}: 脚本缺少必需参数 '{req_param}'"
                )

        # 检查方法特定参数
        if method_info.get('has_signature', False):
            if '-s' not in flat_params and '--signature' not in flat_params:
                errors.append(
                    f"{method_info['name']}: 脚本缺少signature参数"
                )

        # 检查是否使用了已弃用的方法
        if 'sys.frame' in content and 'dirname(sys.frame' in content:
            errors.append(
                f"{method_info['name']}: 脚本使用了不稳定的 sys.frame 方法获取路径"
            )

    except Exception as e:
        errors.append(f"验证参数时出错: {str(e)}")

    return errors


def validate_all_r_params(script_dir: str, methods: List[int]) -> Tuple[List[str], List[str]]:
    """验证所有R脚本的参数一致性"""
    errors = []
    warnings = []

    for m in methods:
        if m in METHOD_MAP:
            method_info = METHOD_MAP[m]
            script_path = os.path.join(script_dir, method_info['script'])

            if os.path.exists(script_path):
                param_errors = validate_r_script_params(script_path, method_info)
                errors.extend(param_errors)

    return errors, warnings


# ==============================================================================
# 质控模块
# ==============================================================================
import csv
import statistics

def check_expression_matrix(file_path: str) -> Tuple[List[str], Dict]:
    """检查表达矩阵质量（纯Python实现，不依赖pandas）"""
    errors = []
    warnings = []
    stats = {}

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            rows = list(reader)

        if len(rows) < 2:
            errors.append("表达矩阵为空或只有表头")
            return errors, stats

        # 检查第一列名称
        first_col = rows[0][0]
        if first_col.upper() != 'SYMBOL':
            # 尝试自动修复：如果第一列是空字符串（索引名），给出明确的修复建议
            if first_col == '' or first_col == '"':
                errors.append(f"表达矩阵格式错误：第一列名为空")
                errors.append(f"  修复方法：将第一列重命名为 'SYMBOL'")
                errors.append(f"  Python命令: df.index.name = 'SYMBOL'; df.to_csv('output.csv')")
            else:
                errors.append(f"表达矩阵第一列名称应为 'SYMBOL'，当前为 '{first_col}'")
                errors.append(f"  修复方法：将第一列重命名为 'SYMBOL'")

        # 获取样本名（列名，排除第一列SYMBOL）
        sample_names = rows[0][1:]
        n_genes = len(rows) - 1  # 减去表头
        n_samples = len(sample_names)
        stats['n_genes'] = n_genes
        stats['n_samples'] = n_samples

        if n_genes < 100:
            warnings.append(f"基因数量过少 ({n_genes})，可能影响分析结果")

        # 检查数据
        missing_count = 0
        negative_count = 0
        gene_means = []

        for row in rows[1:]:
            if len(row) < n_samples + 1:
                continue

            values = []
            for val in row[1:n_samples+1]:
                try:
                    v = float(val)
                    if v is None or (isinstance(v, float) and str(v) == 'nan'):
                        missing_count += 1
                    else:
                        values.append(v)
                        if v < 0:
                            negative_count += 1
                except ValueError:
                    missing_count += 1

            if values:
                gene_means.append(statistics.mean(values))

        total_cells = n_genes * n_samples
        missing_pct = (missing_count / total_cells * 100) if total_cells > 0 else 0
        stats['missing_count'] = missing_count
        stats['missing_pct'] = round(missing_pct, 2)

        if missing_pct > 5:
            warnings.append(f"表达矩阵缺失值比例较高 ({missing_pct:.1f}%)")
        elif missing_pct > 20:
            errors.append(f"表达矩阵缺失值过多 ({missing_pct:.1f}%)，请检查数据")

        if negative_count > 0:
            warnings.append(f"发现 {negative_count} 个负值，TPM 数据不应有负值")

        # 检查常量基因
        constant_genes = 0
        for row in rows[1:]:
            if len(row) < n_samples + 1:
                continue
            values = []
            for val in row[1:n_samples+1]:
                try:
                    v = float(val)
                    if v is not None and not (isinstance(v, float) and str(v) == 'nan'):
                        values.append(v)
                except:
                    pass
            if len(values) > 1 and statistics.stdev(values) == 0:
                constant_genes += 1

        stats['constant_genes'] = constant_genes
        if constant_genes > n_genes * 0.5:
            warnings.append(f"常量基因过多 ({constant_genes})，可能存在技术问题")

    except Exception as e:
        errors.append(f"读取表达矩阵失败: {str(e)}")

    return errors + warnings, stats


def check_group_file(file_path: str, expr_samples: List[str] = None) -> Tuple[List[str], Dict]:
    """检查分组文件质量（纯Python实现）"""
    errors = []
    warnings = []
    stats = {}

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            rows = list(reader)

        if len(rows) < 2:
            errors.append("分组文件为空或只有表头")
            return errors, stats

        # 检查列名
        headers = [h.lower().strip() for h in rows[0]]
        sample_col_idx = None
        group_col_idx = None

        for i, h in enumerate(headers):
            if h == 'sample' or h == 'id':
                sample_col_idx = i
            elif h == 'group' or h == 'group_':
                group_col_idx = i

        if sample_col_idx is None or group_col_idx is None:
            # 检查常见错误
            current_cols = rows[0] if rows else []
            has_samples = any('sample' in str(c).lower() for c in current_cols)
            has_group = any('group' in str(c).lower() for c in current_cols)

            if 'Unnamed' in str(current_cols[0]) if current_cols else False:
                errors.append(f"分组文件格式错误：第一列是索引列（'Unnamed: 0'）")
                errors.append(f"  修复方法：将第一列删除或重命名为 'sample'")
            elif not has_samples and not has_group:
                errors.append(f"分组文件格式错误：必须包含 'sample' 和 'group' 列")
                errors.append(f"  当前列名: {current_cols}")
                errors.append(f"  修复方法：确保CSV包含 sample,group 两列")
            else:
                errors.append(f"分组文件必须包含 'sample' 和 'group' 列，当前列: {current_cols}")
                errors.append(f"  修复方法：将列名改为 sample,group")
            return errors, stats

        # 读取数据
        group_counts = {}
        samples = []

        for row in rows[1:]:
            if len(row) <= max(sample_col_idx, group_col_idx):
                continue
            sample_name = row[sample_col_idx].strip()
            group_name = row[group_col_idx].strip()
            samples.append(sample_name)
            group_counts[group_name] = group_counts.get(group_name, 0) + 1

        stats['n_samples'] = len(samples)
        stats['groups'] = group_counts

        if len(group_counts) < 2:
            warnings.append("分组数量少于2组，无法进行组间比较")

        # 检查样本匹配
        if expr_samples:
            expr_set = set(expr_samples)
            group_set = set(samples)

            missing_in_group = expr_set - group_set
            missing_in_expr = group_set - expr_set

            stats['matched_samples'] = len(expr_set & group_set)
            stats['missing_in_group'] = list(missing_in_group)[:10]
            stats['missing_in_expr'] = list(missing_in_expr)[:10]

            if missing_in_group:
                warnings.append(f"表达矩阵中有 {len(missing_in_group)} 个样本不在分组文件中")
            if missing_in_expr:
                warnings.append(f"分组文件中有 {len(missing_in_expr)} 个样本不在表达矩阵中")

            # 检查分组样本数
            for group, count in group_counts.items():
                matched_count = len([s for s in samples if s in expr_set and
                                   rows[1:][samples.index(s)][group_col_idx].strip() == group])
                if matched_count < 3:
                    warnings.append(f"分组 '{group}' 的有效样本数过少 ({matched_count})")

    except Exception as e:
        errors.append(f"读取分组文件失败: {str(e)}")

    return errors + warnings, stats


def run_quality_control(expression_file: str, group_file: str) -> Tuple[bool, str]:
    """
    运行质控检查（纯Python实现）
    """
    report = []
    report.append("=" * 60)
    report.append("免疫浸润分析 - 数据质控报告")
    report.append("=" * 60)
    report.append("")

    all_passed = True

    # 1. 表达矩阵检查
    report.append("[1] 表达矩阵检查")
    report.append("-" * 40)
    expr_errors, expr_stats = check_expression_matrix(expression_file)

    if expr_errors:
        for err in expr_errors:
            report.append(f"  [错误] {err}")
        all_passed = False
    else:
        report.append(f"  [通过] 基因数: {expr_stats.get('n_genes', 'N/A')}")
        report.append(f"  [通过] 样本数: {expr_stats.get('n_samples', 'N/A')}")

    if expr_stats.get('missing_pct', 0) > 0:
        report.append(f"  [信息] 缺失值: {expr_stats.get('missing_count', 0)} ({expr_stats.get('missing_pct', 0)}%)")

    if expr_stats.get('constant_genes', 0) > 0:
        report.append(f"  [警告] 常量基因: {expr_stats['constant_genes']}")

    report.append("")

    # 2. 分组文件检查
    report.append("[2] 分组文件检查")
    report.append("-" * 40)

    # 获取表达矩阵的样本列表
    try:
        with open(expression_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader)
            expr_samples = header[1:]  # 排除 SYMBOL 列
    except:
        expr_samples = None

    group_errors, group_stats = check_group_file(group_file, expr_samples)

    if group_errors:
        for err in group_errors:
            report.append(f"  [错误] {err}")
        all_passed = False
    else:
        report.append(f"  [通过] 样本数: {group_stats.get('n_samples', 'N/A')}")
        if 'groups' in group_stats:
            for group, count in group_stats['groups'].items():
                report.append(f"  [通过] {group}: {count} 样本")

    if 'matched_samples' in group_stats:
        report.append(f"  [信息] 匹配样本数: {group_stats['matched_samples']}")

    report.append("")

    # 3. 总结
    report.append("=" * 60)
    if all_passed:
        report.append("质控结果: 通过")
    else:
        report.append("质控结果: 失败")
        report.append("请修复上述错误后重新运行")
    report.append("=" * 60)

    report_content = "\n".join(report)

    if not all_passed:
        logging.error("\n" + report_content)
    else:
        logging.info("\n" + report_content)

    return all_passed, report_content


def build_r_command(method_info: Dict, config: configparser.ConfigParser, script_dir: str, force: bool = False) -> List[str]:
    """构建R命令"""
    script_path = os.path.join(script_dir, method_info['script'])
    
    cmd = ['Rscript', script_path]
    
    # 通用参数
    cmd.extend(['-e', config['Input']['expression']])
    cmd.extend(['-g', config['Input']['group']])
    
    # 可选: genes参数
    genes = config['Input']['genes']
    if genes and genes.strip() and method_info['has_genes']:
        cmd.extend(['--genes', genes])
    
    # 可选: signature参数
    signature = config['Input']['signature']
    if signature and signature.strip() and method_info['has_signature']:
        cmd.extend(['--signature', signature])
    
    # 输出目录
    base_dir = config['Output']['base_dir']
    if config.getboolean('Output', 'create_subdirs', fallback=True):
        output_dir = os.path.join(base_dir, f"{method_info['name']}_output")
    else:
        output_dir = base_dir
    cmd.extend(['-o', output_dir])
    
    # arrays参数(默认为false)
    cmd.extend(['--arrays', 'false'])
    
    # control参数(默认为Normal)
    cmd.extend(['--control', 'Normal'])
    
    # cache参数
    cache_dir = config['Cache']['cache_dir']
    cmd.extend(['--cache', cache_dir])
    
    # force参数 (R的logical参数需要TRUE/FALSE)
    if force:
        cmd.extend(['--force', 'TRUE'])

    # 方法特定参数
    if method_info['name'] == 'cibersort':
        perm = config.get('Method1.CIBERSORT', 'perm', fallback='1000')
        cmd.extend(['-p', perm])
    elif method_info['name'] == 'timer':
        tissue = config.get('Method7.TIMER', 'tissue', fallback='LUAD')
        cmd.extend(['--tissue', tissue])
    
    return cmd


def execute_task(task: Tuple[int, Dict, str], verbose: bool, cwd: str = None) -> Dict:
    """执行单个任务"""
    method_num, method_info, cmd = task

    result = {
        'method': method_info['name'],
        'num': method_num,
        'success': False,
        'output': '',
        'error': None
    }

    try:
        if verbose:
            logging.info(f"执行命令: {' '.join(cmd)}")

        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1小时超时
            cwd=cwd  # 设置工作目录
        )

        result['output'] = process.stdout
        if process.stderr:
            result['output'] += f"\nSTDERR:\n{process.stderr}"

        if process.returncode == 0:
            result['success'] = True
            logging.info(f"✓ {method_info['name']} 完成")
        else:
            result['error'] = f"返回码: {process.returncode}"
            # 增强错误日志：打印关键错误信息
            logging.error(f"✗ {method_info['name']} 失败: {result['error']}")

            # 提取错误信息中的关键部分
            if verbose:
                error_lines = [l for l in process.stderr.split('\n') if 'Error' in l or 'error' in l]
                if error_lines:
                    logging.error(f"  关键错误: {error_lines[:3]}")  # 只显示前3条

                # 如果返回码是1，可能是R脚本内部错误
                if process.returncode == 1:
                    # 尝试提取R错误
                    r_errors = [l for l in process.stderr.split('\n')
                               if any(x in l for x in ['Error', '错误', 'cannot', 'not found'])]
                    if r_errors:
                        logging.error(f"  R脚本错误: {r_errors[0]}")

    except subprocess.TimeoutExpired:
        result['error'] = '执行超时(>1小时)'
        logging.error(f"✗ {method_info['name']} 失败: {result['error']}")
    except Exception as e:
        result['error'] = str(e)
        logging.error(f"✗ {method_info['name']} 失败: {result['error']}", exc_info=True)

    return result


def execute_tasks(tasks: List[Tuple[int, Dict, str]], config: configparser.ConfigParser, verbose: bool, cwd: str = None):
    """并行执行任务"""
    results = []

    parallel_enabled = config.getboolean('Execution', 'parallel_enabled', fallback=True)
    max_workers = config.getint('Execution', 'max_workers', fallback=8)

    if not parallel_enabled or max_workers <= 1:
        # 串行执行
        logging.info("串行执行任务...")
        for task in tasks:
            result = execute_task(task, verbose, cwd)
            results.append(result)
    else:
        # 并行执行
        logging.info(f"并行执行任务 (max_workers={max_workers})...")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {
                executor.submit(execute_task, task, verbose, cwd): task
                for task in tasks
            }

            for future in as_completed(future_to_task):
                result = future.result()
                results.append(result)

    return results


def print_summary(results: List[Dict]):
    """打印执行总结"""
    print("\n" + "="*60)
    print("执行总结")
    print("="*60)

    success_count = sum(1 for r in results if r['success'])
    fail_count = len(results) - success_count

    print(f"\n成功: {success_count}/{len(results)}")
    print(f"失败: {fail_count}/{len(results)}")

    if fail_count > 0:
        print("\n失败的方法:")
        for r in results:
            if not r['success']:
                print(f"  - {r['method']}: {r['error']}")


def generate_report(results: List[Dict], config: configparser.ConfigParser,
                    output_dir: str, config_path: str, start_time: datetime) -> str:
    """生成执行报告"""
    end_time = datetime.now()
    duration = end_time - start_time

    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("免疫浸润分析报告")
    report_lines.append("=" * 70)
    report_lines.append("")

    # 基本信息
    report_lines.append("[执行信息]")
    report_lines.append(f"  开始时间: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"  结束时间: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"  总耗时: {duration}")
    report_lines.append(f"  配置文件: {config_path}")
    report_lines.append(f"  输出目录: {output_dir}")
    report_lines.append("")

    # 输入文件
    report_lines.append("[输入文件]")
    report_lines.append(f"  表达矩阵: {config['Input']['expression']}")
    report_lines.append(f"  分组信息: {config['Input']['group']}")
    genes = config['Input'].get('genes', '')
    if genes and genes.strip():
        report_lines.append(f"  基因列表: {genes}")
    signature = config['Input'].get('signature', '')
    if signature and signature.strip():
        report_lines.append(f"  基因签名: {signature}")
    report_lines.append("")

    # 执行统计
    success_count = sum(1 for r in results if r['success'])
    fail_count = len(results) - success_count
    report_lines.append("[执行统计]")
    report_lines.append(f"  总方法数: {len(results)}")
    report_lines.append(f"  成功: {success_count}")
    report_lines.append(f"  失败: {fail_count}")
    report_lines.append("")

    # 方法执行详情
    report_lines.append("[方法执行详情]")
    report_lines.append(f"  {'方法':<12} {'状态':<8} {'说明'}")
    report_lines.append("  " + "-" * 50)
    for r in sorted(results, key=lambda x: x['num']):
        status = "成功" if r['success'] else "失败"
        error_info = r.get('error', '') if not r['success'] else ''
        line = f"  {r['method']:<12} {status:<8}"
        if error_info:
            line += f" {error_info}"
        report_lines.append(line)
    report_lines.append("")

    # 输出文件列表
    if output_dir and os.path.exists(output_dir):
        report_lines.append("[输出文件]")
        for item in sorted(os.listdir(output_dir)):
            item_path = os.path.join(output_dir, item)
            if os.path.isdir(item_path) and not item.startswith('.'):
                # 列出子目录中的文件
                files = [f for f in os.listdir(item_path) if not f.startswith('.')]
                if files:
                    report_lines.append(f"  {item}/")
                    for f in sorted(files)[:10]:  # 最多显示10个文件
                        report_lines.append(f"    - {f}")
                    if len(files) > 10:
                        report_lines.append(f"    ... 共 {len(files)} 个文件")
        report_lines.append("")

    # 失败详情
    if fail_count > 0:
        report_lines.append("[失败详情]")
        for r in results:
            if not r['success']:
                report_lines.append(f"  {r['method']}:")
                report_lines.append(f"    错误: {r.get('error', '未知')}")
                if r.get('output'):
                    # 显示部分输出
                    output_lines = r['output'].strip().split('\n')[-5:]
                    report_lines.append("    输出(最后5行):")
                    for line in output_lines:
                        report_lines.append(f"      {line}")
        report_lines.append("")

    report_lines.append("=" * 70)
    report_lines.append("报告结束")
    report_lines.append("=" * 70)

    return "\n".join(report_lines)


def main():
    parser = argparse.ArgumentParser(
        description='免疫浸润分析流水线 - 批量执行多种免疫浸润分析方法',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例用法:
  # 使用配置文件运行所有方法
  python run_immune_infiltration.py --config config.ini -o /path/to/output

  # 仅运行特定方法
  python run_immune_infiltration.py --config config.ini -o /path/to/output --methods 1,7

  # 使用4核并行运行
  python run_immune_infiltration.py --config config.ini -o /path/to/output --workers 4

  # 强制重新计算
  python run_immune_infiltration.py --config config.ini -o /path/to/output --force

  # 预览将要执行的命令
  python run_immune_infiltration.py --config config.ini -o /path/to/output --dry-run

  # 删除输出目录（只有创建者可以删除）
  python run_immune_infiltration.py --config config.ini -o /path/to/output --rm

可用方法:
  编号  方法        说明
   1    cibersort   CIBERSORT免疫细胞反卷积(支持permutations参数)
   2    epic        EPIC免疫浸润分析
   3    estimate    ESTIMATE免疫评分
   4    ips         IPS免疫表型评分
   5    mcpcounter  MCPcounter免疫细胞计数
   6    ssgsea      ssGSEA基因集富集分析(需要signature文件)
   7    timer       TIMER免疫浸润(需要tissue参数)
   8    xcell       xCell细胞富集评分
        '''
    )
    
    parser.add_argument('--config', required=True, help='配置文件路径(INI格式)')
    parser.add_argument('-o', '--output', dest='output_dir', help='输出目录(覆盖配置文件中的base_dir)')
    parser.add_argument('--methods', help='要执行的方法编号列表,逗号分隔 (如: 1,7 或 1-5)')
    parser.add_argument('--workers', type=int, help='覆盖最大并行任务数')
    parser.add_argument('--force', action='store_true', help='强制重新计算,忽略缓存')
    parser.add_argument('--dry-run', action='store_true', help='仅显示将要执行的命令,不实际运行')
    parser.add_argument('--verbose', action='store_true', help='详细输出模式,显示更多调试信息')
    parser.add_argument('--rm', dest='remove_output', action='store_true', help='删除输出目录')
    parser.add_argument('-i', '--interactive', dest='interactive', action='store_true',
                        help='交互式选择要运行的方法')

    # Quarto 报告参数
    parser.add_argument('--quarto-report', action='store_true', help='生成 Quarto 报告 (需要安装 Quarto)')
    parser.add_argument('--quarto-format', default='docx', choices=['docx', 'html', 'pdf'],
                        help='Quarto 报告格式 [默认: docx]')
    parser.add_argument('--quarto-author', default='IKL', help='Quarto 报告作者 [默认: IKL]')

    args = parser.parse_args()

    # 删除模式
    if args.remove_output:
        if not args.output_dir:
            print("错误: 删除模式需要指定输出目录 -o <目录路径>", file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(args.output_dir):
            print(f"错误: 目录不存在: {args.output_dir}", file=sys.stderr)
            sys.exit(1)
        shutil.rmtree(args.output_dir)
        print(f"已删除: {args.output_dir}")
        sys.exit(0)
    
    # 生成时间戳
    timestamp = get_timestamp()

    # 设置日志 - 使用标准命名格式
    script_dir = Path(__file__).parent
    log_file = script_dir / ".." / "logs" / f"run_immune_infiltration.{timestamp}.log"
    os.makedirs(log_file.parent, exist_ok=True)
    setup_logging(str(log_file), args.verbose)

    # 记录开始时间
    start_time = datetime.now()

    logging.info("免疫浸润分析流水线启动")
    logging.info(f"配置文件: {args.config}")

    # 加载配置
    try:
        config = load_config(args.config)
    except Exception as e:
        logging.error(f"加载配置文件失败: {e}")
        sys.exit(1)

    # 解析方法列表
    try:
        methods = parse_methods(config['Methods']['methods'])
    except Exception as e:
        logging.error(f"解析方法列表失败: {e}")
        sys.exit(1)

    # 写入标准日志头部
    write_log_header(config, methods, timestamp)
    
    # 覆盖配置（如果命令行指定）
    if args.methods:
        config['Methods']['methods'] = args.methods
    if args.workers:
        config['Execution']['max_workers'] = str(args.workers)
    if args.output_dir:
        # 确保输出目录存在
        os.makedirs(args.output_dir, exist_ok=True)
        config['Output']['base_dir'] = args.output_dir
    
    # 全面验证
    logging.info("验证配置...")
    errors, warnings = validate_config(config)
    
    if errors:
        logging.error("配置验证失败:")
        for err in errors:
            logging.error(f"  - {err}")
        sys.exit(1)
    
    if warnings:
        logging.warning("警告:")
        for warn in warnings:
            logging.warning(f"  - {warn}")

    # 获取输入文件路径
    script_dir = Path(__file__).parent
    expression_path = config['Input']['expression']
    group_path = config['Input']['group']

    # 如果是相对路径，尝试相对于脚本目录解析
    if not os.path.isabs(expression_path):
        expression_path = script_dir / expression_path
    if not os.path.isabs(group_path):
        group_path = script_dir / group_path

    # 运行质控检查
    logging.info("运行数据质控...")
    qc_passed, qc_report = run_quality_control(str(expression_path), str(group_path))

    # 质控失败时可以选择退出或继续（这里选择退出）
    if not qc_passed:
        logging.error("数据质控未通过，请检查输入数据")
        # 将质控报告保存到文件
        if args.output_dir:
            qc_report_file = os.path.join(args.output_dir, "quality_control_report.txt")
            os.makedirs(args.output_dir, exist_ok=True)
            with open(qc_report_file, 'w', encoding='utf-8') as f:
                f.write(qc_report)
            logging.info(f"质控报告已保存: {qc_report_file}")
        sys.exit(1)

    # 交互式模式 - 选择输入文件
    if args.interactive:
        print("\n" + "=" * 60)
        print("免疫浸润分析 - 交互式模式")
        print("=" * 60)

        # 选择表达矩阵文件
        print("\n[1] 选择表达矩阵文件")
        print("  请输入表达矩阵文件路径（CSV格式，第一列为SYMBOL）")
        print("  直接回车使用默认: expr_tpm.csv")
        expr_input = input("  > ").strip()
        if not expr_input:
            expr_input = "expr_tpm.csv"
        config['Input']['expression'] = expr_input

        # 选择分组文件
        print("\n[2] 选择分组文件")
        print("  请输入分组文件路径（CSV格式，需包含sample和group列）")
        print("  直接回车使用默认: group.csv")
        group_input = input("  > ").strip()
        if not group_input:
            group_input = "group.csv"
        config['Input']['group'] = group_input

        # 选择输出目录
        print("\n[3] 选择输出目录")
        print("  请输入输出目录路径")
        print("  直接回车使用默认: immune_infiltration_results")
        output_input = input("  > ").strip()
        if output_input:
            config['Output']['base_dir'] = output_input
            os.makedirs(output_input, exist_ok=True)

        # 选择是否进行基因相关性分析
        print("\n[4] 基因相关性分析（可选）")
        print("  是否需要分析免疫细胞与基因的相关性？[y/N]")
        genes_input = input("  > ").strip().lower()
        if genes_input == 'y':
            print("  请输入基因列表文件路径（CSV格式，需包含gene列）")
            genes_file = input("  > ").strip()
            if genes_file:
                config['Input']['genes'] = genes_file

    # 显示genes参数提示
    show_genes_tips(config)

    # 解析方法列表
    try:
        methods = parse_methods(config['Methods']['methods'])
    except ValueError as e:
        logging.error(f"解析方法列表失败: {e}")
        sys.exit(1)

    # 交互式选择算法（如果用户未指定方法或使用 interactive 参数）
    if args.interactive:
        print("\n" + "=" * 60)
        print("请选择要运行的免疫浸润分析方法：")
        print("=" * 60)
        print("\n可用方法：")
        for num, info in METHOD_MAP.items():
            print(f"  {num}. {info['name']:12s} - {'支持基因相关性' if info['has_genes'] else ''}")

        print("\n输入示例：")
        print("  - 1,3,5   : 运行方法1、3、5")
        print("  - 1-4     : 运行方法1到4")
        print("  - all     : 运行所有方法")
        print("  - 直接回车 : 运行所有方法")

        while True:
            user_input = input("\n请输入方法编号 (1-8) 或 'all': ").strip().lower()

            if not user_input or user_input == 'all':
                methods = list(METHOD_MAP.keys())
                break
            elif user_input == 'q' or user_input == 'quit':
                print("已取消")
                sys.exit(0)
            else:
                try:
                    methods = parse_methods(user_input)
                    # 验证方法编号
                    invalid = [m for m in methods if m not in METHOD_MAP]
                    if invalid:
                        print(f"错误: 无效的方法编号 {invalid}，请重新输入")
                        continue
                    break
                except:
                    print("格式错误，请重新输入")

        # 更新配置
        config['Methods']['methods'] = ','.join(map(str, methods))
        print(f"\n已选择方法: {[METHOD_MAP[m]['name'] for m in methods if m in METHOD_MAP]}")

    logging.info(f"将要执行的方法: {', '.join([METHOD_MAP[m]['name'] for m in methods if m in METHOD_MAP])}")
    
    # 获取脚本目录
    script_dir = Path(__file__).parent

    # 验证R脚本
    logging.info("验证R脚本...")
    script_errors, script_warnings = validate_r_scripts(str(script_dir), methods)

    if script_errors:
        logging.error("R脚本验证失败:")
        for err in script_errors:
            logging.error(f"  - {err}")
        sys.exit(1)

    # 验证R脚本参数一致性（增强验证）
    logging.info("验证R脚本参数...")
    param_errors, param_warnings = validate_all_r_params(str(script_dir), methods)

    if param_errors:
        logging.error("R脚本参数验证失败:")
        for err in param_errors:
            logging.error(f"  - {err}")
        # 只记录警告，不中断执行（因为参数验证可能不完美）
        for warn in param_warnings:
            logging.warning(f"  - {warn}")

    # 构建任务 - 输出目录基于当前工作目录（不是配置文件目录）
    base_dir = config['Output']['base_dir']
    if not os.path.isabs(base_dir):
        # 基于当前工作目录解析
        base_dir = os.path.abspath(base_dir)
        config['Output']['base_dir'] = base_dir

    # qs2 缓存统一放在结果目录下 <base_dir>/qs2（覆盖配置文件中的 cache_dir）
    qs2_dir = os.path.join(base_dir, 'qs2')
    os.makedirs(qs2_dir, exist_ok=True)
    config['Cache']['cache_dir'] = qs2_dir
    logging.info(f"qs2 缓存目录（位于结果目录）: {qs2_dir}")

    tasks = []
    for m in methods:
        if m in METHOD_MAP:
            method_info = METHOD_MAP[m]
            cmd = build_r_command(method_info, config, str(script_dir), args.force)
            tasks.append((m, method_info, cmd))
    
    if not tasks:
        logging.error("没有可执行的任务")
        sys.exit(1)
    
    # Dry-run模式
    if args.dry_run:
        print("\n" + "="*60)
        print("Dry-run模式 - 将要执行的命令")
        print("="*60)
        for m, method_info, cmd in tasks:
            print(f"\n{method_info['name']}:")
            print("  " + " ".join(cmd))
        print("\n" + "="*60)
        sys.exit(0)
    
    # 执行任务
    logging.info("开始执行任务...")
    results = execute_tasks(tasks, config, args.verbose, cwd=str(script_dir))
    
    # 打印总结
    print_summary(results)

    # 生成报告
    if args.output_dir:
        report = generate_report(results, config, args.output_dir, args.config, start_time)
        report_file = os.path.join(args.output_dir, "report.txt")
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report)
        logging.info(f"报告已生成: {report_file}")
        print(f"\n报告已保存: {report_file}")

        # 生成 Quarto 报告
        if args.quarto_report:
            logging.info("生成 Quarto 报告...")
            quarto_path = shutil.which('quarto')
            if not quarto_path:
                logging.warning("Quarto 未安装，跳过报告生成。请访问 https://quarto.org 安装。")
            else:
                # 查找报告生成器
                script_dir = Path(__file__).parent
                report_generator = script_dir.parent.parent / "report_templates" / "generate_report.R"

                if report_generator.exists():
                    try:
                        quarto_cmd = [
                            'Rscript', str(report_generator),
                            '--results-dir', args.output_dir,
                            '--format', args.quarto_format,
                            '--project-name', os.path.basename(args.output_dir),
                            '--author', args.quarto_author
                        ]
                        subprocess.run(quarto_cmd, check=True)
                        logging.info("Quarto 报告已生成")
                    except subprocess.CalledProcessError as e:
                        logging.error(f"Quarto 报告生成失败: {e}")
                else:
                    logging.warning(f"报告生成器未找到: {report_generator}")

    # 记录结果统计
    success_count = sum(1 for r in results if r['success'])
    fail_count = len(results) - success_count
    output_dir = config['Output']['base_dir']

    # 获取输出文件列表
    output_files = []
    if os.path.exists(output_dir):
        for method_dir in os.listdir(output_dir):
            method_path = os.path.join(output_dir, method_dir)
            if os.path.isdir(method_path):
                for f in os.listdir(method_path):
                    if f.endswith(('.csv', '.png', '.pdf')):
                        output_files.append(f"{method_dir}/{f}")

    summary = f"""
================================================================================
执行总结
================================================================================
成功: {success_count}/{len(results)}
失败: {fail_count}/{len(results)}

----------------------------------------
输出文件
----------------------------------------
"""
    for f in output_files[:20]:  # 只显示前20个
        summary += f"  {f}\n"
    if len(output_files) > 20:
        summary += f"  ... 还有 {len(output_files) - 20} 个文件\n"

    summary += f"""
----------------------------------------
追溯信息
----------------------------------------
日志文件: run_immune_infiltration.{timestamp}.log
配置文件: {os.path.basename(args.config)}
================================================================================
"""
    if fail_count == 0:
        summary += "分析成功完成\n"
    else:
        summary += "分析部分完成，部分方法失败\n"
    summary += "================================================================================\n"

    logging.info(summary)

    # 将 R sessionInfo() 写入 logs/（与本次运行时间戳一致）
    session_info_path = Path(log_file).parent / f"sessionInfo.{timestamp}.txt"
    wsi_script = script_dir / "write_session_info.R"
    try:
        if wsi_script.exists():
            si_proc = subprocess.run(
                ["Rscript", str(wsi_script), "-o", str(session_info_path)],
                capture_output=True,
                text=True,
                timeout=120,
                cwd=str(script_dir),
            )
            if si_proc.returncode == 0:
                logging.info(f"R sessionInfo 已保存: {session_info_path}")
            else:
                logging.warning(
                    f"sessionInfo 写入返回码 {si_proc.returncode}: "
                    f"{(si_proc.stderr or '')[:500]}"
                )
        else:
            logging.warning(f"未找到 {wsi_script}，跳过 sessionInfo 导出")
    except Exception as e:
        logging.warning(f"sessionInfo 导出失败: {e}")

    # 检查是否有失败的任务
    if not all(r['success'] for r in results):
        sys.exit(1)


if __name__ == '__main__':
    main()
