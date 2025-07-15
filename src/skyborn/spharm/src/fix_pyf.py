#!/usr/bin/env python3
"""
脚本用于修复 _spherepack.pyf 文件，使其与 _spherepack_old.pyf 中相同的 subroutine 保持一致的参数定义。
"""

import re
import os
from typing import Dict, List, Tuple


def parse_pyf_file(file_path: str) -> Dict[str, str]:
    """
    解析 .pyf 文件，返回 subroutine 名称到完整定义的映射
    """
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    # 使用正则表达式匹配 subroutine 定义
    # 匹配从 "subroutine name" 到 "end subroutine name" 的完整定义
    subroutine_pattern = r"subroutine\s+(\w+)\s*\([^)]*\)[^!]*![^!]*!\s*in\s*:[^:]*:([^:]+)\.f.*?end\s+subroutine\s+\1"

    subroutines = {}
    matches = re.finditer(subroutine_pattern, content, re.DOTALL | re.IGNORECASE)

    for match in matches:
        name = match.group(1).lower()
        full_definition = match.group(0)
        subroutines[name] = full_definition

    return subroutines


def extract_header_and_footer(file_path: str) -> Tuple[str, str]:
    """
    提取 .pyf 文件的头部和尾部（非 subroutine 部分）
    """
    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    # 找到第一个 subroutine 的位置
    first_subroutine_match = re.search(
        r"subroutine\s+\w+\s*\([^)]*\)[^!]*![^!]*!\s*in\s*:", content, re.IGNORECASE
    )

    if first_subroutine_match:
        header = content[: first_subroutine_match.start()]

        # 找到最后一个 end subroutine 的位置
        last_end_match = None
        for match in re.finditer(r"end\s+subroutine\s+\w+", content, re.IGNORECASE):
            last_end_match = match

        if last_end_match:
            footer = content[last_end_match.end() :]
        else:
            footer = "\n    end interface \nend python module _spherepack\n"
    else:
        header = content
        footer = ""

    return header, footer


def fix_pyf_file(old_file: str, new_file: str, output_file: str):
    """
    修复 .pyf 文件，使用 old_file 中的定义替换 new_file 中相同的 subroutine
    """
    print(f"正在解析 {old_file}...")
    old_subroutines = parse_pyf_file(old_file)

    print(f"正在解析 {new_file}...")
    new_subroutines = parse_pyf_file(new_file)

    print(f"找到 {len(old_subroutines)} 个老函数，{len(new_subroutines)} 个新函数")

    # 获取文件头部和尾部
    header, footer = extract_header_and_footer(new_file)

    # 创建输出内容
    output_content = header

    # 统计信息
    replaced_count = 0
    new_count = 0

    # 处理所有函数
    processed_functions = set()

    # 首先处理老文件中存在的函数（优先使用老定义）
    for func_name, old_definition in old_subroutines.items():
        if func_name in new_subroutines:
            print(f"✓ 替换函数: {func_name}")
            output_content += old_definition + "\n"
            replaced_count += 1
            processed_functions.add(func_name)
        else:
            print(f"+ 添加老函数: {func_name}")
            output_content += old_definition + "\n"
            new_count += 1
            processed_functions.add(func_name)

    # 然后处理新文件中独有的函数
    for func_name, new_definition in new_subroutines.items():
        if func_name not in processed_functions:
            print(f"○ 保留新函数: {func_name}")
            output_content += new_definition + "\n"
            new_count += 1

    output_content += footer

    # 写入输出文件
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(output_content)

    print(f"\n修复完成！")
    print(f"- 替换了 {replaced_count} 个函数")
    print(f"- 保留了 {new_count} 个新函数")
    print(f"- 输出文件: {output_file}")


def main():
    """主函数"""
    # 文件路径
    script_dir = os.path.dirname(os.path.abspath(__file__))
    old_file = os.path.join(script_dir, "_spherepack_old.pyf")
    new_file = os.path.join(script_dir, "_spherepack.pyf")
    output_file = os.path.join(script_dir, "_spherepack_fixed.pyf")

    # 检查文件是否存在
    if not os.path.exists(old_file):
        print(f"错误: 找不到文件 {old_file}")
        return

    if not os.path.exists(new_file):
        print(f"错误: 找不到文件 {new_file}")
        return

    print("=== PYF 文件修复工具 ===")
    print(f"老文件: {old_file}")
    print(f"新文件: {new_file}")
    print(f"输出文件: {output_file}")
    print()

    # 执行修复
    fix_pyf_file(old_file, new_file, output_file)

    print(f"\n修复完成！请将 {output_file} 重命名为 _spherepack.pyf 并重新构建。")


if __name__ == "__main__":
    main()
