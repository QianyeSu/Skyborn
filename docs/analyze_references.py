#!/usr/bin/env python3
"""
分析functions_classes.rst中的引用，生成需要在API文档中添加:noindex:的列表
"""
import re
from pathlib import Path

def analyze_functions_classes_references():
    """分析functions_classes.rst中的所有:func:、:class:、:meth:引用"""

    functions_classes_path = Path("source/functions_classes.rst")

    if not functions_classes_path.exists():
        print(f"Error: {functions_classes_path} not found")
        return

    with open(functions_classes_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 正则表达式匹配所有的引用
    patterns = {
        'func': re.compile(r':func:`([^`]+)`'),
        'class': re.compile(r':class:`([^`]+)`'),
        'meth': re.compile(r':meth:`([^`]+)`')
    }

    references = {}

    for ref_type, pattern in patterns.items():
        matches = pattern.findall(content)
        references[ref_type] = []

        for match in matches:
            # 处理简化引用（以~开头）
            if match.startswith('~'):
                full_path = match[1:]  # 去掉~
            else:
                full_path = match

            references[ref_type].append(full_path)

    # 按模块分组
    modules = {}

    for ref_type, refs in references.items():
        for ref in refs:
            # 提取模块名
            parts = ref.split('.')
            if len(parts) >= 2:
                module = '.'.join(parts[:2])  # skyborn.calc, skyborn.gridfill等

                if module not in modules:
                    modules[module] = {'funcs': set(), 'classes': set(), 'methods': set()}

                if ref_type == 'func':
                    modules[module]['funcs'].add(ref)
                elif ref_type == 'class':
                    modules[module]['classes'].add(ref)
                elif ref_type == 'meth':
                    modules[module]['methods'].add(ref)

    # 输出分析结果
    print("=" * 60)
    print("FUNCTIONS_CLASSES.RST 引用分析")
    print("=" * 60)

    for module, refs in sorted(modules.items()):
        if any(refs.values()):  # 如果模块有任何引用
            print(f"\n{module}:")

            if refs['funcs']:
                print("  Functions:")
                for func in sorted(refs['funcs']):
                    print(f"    - {func}")

            if refs['classes']:
                print("  Classes:")
                for cls in sorted(refs['classes']):
                    print(f"    - {cls}")

            if refs['methods']:
                print("  Methods:")
                for method in sorted(refs['methods']):
                    print(f"    - {method}")

    # 生成需要修改的API文件列表
    print("\n" + "=" * 60)
    print("需要添加 :noindex: 的API文件:")
    print("=" * 60)

    api_file_mapping = {
        'skyborn.calc': 'api/calculations.rst',
        'skyborn.conversion': 'api/conversion.rst',
        'skyborn.gridfill': 'api/gridfill.rst',
        'skyborn.interp': 'api/interpolation.rst',
        'skyborn.gradients': 'api/gradients.rst',
        'skyborn.plot': 'api/plotting.rst',
        'skyborn.causality': 'api/causality.rst',
        'skyborn.spharm': 'api/spharm.rst',
        'skyborn.windspharm': 'api/windspharm.rst'
    }

    for module, api_file in api_file_mapping.items():
        if module in modules and any(modules[module].values()):
            print(f"\n{api_file}:")
            refs = modules[module]

            if refs['funcs']:
                print("  需要为以下函数添加 :noindex::")
                for func in sorted(refs['funcs']):
                    print(f"    - {func}")

            if refs['classes']:
                print("  需要为以下类添加 :noindex::")
                for cls in sorted(refs['classes']):
                    print(f"    - {cls}")

    return modules

if __name__ == "__main__":
    analyze_functions_classes_references()
