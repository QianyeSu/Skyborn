#!/usr/bin/env python3
"""
Script to create a clean _spherepack.pyf file based on the definitions in _spherepack_old.pyf
"""

import os


def create_clean_pyf_file():
    """Create a clean pyf file based on the old file"""

    # File paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    old_file = os.path.join(script_dir, "_spherepack_old.pyf")
    output_file = os.path.join(script_dir, "_spherepack_clean.pyf")

    # Check if the old file exists
    if not os.path.exists(old_file):
        print(f"Error: File not found {old_file}")
        return

    print("=== Creating clean PYF file ===")
    print(f"Reference file: {old_file}")
    print(f"Output file: {output_file}")
    print()

    # Read the old file
    with open(old_file, "r", encoding="utf-8") as f:
        content = f.read()

    # Write to the new file
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(content)

    print("✓ Successfully created clean PYF file")
    print(f"Please rename {output_file} to _spherepack.pyf and rebuild")
    print()

    # Check if main functions exist
    main_functions = [
        "shaes",
        "shaesi",
        "shaec",
        "shaeci",
        "shags",
        "shagsi",
        "shagc",
        "shagci",
        "shses",
        "shsesi",
        "shsec",
        "shseci",
        "shsgs",
        "shsgsi",
        "shsgc",
        "shsgci",
        "vhaes",
        "vhaesi",
        "vhaec",
        "vhaeci",
        "vhags",
        "vhagsi",
        "vhagc",
        "vhagci",
        "vhses",
        "vhsesi",
        "vhsec",
        "vhseci",
        "vhsgs",
        "vhsgsi",
        "vhsgc",
        "vhsgci",
    ]

    print("Main function check:")
    for func in main_functions:
        if f"subroutine {func}" in content:
            print(f"✓ {func}")
        else:
            print(f"✗ {func} - Not found")


if __name__ == "__main__":
    create_clean_pyf_file()
