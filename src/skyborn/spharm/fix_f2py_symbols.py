#!/usr/bin/env python3
"""
Batch script to fix F2PY symbol generation issues
Moves functions from modules to standalone subroutines to generate correct f2py symbols
"""

import os
import re
from pathlib import Path

# Files and function mappings to process
FILES_TO_PROCESS = {
    "shaec.f90": ["shaec", "shaeci"],
    "shagc.f90": ["shagc", "shagci"],
    "shags.f90": ["shags", "shagsi"],
    "shsec.f90": ["shsec", "shseci"],
    "shses.f90": ["shses", "shsesi"],
    "shsgc.f90": ["shsgc", "shsgci"],
    "shsgs.f90": ["shsgs", "shsgsi"],
    "vhaec.f90": ["vhaec", "vhaeci"],
    "vhaes.f90": ["vhaes", "vhaesi"],
    "vhagc.f90": ["vhagc", "vhagci"],
    "vhags.f90": ["vhags", "vhagsi"],
    "vhsec.f90": ["vhsec", "vhseci"],
    "vhses.f90": ["vhses", "vhsesi"],
    "vhsgc.f90": ["vhsgc", "vhsgci"],
    "vhsgs.f90": ["vhsgs", "vhsgsi"],
}


def extract_subroutine_from_module(content, subroutine_name):
    """Extract subroutine from module"""
    # Match subroutine start and end
    pattern = rf"(\s*!>.*?\n)?\s*subroutine\s+{subroutine_name}\s*\([^)]*\).*?end\s+subroutine\s+{subroutine_name}"
    match = re.search(pattern, content, re.DOTALL | re.IGNORECASE)

    if match:
        return match.group(0)
    return None


def process_file(file_path, functions_to_extract):
    """Process a single file"""
    print(f"Processing: {file_path}")

    with open(file_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Extract module dependencies from header
    module_deps = []
    if "use iso_fortran_env" in content:
        module_deps.append("   use iso_fortran_env, only: real64, int32")
    if "use sphcom_mod" in content:
        module_deps.append("   use sphcom_mod")
    if "use hrfft_mod" in content:
        module_deps.append("   use hrfft_mod")

    # Extract parameter definitions
    param_lines = []
    if "wp = real64" in content:
        param_lines.extend(
            [
                "   integer, parameter :: wp = real64",
                "   integer, parameter :: ip = int32",
            ]
        )

    # Find and extract each function
    extracted_functions = []

    for func_name in functions_to_extract:
        func_content = extract_subroutine_from_module(content, func_name)
        if func_content:
            # Create independent subroutine
            independent_func = f"""
! Independent subroutine for f2py compatibility - must be outside module to generate {func_name}_ symbol
subroutine {func_name}(...) ! Copy complete signature from original function here
{chr(10).join(module_deps)}
   implicit none

{chr(10).join(param_lines)}

   ! Copy parameter declarations and implementation from original function here
   ! {func_content}

end subroutine {func_name}"""
            extracted_functions.append(independent_func)
        else:
            print(f"Warning: Function {func_name} not found in file {file_path}")

    # Modify original file: clear module, add comment, then add independent functions
    module_name = file_path.stem + "_mod"

    new_content = f"""!> @file {file_path.name}
!> @brief SPHEREPACK function - OPTIMIZED for modern Fortran
!> @author SPHEREPACK team, optimized by Qianye Su
!> @date 2025

module {module_name}
{chr(10).join(module_deps)}
   implicit none
   private

   ! Module parameters
{chr(10).join(param_lines)}

   ! Note: {', '.join(functions_to_extract)} are now independent subroutines below for f2py compatibility

end module {module_name}

{chr(10).join(extracted_functions)}
"""

    # Backup original file
    backup_path = file_path.with_suffix(file_path.suffix + ".backup")
    if not backup_path.exists():
        with open(backup_path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"  Backup created: {backup_path}")

    # Write new content
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(new_content)

    print(f"  Processed: {len(functions_to_extract)} functions")


def main():
    """Main function"""
    src_dir = Path("/mnt/d/skyborn/src/skyborn/spharm/src")

    if not src_dir.exists():
        print(f"Error: Source directory does not exist: {src_dir}")
        return

    print("Starting batch processing of F2PY symbol issues...")
    print(f"Total files to process: {len(FILES_TO_PROCESS)}")

    for filename, functions in FILES_TO_PROCESS.items():
        file_path = src_dir / filename

        if file_path.exists():
            process_file(file_path, functions)
        else:
            print(f"Warning: File does not exist: {file_path}")

    print("\n✅ Batch processing complete!")
    print(
        "⚠️  Note: This is a template script, generated functions need manual completion"
    )
    print("Recommendation: Manually modify files one by one to ensure code correctness")


if __name__ == "__main__":
    main()
