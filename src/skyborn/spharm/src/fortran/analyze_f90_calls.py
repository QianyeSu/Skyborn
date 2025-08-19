#!/usr/bin/env python3
"""
Analyze F90 source files to count subroutine call frequencies
"""
import os
import re
from collections import defaultdict
import glob


def analyze_f90_files(directory):
    """Analyze F90 files and count subroutine calls"""
    call_counts = defaultdict(int)

    # Find all .f90 files recursively
    pattern = os.path.join(directory, "**", "*.f90")
    f90_files = glob.glob(pattern, recursive=True)

    print(f"Found {len(f90_files)} F90 files")

    # Regex pattern to match CALL statements
    call_pattern = re.compile(r"^\s*call\s+([a-zA-Z_][a-zA-Z0-9_]*)", re.IGNORECASE)

    for filepath in f90_files:
        try:
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                for line_num, line in enumerate(f, 1):
                    # Skip comment lines
                    if line.strip().startswith("!") or line.strip().startswith("C"):
                        continue

                    # Look for CALL statements
                    match = call_pattern.search(line)
                    if match:
                        subroutine_name = match.group(1)
                        call_counts[subroutine_name] += 1

        except Exception as e:
            print(f"Error reading {filepath}: {e}")

    return call_counts


def print_statistics(call_counts):
    """Print call frequency statistics"""
    print("\n=== Subroutine Call Frequency Analysis ===\n")

    total_calls = sum(call_counts.values())
    print(f"Total subroutine calls found: {total_calls}")
    print(f"Unique subroutines called: {len(call_counts)}")

    print("\n=== Top 20 Most Called Subroutines ===")
    print("Rank | Calls | Subroutine Name")
    print("-" * 40)

    sorted_calls = sorted(call_counts.items(), key=lambda x: x[1], reverse=True)

    for i, (subroutine, count) in enumerate(sorted_calls[:20], 1):
        percentage = (count / total_calls) * 100
        print(f"{i:4d} | {count:5d} | {subroutine} ({percentage:.1f}%)")

    print("\n=== Optimization Candidates ===")
    print("(Subroutines called 10+ times)")
    print("-" * 40)

    optimization_candidates = [
        (name, count) for name, count in sorted_calls if count >= 10
    ]
    for name, count in optimization_candidates:
        percentage = (count / total_calls) * 100
        print(f"{name}: {count} calls ({percentage:.1f}%)")


if __name__ == "__main__":
    # Directory containing F90 source files
    src_directory = r"../src"

    # Analyze the files
    call_counts = analyze_f90_files(src_directory)

    # Print statistics
    print_statistics(call_counts)

    # Save detailed results to file
    output_file = "f90_call_analysis.txt"
    with open(output_file, "w") as f:
        f.write("=== Complete Subroutine Call Analysis ===\n\n")
        f.write(f"Total calls: {sum(call_counts.values())}\n")
        f.write(f"Unique subroutines: {len(call_counts)}\n\n")

        sorted_calls = sorted(call_counts.items(), key=lambda x: x[1], reverse=True)
        for subroutine, count in sorted_calls:
            f.write(f"{subroutine}: {count}\n")

    print(f"\nDetailed results saved to: {output_file}")
