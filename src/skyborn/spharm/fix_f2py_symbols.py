#!/usr/bin/env python3
"""
Script to generate F2PY interface files for Fortran 90 files
Creates a .pyf signature file for all .f90 files in the source directory

Usage:
    # Run from the spharm directory (parent of src/)
    python fix_f2py_symbols.py

    # Or run from any directory
    python /path/to/fix_f2py_symbols.py
"""

import os
import subprocess
import sys
from pathlib import Path


def find_all_f90_files(src_dir):
    """Find all .f90 files in the source directory"""
    f90_files = list(Path(src_dir).glob("*.f"))
    return sorted(f90_files)


def generate_pyf_signature(src_dir, output_file="_spherepack.pyf"):
    """Generate .pyf signature file for all .f90 files"""
    f90_files = find_all_f90_files(src_dir)

    if not f90_files:
        print(f"No .f90 files found in {src_dir}")
        return False

    print(f"Found {len(f90_files)} .f90 files:")
    for f in f90_files:
        print(f"  - {f.name}")

    # Prepare f2py command
    output_path = src_dir / output_file
    cmd = [
        sys.executable,
        "-m",
        "numpy.f2py",
        "-h",
        str(output_path),
        "-m",
        "_spherepack",
    ]

    # Add all .f90 files to the command
    for f90_file in f90_files:
        cmd.append(str(f90_file))

    print(f"\nGenerating signature file: {output_path}")
    print(f"Command: {' '.join(cmd[:7])} ... [files]")

    try:
        # Execute f2py command
        result = subprocess.run(
            cmd, cwd=src_dir, capture_output=True, text=True, timeout=300
        )

        if result.returncode == 0:
            print(f"✅ Successfully generated: {output_path}")
            if result.stdout:
                print("F2PY output:")
                print(result.stdout)
            return True
        else:
            print(f"❌ Error generating signature file!")
            print(f"Return code: {result.returncode}")
            if result.stderr:
                print("Error output:")
                print(result.stderr)
            return False

    except subprocess.TimeoutExpired:
        print("❌ F2PY command timed out (>300s)")
        return False
    except Exception as e:
        print(f"❌ Error running f2py: {e}")
        return False


def verify_pyf_file(pyf_path):
    """Verify the generated .pyf file"""
    if not pyf_path.exists():
        return False

    with open(pyf_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Basic checks
    if "python module" in content and "interface" in content:
        lines = content.count("\n")
        print(f"✅ Generated .pyf file looks valid ({lines} lines)")

        # Count subroutines
        subroutine_count = content.count("subroutine ")
        function_count = content.count("function ")
        print(
            f"   - Found {subroutine_count} subroutines and {function_count} functions"
        )
        return True
    else:
        print("❌ Generated .pyf file appears invalid")
        return False


def main():
    """Main function"""
    # Use relative path - script should be run from the parent directory of src/
    script_dir = Path(__file__).parent
    src_dir = script_dir / "src"

    # Alternative: automatically detect src directory
    if not src_dir.exists():
        # Try current directory + src
        src_dir = Path.cwd() / "src"

    if not src_dir.exists():
        print(f"❌ Error: Source directory not found!")
        print(f"   Tried: {script_dir / 'src'}")
        print(f"   Tried: {Path.cwd() / 'src'}")
        print(
            f"💡 Please run this script from the spharm directory containing 'src/' folder"
        )
        return

    print("🚀 Starting F2PY signature file generation for all .f90 files...")
    print(f"📁 Source directory: {src_dir}")

    # Generate the .pyf signature file
    output_file = "_spherepack.pyf"
    success = generate_pyf_signature(src_dir, output_file)

    if success:
        # Verify the generated file
        pyf_path = src_dir / output_file
        if verify_pyf_file(pyf_path):
            print(f"\n🎉 Success! Generated F2PY signature file: {output_file}")
            print(f"📄 Location: {pyf_path}")
            print("\n📋 Next steps:")
            print("1. Review and edit the .pyf file if needed")
            print("2. Use the .pyf file to compile the Fortran 90 extension:")
            print(f"   python -m numpy.f2py -c {output_file} src/*.f90")
            print("3. Import the compiled module in Python:")
            print("   import _spherepack_f90")
        else:
            print(f"\n❌ Generated file {output_file} appears to be invalid")
    else:
        print(f"\n❌ Failed to generate F2PY signature file")
        print("\n🔧 Troubleshooting tips:")
        print("1. Check that all .f90 files have valid Fortran syntax")
        print("2. Ensure numpy and f2py are properly installed")
        print("3. Try processing individual files to identify issues")


if __name__ == "__main__":
    main()
