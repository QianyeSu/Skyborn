# User-friendly installation script
"""
Intelligent Installation Helper - Provides optimal installation experience for different platforms and scenarios
"""
import platform
import shutil
import subprocess
import sys
from pathlib import Path


class SkyornInstallHelper:
    """Intelligent Installation Helper"""

    PLATFORM_GUIDES = {
        "Linux": {
            "ubuntu/debian": "sudo apt-get install gfortran",
            "centos/rhel": "sudo yum install gcc-gfortran",
            "fedora": "sudo dnf install gcc-gfortran",
            "arch": "sudo pacman -S gcc-fortran",
            "conda": "conda install gfortran_linux-64",
        },
        "Darwin": {  # macOS
            "homebrew": "brew install gcc",
            "macports": "sudo port install gcc12 +gfortran",
            "conda": "conda install gfortran_osx-64",
        },
        "Windows": {
            "conda": "conda install m2w64-toolchain",
            "msys2": "pacman -S mingw-w64-x86_64-gcc-fortran",
            "intel": "Intel oneAPI HPC Toolkit",
        },
    }

    @classmethod
    def detect_platform_details(cls):
        """Detect detailed platform information"""
        system = platform.system()
        machine = platform.machine()

        # macOS special handling
        if system == "Darwin":
            # Detect Apple Silicon
            try:
                import subprocess

                result = subprocess.run(
                    ["sysctl", "-n", "hw.optional.arm64"],
                    capture_output=True,
                    text=True,
                )
                is_apple_silicon = result.stdout.strip() == "1"
            except:
                is_apple_silicon = machine == "arm64"

            return {
                "system": "macOS",
                "arch": "arm64" if is_apple_silicon else "x86_64",
                "is_apple_silicon": is_apple_silicon,
                "machine": machine,
            }

        return {"system": system, "arch": machine, "machine": machine}

    @classmethod
    def get_installation_guide(cls):
        """Get installation guide for current platform"""
        platform_info = cls.detect_platform_details()
        system = platform_info["system"]

        if system not in cls.PLATFORM_GUIDES:
            return cls._generic_guide()

        guides = cls.PLATFORM_GUIDES[system]

        # Generate detailed guide for each platform
        if system == "macOS":
            return cls._macos_guide(platform_info, guides)
        elif system == "Linux":
            return cls._linux_guide(guides)
        elif system == "Windows":
            return cls._windows_guide(guides)
        else:
            return cls._generic_guide()

    @classmethod
    def _macos_guide(cls, platform_info, guides):
        """macOS specific guide"""
        arch = platform_info["arch"]
        is_apple_silicon = platform_info.get("is_apple_silicon", False)

        guide = f"""
🍎 macOS {arch} Installation Guide:

{'⚠️  Apple Silicon Mac Users Note:' if is_apple_silicon else ''}
{'  Ensure compiler supports ARM64 architecture' if is_apple_silicon else ''}

Recommended Methods (by priority):

1️⃣ Homebrew (Recommended):
   {guides['homebrew']}

2️⃣ Conda (if you use Anaconda):
   {guides['conda']}

3️⃣ MacPorts:
   {guides['macports']}

Post-installation verification:
   gfortran --version

Then reinstall:
   pip install skyborn
"""
        return guide

    @classmethod
    def _linux_guide(cls, guides):
        """Linux distribution guide"""
        return f"""
🐧 Linux Installation Guide:

Choose based on your distribution:

Ubuntu/Debian: {guides['ubuntu/debian']}
CentOS/RHEL:   {guides['centos/rhel']}
Fedora:        {guides['fedora']}
Arch Linux:    {guides['arch']}
Conda users:     {guides['conda']}

Retry after installation:
   pip install skyborn
"""

    @classmethod
    def _windows_guide(cls, guides):
        """Windows guide"""
        return f"""
🪟 Windows Installation Guide:

Recommended Methods (by priority):

1️⃣ Anaconda/Miniconda (Easiest):
   {guides['conda']}

2️⃣ MSYS2 (Advanced users):
   {guides['msys2']}

3️⃣ Intel oneAPI (Professional users):
   {guides['intel']}

Note: Windows users are recommended to use Anaconda environment
"""

    @classmethod
    def _generic_guide(cls):
        """Generic guide"""
        return """
❓ Generic Installation Guide:

Your platform may require manual Fortran compiler installation.
Please visit project documentation for detailed instructions:
https://github.com/QianyeSu/Skyborn#installation

Or try source installation:
   git clone https://github.com/QianyeSu/Skyborn
   cd Skyborn
   pip install -e .
"""


def show_install_help():
    """Show installation help"""
    helper = SkyornInstallHelper()
    print(helper.get_installation_guide())


if __name__ == "__main__":
    show_install_help()
