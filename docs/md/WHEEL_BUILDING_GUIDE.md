# Skyborn Wheel Building and Distribution Guide

## Overview

We have established a **comprehensive wheel building and distribution pipeline** for the Skyborn package, enabling seamless installation across multiple platforms and Python environments. This enterprise-grade distribution strategy ensures maximum accessibility for the atmospheric science research community.

## Distribution Architecture

Our multi-tier release strategy provides robust, scalable distribution:

### 1. Local Development Environment
- **Purpose**: Development, testing, and validation
- **Features**: Rapid iteration, debugging capabilities, custom optimization
- **Tools**: Local build scripts with comprehensive testing

### 2. Automated CI/CD Pipeline
- **Purpose**: Multi-platform wheel generation
- **Features**: Cross-platform compatibility, performance optimization, quality assurance
- **Infrastructure**: GitHub Actions with advanced cibuildwheel integration

### 3. PyPI Distribution Network
- **Purpose**: Global package distribution
- **Features**: Automated uploads, trusted publishing, version management
- **Access**: Universal pip installation for end users

## System Requirements and Platform Support

### Development Prerequisites

**Essential Build Dependencies:**
- **Python Environment**: 3.9 or higher with pip package manager
- **Fortran Compiler**: gfortran or Intel Fortran (platform-specific)
- **Build Tools**: meson, ninja, wheel, build (automatically managed)
- **Scientific Computing**: NumPy 1.21.0+ for array operations

### Supported Platform Matrix

| Platform | Architecture | Status | Python Versions |
|----------|-------------|--------|----------------|
| **Linux** | x86_64 | ✅ Production Ready | 3.9, 3.10, 3.11, 3.12 |
| **macOS** | Intel x86_64 | ✅ Production Ready | 3.9, 3.10, 3.11, 3.12 |
| **macOS** | Apple Silicon | ✅ Production Ready | 3.9, 3.10, 3.11, 3.12 |
| **Windows** | x86_64 | ✅ Production Ready | 3.9, 3.10, 3.11, 3.12 |

## Local Development and Testing

### Streamlined Build Process

We provide an automated build script for efficient local development:

```bash
# Complete build and test cycle
python build_wheel.py
```

**Automated Process Includes:**
- Dependency verification and installation
- Clean build environment preparation
- Optimized wheel compilation
- Comprehensive functionality testing
- Virtual environment validation

### Manual Build Process

For advanced users requiring granular control:

```bash
# Environment preparation
python -m pip install build wheel meson ninja numpy

# Clean previous builds
rm -rf build/ dist/ *.egg-info/

# Generate optimized wheel
python -m build --wheel --no-isolation

# Installation validation
pip install dist/*.whl
python -c "from skyborn.spharm import Spharmt; print('✅ Installation successful')"
```

## Advanced CI/CD Pipeline

### Automated Build Infrastructure

Our sophisticated GitHub Actions workflow leverages **cibuildwheel** for comprehensive multi-platform support:

**Trigger Conditions:**
- **Release Tags**: Automatic building on version tags (e.g., `v0.3.8`)
- **Platform Coverage**: Linux, macOS (Intel/ARM), Windows
- **Python Matrix**: Complete testing across Python 3.9-3.12
- **Architecture**: Optimized x86_64 builds (64-bit only)

### Pipeline Stages

#### 1. Build Environment Configuration
- **Compiler Setup**: Platform-specific Fortran compiler installation
- **Python Environments**: Isolated virtual environments for each Python version
- **Dependency Management**: Automated scientific computing stack preparation

#### 2. Performance-Optimized Compilation
- **Compiler Optimization**: Advanced flags for maximum performance
- **Multi-Version Building**: Parallel compilation for all supported Python versions
- **Quality Assurance**: Automated testing of each generated wheel

#### 3. Distribution and Release
- **Artifact Management**: Consolidated wheel collection from all platforms
- **GitHub Integration**: Automatic artifact storage and versioning
- **PyPI Publishing**: Trusted publishing with secure, token-less uploads

### Configuration Architecture

**Core Configuration Files:**
- **`.github/workflows/build-wheels.yml`**: Complete CI/CD orchestration
- **`pyproject.toml`**: Build system specifications and metadata
- **`meson.build`**: Fortran compilation optimization settings

## Professional Release Management

### Version Coordination Process

```bash
# Synchronized version update
python set_version.py 0.3.8

# Verification of version consistency
grep -n "version.*=" pyproject.toml src/skyborn/_version.py
```

### Comprehensive Testing Protocol

**Pre-Release Validation:**
```bash
# 1. Local build verification
python build_wheel.py

# 2. Isolated environment testing
python -m venv test_env
source test_env/bin/activate  # Linux/macOS
# test_env\Scripts\activate   # Windows
pip install dist/*.whl
python -c "from skyborn.spharm import Spharmt; sht = Spharmt(8,6); print('✅ Validation successful')"
deactivate
```

### Release Execution

**Professional Release Workflow:**
```bash
# 1. Commit preparation
git add .
git commit -m "Release v0.3.8: Enhanced wheel building with cross-platform optimization"

# 2. Version tagging and distribution
git tag v0.3.8
git push origin main --tags

# 3. Automated monitoring
# GitHub Actions: https://github.com/QianyeSu/Skyborn/actions
```

## Distribution Strategy and User Experience

### Primary Distribution: Pre-Compiled Wheels

**User Experience Benefits:**
- **Installation Speed**: `pip install skyborn` (sub-minute installation)
- **Zero Dependencies**: No compiler requirements for end users
- **Universal Compatibility**: Supports all major platforms and Python versions
- **Performance**: Pre-optimized with platform-specific compiler flags

### Fallback: Source Distribution

**Advanced User Options:**
- **Custom Compilation**: `pip install skyborn --no-binary=skyborn`
- **Platform Flexibility**: Support for specialized architectures
- **Development Mode**: `python -m pip install -e .` for contributors

### Installation Examples

```bash
# Recommended: Optimized wheel installation
pip install skyborn

# Advanced: Source compilation with custom flags
pip install skyborn --no-binary=skyborn

# Development: Editable installation
git clone https://github.com/QianyeSu/Skyborn.git
cd Skyborn
python -m pip install -e .
```

## Performance Engineering

### Compiler Optimization Strategy

**Platform-Specific Optimizations:**
- **Linux/macOS**: `-O3 -march=x86-64 -mtune=generic` for maximum compatibility
- **Windows**: `/O2` with MSVC optimization
- **Advanced Features**: OpenMP parallelization, automatic vectorization

### Distribution Optimization

**Wheel Size and Performance:**
- **Debug Symbol Removal**: Minimized distribution size
- **Generic CPU Targeting**: Broad compatibility without sacrificing performance
- **Selective Packaging**: Exclusion of development and testing artifacts

## Troubleshooting and Support

### Common Resolution Strategies

#### Fortran Compiler Issues
```bash
# Ubuntu/Debian
sudo apt-get update && sudo apt-get install gfortran

# macOS
brew install gcc

# Windows (Anaconda)
conda install -c conda-forge fortran-compiler
```

#### Import and Compatibility Issues
- **Python Version**: Verify Python 3.9+ compatibility
- **Architecture Matching**: Ensure wheel architecture compatibility
- **Fallback Strategy**: Source installation for unsupported configurations

#### Build Environment Issues
```bash
# Tool chain updates
pip install --upgrade meson ninja build wheel

# Dependency verification
pip install "numpy>=1.21.0"

# Clean build environment
rm -rf build/ dist/ *.egg-info/
```

## Community and Professional Support

### Support Channels
- **Issue Tracking**: [GitHub Issues](https://github.com/QianyeSu/Skyborn/issues) for bug reports and feature requests
- **Community Discussion**: [GitHub Discussions](https://github.com/QianyeSu/Skyborn/discussions) for user collaboration
- **Documentation**: [Comprehensive Documentation](https://skyborn.readthedocs.io/) with examples and tutorials

## Future Development Roadmap

### Advanced Distribution Features
- **Extended Architecture Support**: ARM64 Linux and additional platforms
- **Distribution Optimization**: Further wheel size reduction and performance enhancement
- **Continuous Performance Monitoring**: Automated benchmarking in CI pipeline
- **Alternative Distribution**: conda-forge ecosystem integration

### Infrastructure Enhancements
- **Advanced Testing**: Extended platform and environment coverage
- **Security Hardening**: Enhanced build security and supply chain protection
- **Performance Analytics**: Real-world usage and performance monitoring

---

*Our comprehensive wheel building and distribution system represents a commitment to providing the atmospheric science community with accessible, high-performance computational tools that meet the highest professional standards.*
