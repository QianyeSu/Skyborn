# Skyborn Version Management System

## Overview

We have implemented a sophisticated **centralized version management system** for the Skyborn project to ensure consistent version tracking across all build configurations and distribution channels. This enterprise-grade approach eliminates version drift and maintains synchronization across complex build environments.

## Architecture and Design

### Single Source of Truth Philosophy

Our version management architecture is built on the principle of maintaining a **single authoritative source** for version information:

- **Primary Source**: `src/skyborn/_version.py` serves as the definitive version repository
- **Automatic Propagation**: `src/skyborn/__init__.py` dynamically imports version information
- **Build Integration**: `pyproject.toml` maintains synchronized version data for build systems
- **Cross-Platform Consistency**: Version information remains consistent across all supported platforms

## Version Update Procedures

### Automated Version Management (Recommended)

We provide a streamlined version update utility for seamless version management:

```bash
python set_version.py 0.3.8
```

**Automated Updates Include:**
- `src/skyborn/_version.py` - Core version definition
- `pyproject.toml` - Build system configuration
- **Validation**: Automatic verification of version synchronization
- **Error Checking**: Comprehensive validation of version format compliance

### Manual Version Management

For advanced users requiring direct control:

1. **Update Primary Source**: Modify `__version__` in `src/skyborn/_version.py`
2. **Synchronize Build System**: Update `version =` field in `pyproject.toml` `[project]` section
3. **Validation**: Verify exact version matching across all configuration files
4. **Testing**: Execute build verification to ensure compatibility

## Version Information Distribution

### Core Configuration Files

| File | Purpose | Update Method |
|------|---------|---------------|
| **`src/skyborn/_version.py`** | Primary version authority | Automated script or manual |
| **`pyproject.toml`** | Build system integration | Synchronized with primary |
| **`src/skyborn/__init__.py`** | Runtime version access | Automatic import from primary |
| **`environment.yml`** | Conda environment specification | Manual coordination |

## Platform and Environment Requirements

### Python Version Policy

Skyborn maintains **Python 3.9+ compatibility** as our minimum supported version. This requirement is consistently enforced across:

- **Build Configuration**: `pyproject.toml` - `requires-python = ">=3.9"`
- **Environment Specification**: `environment.yml` - `python>=3.9`
- **Continuous Integration**: GitHub Actions testing matrix covers Python 3.9, 3.10, 3.11, 3.12
- **Distribution**: PyPI package metadata enforces compatibility requirements

### Cross-Platform Version Consistency

Our version management system ensures identical version reporting across:
- Linux distributions (x86_64)
- macOS (Intel and Apple Silicon)
- Windows (x86_64)
- Container environments and cloud platforms

## Technical Implementation Notes

### Build System Compatibility

**Meson Build System Considerations:**
- Meson lacks native dynamic versioning support
- Our hybrid approach bridges this limitation through coordinated static configuration
- Build-time version validation ensures consistency

### Best Practices

1. **Consistency First**: Always utilize the `set_version.py` script to maintain synchronization
2. **Semantic Versioning**: Adhere to `X.Y.Z` format (MAJOR.MINOR.PATCH)
3. **Pre-Release Testing**: Validate version consistency before distribution
4. **Documentation Updates**: Coordinate version changes with documentation releases

## Quality Assurance

### Version Validation Process

- **Automated Verification**: Script-based validation of version synchronization
- **Build Testing**: Integration testing across all supported platforms
- **Distribution Verification**: PyPI upload validation and testing
- **Documentation Alignment**: Version consistency in all documentation

---

*This version management system reflects our commitment to professional software development practices and ensures reliable, consistent releases for the atmospheric science community.*
