#!/usr/bin/env python3
"""
Dependency checker for Sanger sequencing pipeline
Checks all required software and Python packages
"""

import sys
import subprocess
import shutil
import platform

def check_command(cmd):
    """Check if a command exists in PATH"""
    return shutil.which(cmd) is not None

def get_version(cmd, version_flag="--version", line_num=0, stderr=False):
    """Get version string from a command"""
    try:
        result = subprocess.run(
            [cmd, version_flag], 
            capture_output=True, 
            text=True, 
            timeout=5
        )
        output = result.stderr if stderr else result.stdout
        lines = output.strip().split('\n')
        return lines[line_num] if lines else "Unknown version"
    except:
        return "Error getting version"

print("=" * 60)
print("Sanger Pipeline Dependency Checker")
print("=" * 60)
print()

# System information
print("SYSTEM INFORMATION")
print("-" * 40)
print(f"Operating System: {platform.system()} {platform.release()}")
print(f"Machine: {platform.machine()}")
print(f"Python executable: {sys.executable}")
print(f"Python version: {sys.version.split()[0]}")
print()

# Check Python version
print("PYTHON VERSION CHECK")
print("-" * 40)
python_version = sys.version_info
if python_version >= (3, 7):
    print(f"✓ Python {python_version.major}.{python_version.minor}.{python_version.micro} - OK")
else:
    print(f"✗ Python {python_version.major}.{python_version.minor}.{python_version.micro} - Too old (need 3.7+)")
print()

# Check Python packages
print("PYTHON PACKAGES")
print("-" * 40)

required_packages = {
    "Bio": "biopython",
    "pandas": "pandas",
    "numpy": "numpy"
}

all_packages_ok = True
for module_name, package_name in required_packages.items():
    try:
        module = __import__(module_name)
        version = getattr(module, "__version__", "Unknown version")
        print(f"✓ {package_name}: {version}")
    except ImportError:
        print(f"✗ {package_name}: NOT INSTALLED")
        print(f"  Install with: pip install {package_name}")
        all_packages_ok = False

print()

# Check external tools
print("EXTERNAL TOOLS")
print("-" * 40)

# Check MAFFT
if check_command("mafft"):
    version = get_version("mafft", "--version", stderr=True)
    print(f"✓ MAFFT: {version}")
else:
    print("✗ MAFFT: NOT FOUND")
    print("  Install with: brew install mafft")

# Check IQ-TREE (both versions)
iqtree_found = False
if check_command("iqtree2"):
    version = get_version("iqtree2", "--version")
    print(f"✓ IQ-TREE2: {version}")
    iqtree_found = True
elif check_command("iqtree"):
    # Test if iqtree works (might have library issues)
    try:
        result = subprocess.run(
            ["iqtree", "--version"], 
            capture_output=True, 
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            version = result.stdout.strip().split('\n')[0]
            print(f"✓ IQ-TREE: {version}")
            iqtree_found = True
        else:
            print("✗ IQ-TREE: Installed but not working")
            print(f"  Error: {result.stderr}")
            if "libgomp" in result.stderr:
                print("  Fix: brew install gcc")
    except:
        print("✗ IQ-TREE: Installed but not working")
else:
    print("✗ IQ-TREE: NOT FOUND (optional)")
    print("  Install with: brew tap brewsci/bio && brew install brewsci/bio/iqtree")
    print("  Or download from: http://www.iqtree.org")

print()

# Check pipeline script
print("PIPELINE SCRIPT")
print("-" * 40)
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
pipeline_script = os.path.join(script_dir, "sanger_pipeline.py")

if os.path.exists(pipeline_script):
    size = os.path.getsize(pipeline_script)
    print(f"✓ sanger_pipeline.py found ({size:,} bytes)")
    
    # Check if it's executable
    if os.access(pipeline_script, os.X_OK):
        print("✓ Script is executable")
    else:
        print("! Script is not executable")
        print(f"  Fix with: chmod +x {pipeline_script}")
else:
    print("✗ sanger_pipeline.py NOT FOUND")

print()

# Summary
print("SUMMARY")
print("-" * 40)

can_run_basic = all_packages_ok and check_command("mafft")
can_run_full = can_run_basic and iqtree_found

if can_run_full:
    print("✓ All dependencies satisfied - pipeline fully functional")
elif can_run_basic:
    print("✓ Basic dependencies satisfied - pipeline functional (no tree building)")
    print("  Use --no_tree flag when running the pipeline")
else:
    print("✗ Missing critical dependencies - please install required packages")

print()
print("=" * 60)

# Test import of pipeline
if os.path.exists(pipeline_script):
    print("\nTesting pipeline import...")
    try:
        # Add script directory to path
        sys.path.insert(0, script_dir)
        import sanger_pipeline
        print("✓ Pipeline script imports successfully")
    except Exception as e:
        print(f"✗ Error importing pipeline: {e}")

print("\nDependency check complete!")
