#!/usr/bin/env python3
"""
Batch rename AB1 files to match pipeline format
"""
import os
import re
import shutil
from pathlib import Path

def clean_filename(old_name):
    """Convert filename to pipeline format"""
    # Remove .ab1 extension for processing
    name = old_name.replace('.ab1', '')
    
    # Extract the main parts using the pattern in your files
    # Pattern: "Sample info_RBCLA F/R conrad'_code.ab1"
    
    # Replace RBCLA with RBCL (standardize)
    name = name.replace('RBCLA', 'RBCL')
    
    # Find F or R for direction
    if ' F ' in name:
        direction = 'F'
        parts = name.split(' F ')
    elif ' R ' in name:
        direction = 'R'
        parts = name.split(' R ')
    else:
        return None
    
    if len(parts) >= 1:
        # Get sample name and clean it
        sample = parts[0].strip()
        # Remove spaces and special characters
        sample = sample.replace(' ', '-')
        sample = sample.replace("'", "")
        
        # Create new filename
        new_name = f"{sample}_RBCL_{direction}.ab1"
        return new_name
    
    return None

def main():
    input_dir = Path(os.path.expanduser("~/Analysis/ab1_files"))
    
    # Create backup directory
    backup_dir = input_dir / "original_names_backup"
    backup_dir.mkdir(exist_ok=True)
    
    files = list(input_dir.glob("*.ab1"))
    
    print(f"Found {len(files)} AB1 files")
    print("=" * 70)
    
    renamed_count = 0
    
    for old_path in files:
        new_name = clean_filename(old_path.name)
        
        if new_name:
            new_path = input_dir / new_name
            
            print(f"OLD: {old_path.name}")
            print(f"NEW: {new_name}")
            
            # Copy to backup first
            backup_path = backup_dir / old_path.name
            shutil.copy2(old_path, backup_path)
            
            # Rename the file
            old_path.rename(new_path)
            renamed_count += 1
            print("✓ Renamed successfully")
            print("-" * 70)
    
    print(f"\nSummary:")
    print(f"- Renamed {renamed_count} files")
    print(f"- Original files backed up in: {backup_dir}")
    print(f"\nYour files are now ready for the pipeline!")

if __name__ == "__main__":
    main()
