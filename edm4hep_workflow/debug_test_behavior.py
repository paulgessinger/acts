#!/usr/bin/env python3

import tempfile
import shutil
from pathlib import Path
from src.colliderml.cli.madgraph import parse_pythia8_card, update_pythia8_card

def debug_test_behavior():
    """Debug what actually happens in the test."""
    
    # Create test content similar to real file
    test_content = """# Pythia8 test
Main:numberOfEvents = -1
HEPMCoutput:file = hepmc.gz  
JetMatching:qCut = -1.0
JetMatching:doShowerKt = off
Merging:Process = <set_by_user>
SysCalc:fullCutVariation = off
"""
    
    # Test with both old and new regex patterns
    patterns = [
        r"^(\s*?)(\w+)(\s*=\s*)(.*?)(\s?#\s.*)",      # Old (broken)
        r"^(\s*?)([\w:]+)(\s*=\s*)(.*?)(\s?#\s.*)"    # New (fixed)
    ]
    
    for i, pattern in enumerate(patterns, 1):
        print(f"\n=== Testing Pattern {i}: {'OLD (broken)' if i == 1 else 'NEW (fixed)'} ===")
        print(f"Pattern: {pattern}")
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            test_file = tmp_path / "test_pythia8.dat"
            test_file.write_text(test_content)
            
            print(f"\nOriginal file:")
            print(test_file.read_text())
            
            # Parse original
            original = parse_pythia8_card(test_file)
            print(f"Original parsed: {original}")
            
            # Updates including colon keys
            updates = {
                "Main:numberOfEvents": "1000",      # Should UPDATE existing
                "JetMatching:qCut": "20.0",         # Should UPDATE existing  
                "Beams:eCM": "13000.0",             # Should ADD new
            }
            
            print(f"\nApplying updates: {updates}")
            
            # Temporarily patch the regex pattern to test the old behavior
            import re
            from src.colliderml.cli import madgraph
            
            # Save original function
            original_func = madgraph.update_pythia8_card
            
            def patched_update(file_path, updates_dict, log_base):
                raw = file_path.read_text()
                existing = set()
                
                def repl(m):
                    prefix, key, infix, val, suffix = m.groups()
                    existing.add(key)
                    if key in updates_dict:
                        print(f"  UPDATING: {key} = {val} -> {updates_dict[key]}")
                        val = updates_dict[key]
                    return prefix + key + infix + str(val) + suffix
                
                updated = re.sub(pattern, repl, raw, flags=re.MULTILINE)
                
                for key in set(updates_dict.keys()) - existing:
                    print(f"  ADDING: {key} = {updates_dict[key]}")
                    updated += f"\n{key} = {updates_dict[key]} # added after the fact\n"
                
                file_path.write_text(updated)
            
            # Use patched function
            patched_update(test_file, updates, tmp_path)
            
            print(f"\nUpdated file:")
            updated_content = test_file.read_text()
            print(updated_content)
            
            # Parse updated
            updated = parse_pythia8_card(test_file)
            print(f"Updated parsed: {updated}")
            
            # Check if values are correct
            print(f"\nValidation:")
            for key, expected in updates.items():
                actual = updated.get(key, "<MISSING>")
                status = "✅" if actual == expected else "❌"
                print(f"  {status} {key}: expected='{expected}', actual='{actual}'")

if __name__ == "__main__":
    debug_test_behavior()