"""Parser for MadGraph run card format."""

import re
from pathlib import Path
from typing import Any, Dict


def parse_run_card(file_path: Path) -> Dict[str, Any]:
    """
    Parse a MadGraph run_card.dat file into a dictionary.
    
    The format is: value = variable ! comment
    
    Args:
        file_path: Path to the run_card.dat file
        
    Returns:
        Dictionary with variable names as keys and parsed values
    """
    result = {}
    
    # Pattern to match: optional whitespace, value, =, variable, optional comment
    pattern = re.compile(r'^\s*(.+?)\s*=\s*(\w+)\s*(?:!\s*(.*))?$')
    
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Skip empty lines and comments
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            match = pattern.match(line)
            if match:
                value_str, variable, comment = match.groups()
                value_str = value_str.strip()
                
                # Try to convert value to appropriate type
                value = _parse_value(value_str)
                
                result[variable] = {
                    'value': value,
                    'comment': comment.strip() if comment else None,
                    'line': line_num,
                    'raw_value': value_str
                }
    
    return result


def _parse_value(value_str: str) -> Any:
    """Parse a value string to the appropriate Python type."""
    value_str = value_str.strip()
    
    # Handle boolean values
    if value_str.lower() in ('true', '.true.'):
        return True
    elif value_str.lower() in ('false', '.false.'):
        return False
    
    # Try integer
    try:
        return int(value_str)
    except ValueError:
        pass
    
    # Try float
    try:
        return float(value_str)
    except ValueError:
        pass
    
    # Return as string
    return value_str


def get_simple_dict(parsed_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract just the values from the parsed data for simpler access.
    
    Args:
        parsed_data: Output from parse_run_card()
        
    Returns:
        Dictionary with just variable: value pairs
    """
    return {key: data['value'] for key, data in parsed_data.items()}


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python madgraph_parser.py <run_card.dat>")
        sys.exit(1)
    
    file_path = Path(sys.argv[1])
    parsed = parse_run_card(file_path)
    
    print("Parsed run card parameters:")
    for var, data in parsed.items():
        print(f"{var:20} = {data['value']:15} # {data['comment'] or ''}")
    
    print("\nSimple dictionary:")
    simple = get_simple_dict(parsed)
    for key, value in simple.items():
        print(f"{key}: {value}")