#!/usr/bin/env python3
# Main.py
import sys, os
from pathlib import Path

SCRIPTS = {
    "LWE":  "LWE.py",
    "RLWE": "RLWE.py",
    "MLWE": "MLWE.py",
    "SIS":  "SIS.py",
    "RSIS": "RSIS.py",
    "MSIS": "MSIS.py",
}

def main():
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help", "help"}:
        print("Usage: python3 finder.py <PROBLEM_NAME> [args...]")
        print("Problems:", ", ".join(SCRIPTS))
        print("Example: python3 finder.py LWE -bit 40 -n n -q q -t t -c 'n % 2 == 0'")
        sys.exit(0)

    problem = sys.argv[1].upper()
    script_name = SCRIPTS.get(problem)
    if not script_name:
        print(f"Unknown problem '{sys.argv[1]}'. Choose from: {', '.join(SCRIPTS)}")
        sys.exit(2)

    script_path = Path(__file__).resolve().parent / script_name
    if not script_path.exists():
        print(f"Script not found: {script_path}")
        sys.exit(2)

    # Replace this process with the target script
    os.execv(sys.executable, [sys.executable, str(script_path), *sys.argv[2:]])

if __name__ == "__main__":
    main()
