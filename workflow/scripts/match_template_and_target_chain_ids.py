import re
import argparse
import string

def parse_stoichiometry(basename):
    parts = re.findall(r'(\d+)([A-Za-z]+)', basename)

    output1 = [label for count, label in parts for _ in range(int(count))]
    used_chars = set("".join(output1))
    available_chars = [c for c in string.ascii_letters if c not in used_chars]

    total_needed = sum(int(count) for count, _ in parts)
    if len(available_chars) < total_needed:
        raise ValueError("Not enough unique characters for output2 without overlap.")

    output2 = available_chars[:total_needed]
    return output1, output2

def main():
    parser = argparse.ArgumentParser(description="Generate chain identifiers from a stoichiometry string")
    parser.add_argument("basename", help="Input like '1A1B2C'")
    args = parser.parse_args()

    output1, output2 = parse_stoichiometry(args.basename)
    print(" ".join(output1) + " - " + " ".join(output2))

if __name__ == "__main__":
    main()

