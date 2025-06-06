import re
import sys
import argparse

def parse_stoichiometry(basename):
    # Extract pattern like 2A, 3B, etc. allowing arbitrary text between
    parts = re.findall(r'(\d+)([A-Za-z]+)', basename)
    output1 = []
    output2 = []
    chain_ord = ord("A")
    for count_str, label in parts:
        count = int(count_str)
        output1.extend([label] * count)
        output2.extend([chr(chain_ord + i) for i in range(count)])
        chain_ord += count
    return output1, output2

def main():
    parser = argparse.ArgumentParser(description="Generate chain identifiers from a stoichiometry string")
    parser.add_argument("basename", help="Input like '1A1B2C'")

    args = parser.parse_args()

    output1, output2 = parse_stoichiometry(args.basename)
    print(" ".join(output1) + "-" +  " ".join(output2))
#    print(" ".join(output2))

if __name__ == "__main__":
    main()
