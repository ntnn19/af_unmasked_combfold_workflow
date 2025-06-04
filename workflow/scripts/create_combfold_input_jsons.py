import click
import re
import json
from Bio import SeqIO
from collections import defaultdict
from itertools import product
from pathlib import Path
import string

def parse_stoichiometry(stoich_str):
    """Parse a stoichiometry string like '5A4B3C' into {'A': 5, 'B': 4, 'C': 3}"""
    return {m[1]: int(m[0]) for m in re.findall(r'(\d+)([A-Z])', stoich_str)}

def generate_combinations(max_counts, uniform=False):
    """Generate combinations from 1 to max for each ID"""
    keys = sorted(max_counts)
    ranges = [range(1, max_counts[k] + 1) for k in keys]
    
    for combo in product(*ranges):
        if uniform and len(set(combo)) != 1:
            continue  # skip non-uniform combos
        yield dict(zip(keys, combo))

def chain_name_generator():
    """Generate chain names like A, B, ..., Z, AA, AB, ..."""
    letters = string.ascii_uppercase
    for l in letters:
        yield l
    for i in range(2, 4):  # A-Z, then AA-ZZ, AAA-ZZZ if needed
        for combo in product(letters, repeat=i):
            yield ''.join(combo)

@click.command()
@click.argument("fasta_file", type=click.Path(exists=True))
@click.argument("stoichiometry", type=str)
@click.option("--output-dir", "-d", default="json_output", type=click.Path(file_okay=False), show_default=True)
@click.option("--uniform-copies", "-u", is_flag=True, help="Only generate combinations where all chains have the same number of copies.")
def main(fasta_file, stoichiometry, output_dir, uniform_copies):
    """Generate JSON files up to a given STOICHIOMETRY based on FASTA input"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Parse FASTA and assign letters A, B, C... to the records
    records = list(SeqIO.parse(fasta_file, "fasta"))
    chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    if len(records) > len(chain_letters):
        raise ValueError("Too many sequences in FASTA file for automatic chain labeling (max 26).")

    # Map chain letter to sequence record and ID
    chain_to_record = {chain_letters[i]: r for i, r in enumerate(records)}
    chain_to_id = {k: v.id for k, v in chain_to_record.items()}

    max_stoich = parse_stoichiometry(stoichiometry)
    combos = list(generate_combinations(max_stoich, uniform=uniform_copies))

    for combo in combos:
        json_data = {}
        record_id_counts = defaultdict(int)
        chain_names = chain_name_generator()

        for chain_letter, count in combo.items():
            record = chain_to_record.get(chain_letter)
            if not record:
                raise ValueError(f"No record for chain {chain_letter}")

            record_id = record.id
            assigned_chains = [next(chain_names) for _ in range(count)]
            if record_id not in json_data:
                json_data[record_id] = {
                    "name": record_id,
                    "chain_names": [],
                    "start_res": 1,
                    "sequence": str(record.seq)
                }
            json_data[record_id]["chain_names"].extend(assigned_chains)
            record_id_counts[record_id] += count

        job_name = "_".join(f"{v}{k}" for k, v in sorted(record_id_counts.items()))
        output_file = output_path / f"{job_name}.json"

        with open(output_file, "w") as f:
            json.dump(json_data, f, indent=2)

    click.echo(f"✅ Generated {len(combos)} JSON files in {output_path}")

if __name__ == "__main__":
    main()
