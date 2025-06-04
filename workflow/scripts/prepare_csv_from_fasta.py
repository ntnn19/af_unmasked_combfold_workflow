import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import combinations_with_replacement
from collections import Counter
import pandas as pd
from pathlib import Path
import re

def get_template_names(templates_dir):
    templates_path = Path(templates_dir)
    if not templates_path.exists():
        return ["no_template"]
    return sorted([
        re.sub(r'\W+', '_', f.stem)
        for f in templates_path.glob("*.pdb")
    ]) or ["no_template"]

@click.command()
@click.argument("fasta_file", type=click.Path(exists=True))
@click.option("--min-size", default=2, show_default=True, help="Minimum combination size")
@click.option("--max-size", default=6, show_default=True, help="Maximum combination size")
@click.option("--output", "-o", default="unique_combinations.tsv", show_default=True, help="Base name for TSV output file")
@click.option("--output-dir", "-d", default=".", type=click.Path(file_okay=False), show_default=True, help="Top-level output directory")
@click.option("--templates-dir", "-t", default="templates", show_default=True, help="Directory containing template PDB files")
def generate_combinations(fasta_file, min_size, max_size, output, output_dir, templates_dir):
    """Generate FASTA files for unique sequence combinations per template, saved in template-named subdirectories."""

    base_output_path = Path(output_dir)
    base_output_path.mkdir(parents=True, exist_ok=True)
    template_names = get_template_names(templates_dir)

    records = list(SeqIO.parse(fasta_file, "fasta"))
    id_to_record = {r.id: r for r in records}
    ids = list(id_to_record)

    unique_combos = set()
    for k in range(min_size, max_size + 1):
        for combo in combinations_with_replacement(ids, k):
            count = Counter(combo)
            signature = tuple(sorted(count.items()))
            unique_combos.add(signature)

    unique_combos = sorted(unique_combos)

    for template in template_names:
        template_dir = base_output_path / template
        template_dir.mkdir(parents=True, exist_ok=True)

        rows = []
        for combo in unique_combos:
            base_job_name = "_".join(f"{v}{seq_id}" for seq_id, v in combo)
            job_name = f"{base_job_name}_{template}"

            seq_records = []
            for seq_id, count in combo:
                record = id_to_record[seq_id]
                for i in range(count):
                    record_id = f"{job_name}_{seq_id}_{i+1}" if count > 1 else f"{job_name}_{seq_id}"
                    rows.append({
                        "job_name": job_name,
                        "type": "protein",
                        "id": record_id,
                        "sequence": str(record.seq),
                    })
                    seq_records.append(SeqRecord(
                        seq=record.seq,
                        id=record_id,
                        description=record_id
                    ))

            # Save FASTA in this template subdirectory
            fasta_path = template_dir / f"{job_name}.fasta"
            SeqIO.write(seq_records, fasta_path, "fasta")

        # Save TSV in the same subdirectory
        df = pd.DataFrame(rows)
        df_out = template_dir / f"{Path(output).stem}_{template}.tsv"
        df.to_csv(df_out, sep="\t", index=False)

    click.echo(f"✅ Wrote {len(unique_combos)} combinations per template to subfolders in: {base_output_path}")

if __name__ == "__main__":
    generate_combinations()
