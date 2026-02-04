#!/usr/bin/env python3
import logging
import json
from typing import Set, Tuple

from Bio import SeqIO
import click
import pandas as pd

from filter_reformat import read_strain_accessions

REFERENCE_ACCESSION = "NC_002031.1"

logger = logging.getLogger(__name__)
logging.basicConfig(
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)8s (%(filename)20s:%(lineno)4d) - %(message)s ",
    datefmt="%H:%M:%S",
)


def load_clades(clade_json: str) -> pd.DataFrame:
    with open(clade_json, "r") as f:
        data = json.load(f)

    df_clade = {
        "Accession": [],
        "Clade": [],
    }
    for k, v in data["nodes"].items():
        if k.startswith("NODE_"):
            continue
        df_clade["Accession"].append(k)
        df_clade["Clade"].append(v["clade_membership"])

    return pd.DataFrame(df_clade)


@click.command(
    help="""
            Generate an example dataset by taking a random subsample of `-n` sequences from `--input-fasta`. The
            sample will be weighted absed on the size of clades in `--clade-membership` so that each clade should 
            be represented by roughly the same amount of sequences in the output sample.

            Optionally, you can specify files containing manually selected sequence accessions to keep using 
            the `--bad-strains` and `--mediocre-strains` arguments. Sequences listed in this file are added to
            the randomly selected sample.

            The script errors if `n` is larger than the number of sequences in the dataset.
         """
)
@click.option("--input-fasta", required=True, type=click.Path(exists=True))
@click.option("-n", required=True, type=int)
@click.option("--clade-membership", required=True, type=click.Path(exists=True))
@click.option("--seed", default=42, type=int)
@click.option("--output", required=True, type=str)
@click.option("--bad-strains", required=False, type=click.Path(exists=True))
@click.option("--mediocre-strains", required=False, type=click.Path(exists=True))
@click.option(
    "--log-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
)
def main(
    input_fasta: str,
    n: int,
    clade_membership: str,
    seed: int,
    output: str,
    bad_strains: str,
    mediocre_strains: str,
    log_level: str,
):
    logger.setLevel(log_level)

    df_clades = load_clades(clade_membership)
    if df_clades.shape[0] < n:
        raise ValueError(
            f"cannot sample {n} sequences when there are only {df_clades.shape[0]} listed in {clade_membership}"
        )

    # DataFrame.sample() automatically normalizes weights, so don't need to sum to 1.0
    weights = 1 / df_clades.groupby("Clade")["Clade"].transform("size")
    sample: set[str] = set(
        df_clades.sample(n, weights=weights, random_state=seed)["Accession"].to_list()
    )

    if mediocre_strains:
        mediocre = set(read_strain_accessions(mediocre_strains))
        sample = sample.union(set(mediocre))
    if bad_strains:
        bad = set(read_strain_accessions(bad_strains))
        sample = sample.union(set(bad))

    sample.remove(REFERENCE_ACCESSION)

    with open(output, "w") as o:
        with open(input_fasta):
            records = SeqIO.parse(input_fasta, "fasta")
            to_keep = [record for record in records if record.id in sample]
            SeqIO.write(to_keep, o, "fasta")
    logger.info(f"Sampled {len(to_keep)} sequences and wrote them to {output}")


if __name__ == "__main__":
    main()
