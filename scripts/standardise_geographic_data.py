#!/usr/bin/env python3
import logging

import click
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)8s (%(filename)20s:%(lineno)4d) - %(message)s ",
    datefmt="%H:%M:%S",
)


@click.command(
    help="""
            Rename the 'Geographic Region' column to 'Region'.

            Also strip region specification from `--location-with-region-col` entries in the metadata file and add them as a new `Country` column.
            The reason for this is that the 'Geographic Location' column in NCBI metadata may or may not contain a region specification, e.g.,
            Brazil, Brazil: Brazil: Afua, Par, Brazil: Amorinopolis, GO, etc. 

            To add less fine-grained, country-level information, we extract the part before the ':' and add this as a new column. It is recommended to
            check whether this approach is valid for your current metadata file before running this script.
         """
)
@click.option("--metadata", required=True, type=click.Path(exists=True))
@click.option("--location-with-region-col", required=True, type=str)
@click.option("--output", required=True, type=str)
@click.option(
    "--log-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
)
def main(
    metadata: str,
    location_with_region_col: str,
    output: str,
    log_level: str,
):
    logger.setLevel(log_level)

    df = pd.read_csv(metadata, sep="\t", encoding="utf-8").rename(
        columns={"Geographic Region": "Region"}
    )
    df["Country"] = pd.Series(index=range(0, df.shape[0]), dtype=object)

    expression = {
        "Country": lambda x: [i.split(":")[0] for i in x[location_with_region_col]]
    }
    df[~df[location_with_region_col].isna()] = df[
        ~df[location_with_region_col].isna()
    ].assign(**expression)

    df.to_csv(output, sep="\t", index=False)

    logger.info(
        f"Read metadata file at {metadata}, added column 'Country' with country information, saved updated metadata to {output}"
    )


if __name__ == "__main__":
    main()
