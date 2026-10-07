#!/usr/bin/env python3

"""
Creates a provisional isoform file for BED-to-GTF conversion
"""

from collections import defaultdict
from typing import Any, Dict, List

import click

CONTEXT_SETTINGS: Dict[str, Any] = {
    "help_option_names": ["-h", "-help", "--help"],
    "ignore_unknown_options": True,
    "allow_extra_args": True,
    "max_content_width": 150,
}

HEADER: str = "query_gene"


def segment_base(projection: str) -> str:
    """Returns the fragmented projection's name, ignoring the fragment number"""
    return projection.split("$")[0]


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.argument("isoform_file", type=click.Path(exists=True), metavar="IN_ISOFORM_FILE")
@click.argument("bed_file", type=click.Path(exists=True), metavar="FINAL_BED_FILE")
@click.argument("output", type=click.Path(exists=False), metavar="OUT_ISOFORM_FILE")

def main(isoform_file: click.Path, bed_file: click.Path, output: click.Path) -> None:
    ## record the fragmented projections' names
        fragmented_names: Dict[str, int] = defaultdict(int)
        with open(bed_file, "r") as h:
            for i, line in enumerate(h, start=1):
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                if len(data) != 12:
                    raise Exception(
                        (
                            "Improper formatting at BED file line %i; "
                            "expected 12 fields, got %i"
                        )
                        % (i, len(data))
                    )
                name: str = data[3]
                if "$" not in name:
                    continue
                basename: str = segment_base(name)
                fragmented_names[basename] += 1

        ## read the isoform file, check each line in it;
        ## fragmented projection mapping is to be expanded to accommodate for each fragment,
        ## the rest of the lines are written as-is
        with open(isoform_file, "r") as ih, open(output, "w") as oh:
            for i, line in enumerate(ih, start=1):
                data: List[str] = line.strip().split("\t")
                if not data or not data[0]:
                    continue
                if len(data) != 2:
                    raise Exception(
                        (
                            "Improper formatting at isoform file line %i; "
                            "expected 2 fields, got %i"
                        )
                        % (i, len(data))
                    )
                gene, proj = data
                if gene == HEADER:
                    continue
                if proj not in fragmented_names:
                    oh.write(line)
                    continue
                for num in range(1, fragmented_names[proj] + 1):
                    fragm_name: str = f"{proj}${num}"
                    oh.write(gene + "\t" + fragm_name + "\n")


if __name__ == "__main__":
    main()