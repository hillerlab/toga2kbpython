#!/usr/bin/env python3

"""
post-kbpython cleanup: 
1) Merge the quantifications of fragmented/duplicated genes. 
2) Change gene IDs to TOGA gene symbols.
The output cleaned files can be used in inter-species comparisons: TPM for visualization, counts for DESeq2.
"""


__author__ = "Yury V. Malovichko"
__credits__ = ("Xueling Yi",)

from typing import Dict, List

import click
import os

ABUNDANCE: str = "abundance.gene.tsv"
NOT_FOUND: str = "GENE_NOT_FOUND"
GENE_ID: str = "gene_id"
COUNT_FILE: str = "toga.kb.counts.tsv"
TMP_FILE: str = "toga.kb.tpm.tsv"
COUNT_HEADER: str = "gene\test_counts\n"
TPM_HEADER: str = "gene\ttpm\n"
OUT_LINE: str = "{gene}\t{value}\n"

@click.command(no_args_is_help=True)
@click.option(
    "--input_dir",
    "-i",
    type=click.Path(exists=True),
    metavar="PATH",
    required=True,
    default=None,
    show_default=True,
    help="Path to the `kb-python count` output",
)
@click.option(
    "--t2g",
    "-t",
    type=click.Path(exists=True),
    metavar="PATH",
    required=True,
    default=None,
    show_default=True,
    help="Path to the (corrected) transcript2gene mapping file",
)
@click.option(
    "--output_dir",
    "-o",
    type=click.Path(exists=False),
    metavar="PATH",
    required=True,
    default=None,
    show_default=True,
    help="Path to the output directory",
)

def main(input_dir: click.Path, t2g: click.Path, output_dir: click.Path) -> None:
    ## read the transcript-to-gene mapping
    g2t: Dict[str, str] = {}
    with open(t2g, "r") as h:
        for i, line in enumerate(h, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 3:
                raise Exception(
                    (
                        "Improper t2g file formatting at line %i: expected three fields, got %i"
                    )
                    % (i, len(data))
                )
            g2t[data[1]] = data[0]
    ## read the output abundance file
    abundance_file: str = os.path.join(input_dir, ABUNDANCE)
    if not os.path.exists(abundance_file):
        raise FileNotFoundError("Abundance file %s does not exist" % abundance_file)
    gene2abundance: Dict[str, float] = {}
    gene2tmp: Dict[str, float] = {}
    with open(abundance_file, "r") as h:
        for i, line in enumerate(h, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 3:
                raise Exception(
                    (
                        "Improper abundance file formatting at line %i: expected three fields, got %i"
                    ) % (i, len(data))
                )
            if data[0] == GENE_ID:
                continue
            ## field 0 - gene_id, field 1 - gene_name, field 2 - est_counts, field 3 - TMP
            gene: str = data[0]
            try:
                est_count: float = float(data[2])
                tpm: float = float(data[3])
            except ValueError:
                raise Exception(
                    (
                        "Improper abundance file formatting at line %i: "
                        "non-numeric values in `est_count` or `tpm` fields"
                    ) % i
                )
            if gene in gene2abundance:
                gene2abundance[gene] += est_count
                gene2tmp[gene] += tpm
            else:
                gene2abundance[gene] = est_count
                gene2tmp[gene] = tpm
    if not gene2abundance:
        raise Exception(
            "No genes found in the output abundance file"
        )
    ## sum up counts and tpm of all genes whose IDs mapped to the same gene symbol (fragmented, one2many, many2many)
    ## many2one is already one unique quantification

    ## write the output files
    os.makedirs(output_dir, exist_ok=True)
    count_file: str = os.path.join(output_dir, COUNT_FILE)
    tpm_file: str = os.path.join(output_dir, TMP_FILE)
    with open(count_file, "w") as ch, open(tpm_file, "w") as th:
        ch.write(COUNT_HEADER)
        th.write(TPM_HEADER)
        for gene in gene2abundance.keys():
            new_name: str = g2t.get(gene)
            if new_name is None:
                raise ValueError("Gene %s does not have a corresponding updated name" % gene)
            count_line: str = OUT_LINE.format(gene=new_name, value=gene2abundance[gene])
            tpm_line: str = OUT_LINE.format(gene=new_name, value=gene2tmp[gene])
            ch.write(count_line)
            th.write(tpm_line)


if __name__ == "__main__":
    main()