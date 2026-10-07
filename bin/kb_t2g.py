#!/usr/bin/env python3

"""
"""


__author__ = "Yury V. Malovichko"
__credits__ = ("Xueling Yi", )

from typing import List, Set

import click
import networkx as nx

NOT_FOUND: str = "GENE_NOT_FOUND"
FRAGMENT: str = "_FRAGMENT."
FROM: str = "#f#"
TO: str = "#t#"
HEADER: str = "transcript\tgene_id\tTOGA_ref\n"


def get_connected_components(graph: nx.Graph) -> List[nx.Graph]:
    """ """
    ## check the NetworkX version
    nx_v: str = nx.__version__
    v_split: List[str] = [x for x in nx_v.split(".") if x.isnumeric()]
    if len(v_split) > 1:
        NX_VERSION: float = float(f"{v_split[0]}.{v_split[1]}")
    else:
        NX_VERSION: float = float(v_split[0])
    if NX_VERSION < 2.4:
        raw_components = list(nx.connected_component_subgraphs(graph))
    else:
        raw_components = [graph.subgraph(c) for c in nx.connected_components(graph)]
    return raw_components


@click.command(no_args_is_help=True)
@click.option(
    "--input",
    type=click.Path(exists=True),
    required=True,
    default=None,
    help="Path to the input t2g file from kb-ref output"
)
@click.option(
    "--output",
    type=click.Path(exists=False),
    required=True,
    default=None,
    help="Path to the output file"
)
def main(input: click.Path, output: click.Path) -> None:
    """
    \b
    For a reference assembly used in kb-python, generate the file that matches 
    each gene ID with the TOGA2-output gene symbols based on the kb-python index t2g.txt, 
    with gene symbols merged for fragmented and/or duplicated genes.

    \b
    WARNING: As of now, the script relies on the isoforms to have a gene name as the second sharp-delimited field 
    (e.g., TRANSCRIPT_ID#GENE_ID), as represented in the official Hiller Lab reference annotations and prepare-input 
    output files if GTF/GFF3-format reference annotation was used.
    """
    # tr2gene: Dict[str, str] = {}
    # gene2source: Dict[str, Set[str]] = defaultdict(set)
    graph: nx.Graph = nx.Graph()
    with open(input, "r") as h:
        for i, line in enumerate(h, start=1):
            data: List[str] = line.strip().split("\t")
            if not data or not data[0]:
                continue
            if len(data) < 2:
                raise Exception(
                    "Wrong file formatting at line %i: expected at least 2 fields, got %i"
                    % (i, len(data))
                )
            tr, gene = data[:2]
            ## ignore the GENE_NOT_FOUND items
            if gene == NOT_FOUND:
                continue
            ## collapse the fragments into a single gene
            gene = gene.split(FRAGMENT)[0]
            ## retrieve the gene ID
            ## WARNING: See the caveat in the help message
            gene_symbol: str = tr.split("#")[1]
            # gene2source[gene].add(gene_symbol)
            ## add the prefices
            gene = FROM + gene
            gene_symbol = TO + gene_symbol
            graph.add_edge(gene, gene_symbol)
    ## write the fixed output names
    with open(output, "w") as h:
        h.write(HEADER)
        for component in get_connected_components(graph):
            all_from: Set[str] = {x[len(FROM):] for x in component.nodes() if x.startswith(FROM)}
            all_to: List[str] = sorted({x[len(TO):] for x in component.nodes() if x.startswith(TO)})
            ref_gene_num: int = len(all_to)
            combined_symbol: str = ",".join(all_to)
            for toga_gene in all_from:
                h.write(f"{combined_symbol}\t{toga_gene}\t{ref_gene_num}\n")


if __name__ == "__main__":
    main()