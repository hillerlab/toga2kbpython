<p align="center">
  <span>
    <h1 align="center">
        toga2kbpython
    </h1>
    <h3 align="center">
        Part of TOGA2 suite
    </h3>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/TOGA2" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/TOGA2?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> Gene expression quantification pipeline based on kb-python compatible with TOGA2 output structure </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://github.com/hillerlab/TOGA2/wiki">docs</a> .
        <a href="https://github.com/hillerlab/TOGA2?tab=readme-ov-file#Installation">install</a> .
        <a href ="https://www.biorxiv.org/content/10.64898/2026.06.30.735536v1">preprint</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

>[!WARNING]
This pipeline is currently in early access.

# Description
The pipeline uses TOGA2 output and paired-end RNAseq reads to quantify gene expression with [kb-python](https://github.com/pachterlab/kb_python) (Melsted, P., Booeshaghi, A.S., et al., 2021). The pipeline consists of the following steps:
1. Each unique genome assembly file is indexed with [samtools faidx](https://www.htslib.org/doc/samtools-faidx.html)
2. Each unique TOGA2 output annotation directory is used to create a `kb-python`-compatible GTF file using [bed2gtf](https://github.com/alejandrogzi/bed2gtf).
2. For each unique genome-TOGA2 output combination, the pipeline produces a decoy file.
3. Input genome files are indexed with `kb ref` for each unique decoy file (genome-TOGA2 combination).
4. TOGA2 genes from the GTF files are renamed for each unique TOGA2 output directory, with many:1/many:many genes collapsed into single units (see "Caveats" for more information).
5. For each RNAseq sample in each species, pairwise RNAseq reads are pseudoligned to the indexed genome using the index and the decoy files obtained at the previous steps.
6. Genes in all output directories are filtered and assigned the new names as established at step 4.

# Running toga2kbpython
Starting from `v0.2`, `toga2kbpython` requires only one mandatory argument, `--table`, standing for the input tab-separated file describing input for every unique `kb-python` run. See [example file](https://github.com/hillerlab/toga2kbpython/tree/main/assets/sample_input/table.tsv) and [table file README](https://github.com/hillerlab/toga2kbpython/tree/main/assets/sample_input/README.md) for additional clarification.  
You can specify output directory with `--output` parameter. The `--include_utr` parameter controls using TOGA2-predicted untranslated region sequences in the pseudoalignment step (see *Caveats* section below).

# Caveats
> [!NOTE]
> *
> * This pipeline uses Conda for dependency resolution.
> * The following conventions are applied to TOGA2-inferred genes and their names for the purpose of `kb-python` compatibility:
>    * many:many genes are treated as single transcription units. For example, if a many:many orthology group contains five genes in the query, they are represented by a single gene entry in the output files, with transcription count/transcript per million values representing the sum of respective values for all five genes.
>    * Likewise, for query genes containing fragmented projections, their transcription values are summed across all fragments.
>    * To facilitate cross-species comparison, all genes are assigned expanded gene names containing all refeference gene symbols mapping to the respective locus (or loci, in case of many:many orthologs). See `toga.gene_names.txt` for kb-python-to-TOGA2 gene mapping.
>    * **Important**: The current implementation (v1.0) relies on the Hiller Lab transcript naming convention (TRANSCRIPT#GENE for reference transcripts, TRANSCRIPT#GENE#CHAIN for query projections) to create gene naming mapping.
> * It is highly recommended **not** to include the untranslated regions (UTRs) in the input annotation since a) the UTR projections from TOGA2 may be less accurate and vary between tissues, and b) the general aim of `toga2kbpython` is to compare species rather than to measure the exact read counts. As of such, the pipeline uses the bare CDS annotation version from TOGA2 (`query_annotation.bed`) by default. 
>    * If you want to include the UTR sequence into your annotation, add the `--include_utr` flag or set the respective config parameter to `true`. In this case, the entire input genome is used as decoy, with no CDS masking.
> * If UTR is excluded from the annotation (current default behavior), the CDS-masked genome fasta is passed to `kb-python` as the decoy file (`--d-list`). This makes the decoy and the transcript mutually exclusive, so that reads that map to both UTR and CDS will be counted (otherwise they will be excluded by decoy).
> * The pipeline (currently) excludes the retrogene candidates from the input annotation. While retrogene candidates are expected to be functionally intact by default (loss status of *FI* or *I* being a prerequisite for a processed pseudogene to be classified as a retrogene candidate), the false positive retrogene predictions may draw reads mapping to the orthologous loci, resulting in underestimation of their expression.
>    * Currently there is **no** option to keep the retrogene candidates in the input annotation; this may change in the future.
> * The pipeline **does not** filter the query genes based on their loss status and/or orthology class.
