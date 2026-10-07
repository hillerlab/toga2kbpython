nextflow.enable.types = true

process create_gtf_for_kbpython {

    conda "anaconda::click bioconda::bed2gtf=1.9.3"

    input:
    toga_dir: Path
    include_utr: Boolean

    output:
    record(
        gtf: file(gtf_file),
        _toga_dir: toga_dir_key,
    )

    script:
    toga_dir_key = "${toga_dir}".tokenize('/')[-1]
    bed_file = include_utr ? "query_annotation.with_utrs.bed" : "query_annotation.bed"
    gtf_file = "${toga_dir_key}.query_annotation.gtf"

    """
    ## step 1: create a provisional isoform file
    isoforms_for_gtf.py ${toga_dir}/query_genes.tsv ${toga_dir}/${bed_file} prov_isoforms.tsv

    ## step 2: create a GTF file
    bed2gtf -b ${toga_dir}/${bed_file} -i prov_isoforms.tsv -o ${gtf_file}

    ## step 3: remove retrogenes as suggested
    sed -i '/#retro/d; /retro_/d' ${gtf_file}

    ## step 3: clean the provisional file
    rm prov_isoforms.tsv
    """
}
