nextflow.enable.types = true

process fix_t2g_names {
    conda "anaconda::click anaconda::networkx"

    input:
    record(
        t2g_in: Path,
        key: String,
    )

    output:
    record(
        t2g_out: file(t2g_out),
        key: key,
    )

    script:
    t2g_in_name = "${t2g_in}".tokenize('/')[-1]
    t2g_out = "${t2g_in_name}.fixed_names"

    """
    kb_t2g.py --input ${t2g_in} --output ${t2g_out}
    """
}
