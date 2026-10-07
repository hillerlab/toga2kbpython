/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.

Author: Yury V. Malovichko
Credits: Xueling Yi, Alejandro Gonzales-Iribarren
Version: v0.3
*/

include { process_input_line } from "./modules/types.nf"
include { twoBitToFa } from "./modules/twoBitToFa.nf"
include { samtools_faidx } from "./modules/samtools_faidx.nf"
include { create_gtf_for_kbpython } from "./modules/create_gtf.nf"
include { create_decoy } from "./modules/create_decoy.nf"
include { kb_ref } from "./modules/kb_ref.nf"
include { fix_t2g_names } from "./modules/fix_t2g_names.nf"
include { kb_count } from "./modules/kb_count.nf"
include { post_count } from "./modules/post_kbcount.nf"


params {
    table: Path
    out_dir: String = "output"
    include_utr = false
    strand: String = "unstranded"
    pairwise = false
}

workflow{

    main:

    input_channel = channel
        .fromPath(params.table)
        .splitCsv(sep: '\t')
        .map{ line -> process_input_line(line) }

    // part 1: kb_index_final.sh content
    prep_channel = input_channel
        .map{line -> record(key: "${line.genome}.${line.toga_dir}", genome: line.genome, toga_dir: line.toga_dir) }
        .unique()

    // convert the 2bit genome to Fasta format
    twoBitToFa(prep_channel.map{ x -> x.genome }.unique()) // emits record(fasta, _genome)

    // create a Samtools index
    samtools_faidx(twoBitToFa.out.map {x -> x.fasta} ) // does not need any binding keys

    // create a GTF file for UTR-less annotation
    create_gtf_for_kbpython(
        prep_channel.map{x -> x.toga_dir},
        params.include_utr,
    ) // emits Gtf output and toga_dir key

    // combine indexed Fasta files and the progenitor channel into a single channel
    decoy_channel = prep_channel
        .map {x -> record(genome: x.genome, toga_dir: x.toga_dir, _genome: file(x.genome).name)}
        .join(twoBitToFa.out, by: '_genome')

    // create a decoy file
    create_decoy(
        decoy_channel.map {x -> record(fasta: x.fasta, toga_dir: x.toga_dir)},
        params.include_utr
    ) // emits a decoy Fasta file and a (genome, TOGA2) key

    // join the results of all previous steps
    kbref_channel = create_decoy
        .out
        .join(twoBitToFa.out, by: '_genome')
        .join(create_gtf_for_kbpython.out, by: '_toga_dir')
        .map {
            x -> record(
                fasta: x.fasta,
                gtf: x.gtf,
                decoy: x.decoy,
                key: "${x._genome}.${x._toga_dir}"
            )
        }

    // run kb ref
    kb_ref(kbref_channel)

    // modify the t2g.txt file for all genome-toga pairs
    fix_t2g_names(kb_ref.out.map{ x -> record(t2g_in: x.t2g, key: x.key) })

    // part 2: kb-count
    // create an input channel
    kbcount_channel = input_channel
        .map {
            x -> record(
                species: x.species,
                sample: x.sample,
                key: "${file(x.genome).name}.${file(x.toga_dir).name}",
                read_dir: x.read_dir,
                reads: x.reads,
            )
        }
        .join(kb_ref.out, by: 'key') // do formats match here?
        .join(fix_t2g_names.out, by: 'key')

    // run kb-count
    kb_count(
        kbcount_channel,
        params.strand,
        params.pairwise
    ).view()

    // post-count table fix for each kb-count run
    post_count(kb_count.out)

    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }

    publish:
    out_dir = post_count.out.map{ x -> x.out_dir }//kb_count.out.map{ x -> x.out_dir }
}

output {
    out_dir {
        path { params.out_dir.startsWith("/") ? "${params.out_dir}" : "../${params.out_dir}" }
        mode "copy"
    }
}
