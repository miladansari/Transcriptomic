
process MERGE_STRINGTIE {
    label 'process_medium'
    cpus 8
    memory '35 GB'
    publishDir "${params.outdir}/merge_stringtie_", mode: 'copy'
    input:
    tuple val(sample),  path (stringtie_gtf)
    path gtf

    output:
    path "stringtie.merged.gtf", emit: merged_

    script:
    if (params.merge_with)
    """
    stringtie \\
        --merge $stringtie_gtf \\
        -G $gtf \\
        -o stringtie.merged.gtf \\
        -i \\
         R \\
        -f 0 \\
    """
    else
    """
    stringtie \\
        --merge $stringtie_gtf \\
        -o stringtie.merged.gtf \\
        -i \\
        -R \\
        -f 0 \\
    """
}
