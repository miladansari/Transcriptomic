process MULTIQC {
    cpus 2
    memory '2 GB'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    conda (params.enable_conda ? params.conda_env ?: "bioconda::multiqc=1.9" : null)

    input:
    path '*'
    path config

    output:
    path 'multiqc_*'

    script:
    """
    multiqc . --config ${config}
    """
}
