process FASTQSCREEN {
    cpus 4
    memory '8 GB'

    publishDir "${params.outdir}/${sample}/fastqscreen", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path fastqscreen_conf

    output:
    path '*'
    path '*_screen.txt', emit: mqc

    script:
    """
    fastq_screen --conf ${fastqscreen_conf} ${reads}
    """
}
