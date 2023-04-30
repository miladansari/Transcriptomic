
process CONCAT_READS_SE {
    tag "${name}"
    cpus 2
    memory '100 MB'

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("${name}.fastq.gz")

    script:

    if (reads instanceof List && reads.size() > 1)
    """
    zcat ${reads} | gzip > ${name}.fastq.gz
    """

    else
    """
    mv ${reads} ${name}.fastq.gz
    """
}

process CONCAT_READS {
    cpus 2
    memory '1 GB'
    publishDir "${params.outdir}/${sample}/concated", mode: 'copy'

    input:
    tuple val(name), path(reads), path(mates)

    output:
    tuple val(name), path("${name}_R{1,2}.fastq") , emit : concated

    script:
    if (reads instanceof List && reads.size() > 1)
    """
    zcat ${reads} | gzip > ${name}_R1.fastq
    zcat ${mates} | gzip > ${name}_R2.fastq
    """

    else
    """
    mv ${reads} ${name}_R1.fastq
    mv ${mates} ${name}_R2.fastq
    """
}

