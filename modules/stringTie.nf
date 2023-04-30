
process STRINGTIE {
  tag "DeNovo"
    cpus 8
    memory '35 GB'
    publishDir "${params.outdir}/${sample}/stringtied_", mode: 'copy'

    input:
    tuple val(sample), path(paired_)
    path  gtf

    output:
    // tuple val(sample), path("*.coverage.gtf")   , emit: coverage_gtf
    //tuple val(sample), path("*.transcripts.gtf"), emit: transcript_gtf
    //tuple val(sample), path("*.abundance.txt")  , emit: abundance
    //tuple val(sample), path("*.ballgown")       , emit: ballgown
    tuple val(sample), path("*stringtie.gtf"), emit: guided_

    script: 

    """
    stringtie \\
    -i ${paired_} \\
    -o ${sample}_GUIDED_stringtie.gtf \\
    -G ${gtf} \\
    --rf\\
    """
}
