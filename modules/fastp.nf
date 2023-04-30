// vim: set filetype=groovy:
process FASTP {
    tag "${sample}"
    cpus 2
    memory '4 GB'
    publishDir "${params.outdir}/${sample}/fastp", mode: 'copy'
    conda (params.enable_conda ? 'bioconda::fastp=0.23.2' : null)
    container 'https://depot.galaxyproject.org/singularity/fastp%3A0.23.2--hb7a2d85_2'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*.trimmed.fastq.gz'), emit: trimmed
    path '*.json', emit: json
    path '*.html', emit: html
    path 'fastp.version.txt', emit: version

    script:
    def additional_params = params.fastp_params ?: ''

    if (params.single_end)
    """
     fastp --version 2>&1 | sed -e "s/fastp //g" > fastp.version.txt
     fastp \\
      --in1 ${reads} \\
      --out1 ${sample}_R1.trimmed.fastq.gz \\
      --thread ${task.cpus} \\
      ${additional_params}
    """

    else
    """
     fastp --version 2>&1 | sed -e "s/fastp //g" > fastp.version.txt
     fastp \\
      --in1 ${reads[0]} \\
      --in2 ${reads[1]} \\
      --out1 ${sample}_R1.trimmed.fastq.gz \\
      --out2 ${sample}_R2.trimmed.fastq.gz \\
      --thread ${task.cpus} \\
      ${additional_params}
    """
}
