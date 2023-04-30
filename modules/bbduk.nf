process BBDUK {

    publishDir "${params.outdir}/${sample}/bbduk/", mode:'copy',
                pattern: "*_rRNA_stats.txt"

    input:
    tuple val(sample), path(reads)
    path rRNAs

    output:
    tuple val(sample), path('*_BB.fastq.gz'), emit: bbduk_out  
    path "${sample}_rRNA_stats.txt", emit: mqc

    script:
    """
    bbduk.sh -Xx1g in=${reads[0]} in2=${reads[1]} \\
    out=${sample}_1_BB.fastq.gz out2=${sample}_2_BB.fastq.gz \\
    stats=${sample}_rRNA_stats.txt \\
    ref=${rRNAs} t=10 k=31 hdist=1
    """
}
