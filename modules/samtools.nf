 process SAMTOOLS{
    /*
    #PBS -N BAM_filtering2properly_paired
    #PBS -q ingmq
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=20gb
    #PBS -d /mnt/projects/labs/GEBI/Lymphocytes/panepuccia/Final_DeNovo/Naive_CHR/
    #PBS -V 
    */

    tag "SAMTOOLS: PROPER PAIRING "
    cpus 4
    memory '2 GB'
    publishDir "${params.outdir}/${sample}/samtools_out", mode: 'copy'

    input:
    tuple val(sample), path(pass)

    output:
    tuple val(sample), path('*_sorted.bam') , emit: paisortprop


    script:
    """
    samtools view  -H ${pass} >> ${sample}_paired.sam;
    samtools view  ${pass}  >> ${sample}_paired.sam ;
    samtools view -47 ${pass} >> ${sample}_paired.sam;
    samtools view 83 ${pass} >> ${sample}_paired.sam;
    samtools view -f13 ${pass} >> ${sample}_paired.sam;
    samtools sort -O bam -o ${sample}_properlypaired_sorted.bam  ${sample}_paired.sam;
    
    """
    }

