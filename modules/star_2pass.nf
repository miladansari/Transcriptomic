
process STAR_2PASS {
  tag "${sample}"
  cpus 16
  memory '40 GB'
  publishDir "${params.outdir}/${sample}/star2pass", mode: 'copy'
  conda (params.enable_conda ? params.conda_env ?: "bioconda::star=2.7.6a" : null)

  

  input :
  path star
  path gtf 
  path SJs
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path('*.out.bam') , emit: pass2
  path '*.out.bam' , emit : tanha
  path '*Log.final.out', emit: mqc
  path '*.{out,out.tab}'

  script:

  """
  STAR \\
  --runMode alignReads \\
  --runThreadN ${task.cpus} \\
  --genomeDir ${star} \\
  --sjdbGTFfile ${gtf} \\
  --sjdbFileChrStartEnd ${SJs} \\
  --readFilesCommand zcat \\
  --readFilesIn ${reads} \\
  --genomeLoad NoSharedMemory \\
  --outSAMtype BAM SortedByCoordinate \\
  --outFileNamePrefix ${sample}_ \\
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \\
  --outFilterScoreMinOverLread 0 \\
  --outFilterMatchNminOverLread 0.3 \\
   \\
  --winAnchorMultimapNmax 080 \\
  --outFilterMultimapNmax 20
  """
}
