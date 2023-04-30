
process STAR_1PASS {
  tag "${sample}"
  cpus 8
  memory '35 GB'
  publishDir "${params.outdir}/${sample}/star", mode: 'copy'
  conda (params.enable_conda ? params.conda_env ?: "bioconda::star=2.7.6a" : null)

  input :
  path star
  path gtf 
  tuple val(sample), path(reads)


  output:
  tuple val(sample), path('*.out.bam'), emit: alignments
  path('*SJ.out.tab'), emit: sj
  path '*Log.final.out', emit: mqc
  path '*.{out,out.tab}'
   
 
  script:

  // additional_params = params.star_params ?: "" 

  """
  STAR \\
  --runMode alignReads \\
  --runThreadN ${task.cpus} \\
  --genomeDir ${star} \\
  --sjdbGTFfile ${gtf} \\
  --readFilesCommand zcat \\
  --readFilesIn ${reads} \\
  --genomeLoad NoSharedMemory \\
  --outSAMtype BAM SortedByCoordinate \\
  --outFileNamePrefix \\
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \\
  --outFilterScoreMinOverLread 0.4 \\
  --outFilterMatchNminOverLread 0.7 \\
  --outFilterMatchNmin 0 \\
  --outFilterMismatchNmax 2 \\
  --winAnchorMultimapNmax 10 \\
  --outFilterMultimapNmax 200 \\
  --outSAMstrandField intronMotif \\
  --outSAMattributes nM NM MD jM jI XS MC ch \\

  """

}
